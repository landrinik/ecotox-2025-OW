# ============================================================================
#                                  SHINY APP
# ============================================================================

# Source all modules (loaded once at startup)
source("global.R", local = TRUE)      # Libraries, data, config
source("helpers.R", local = TRUE)     # Data processing functions  
source("modeling.R", local = TRUE)    # Statistical modeling functions

# ---------------------------------------------------------------------------
# INLINE EXPORT & STYLE HELPERS (place near top of app.R after libraries)
# ---------------------------------------------------------------------------

# Central export defaults (publication quality)
plot_export_settings <- list(
  dpi   = 300,   # default DPI
  width = 8,     # inches
  height= 6,     # inches
  units = "in"   # ggsave units
)

# Simple, consistent theme for all ggplots
theme_publication <- function(
    base_size   = 14,
    base_family = ""
) {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 12, color = "#666666"),
      axis.title    = ggplot2::element_text(size = 14, face = "bold"),
      axis.text     = ggplot2::element_text(size = 11),
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title    = ggplot2::element_text(size = 12, face = "bold"),
      legend.text     = ggplot2::element_text(size = 11),
      panel.grid.major = ggplot2::element_line(color = "#E0E0E0"),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#5bc0de"),
      strip.text       = ggplot2::element_text(size = 13, face = "bold", color = "white")
    )
}

# Minimal palettes (edit here anytime)
get_fiber_palette <- function() {
  # Accept both Title and lower case keys
  vals <- c("#2E7D32", "#1565C0", "#757575")
  nm_title <- c("Cotton", "PET", "Control")
  nm_lower <- tolower(nm_title)
  pal <- setNames(rep(vals, 2), c(nm_title, nm_lower))
  pal
}

get_treatment_palette <- function() {
  vals <- c("#D32F2F", "#1976D2", "#757575")
  nm_title <- c("Treated", "Untreated", "Control")
  nm_lower <- tolower(nm_title)
  pal <- setNames(rep(vals, 2), c(nm_title, nm_lower))
  pal
}

# Central export defaults (keep your existing values)
plot_export_settings <- list(dpi = 300, width = 8, height = 6, units = "in")

# SAFER factory: adds SVG with correct MIME; only passes dpi for raster
create_plot_download_handler <- function(
    plot_reactive,
    base_filename = "plot",
    format = c("png","pdf","tiff","svg"),
    width  = plot_export_settings$width,
    height = plot_export_settings$height,
    dpi    = plot_export_settings$dpi
) {
  format <- match.arg(format)
  mime <- switch(format,
                 png  = "image/png",
                 tiff = "image/tiff",
                 pdf  = "application/pdf",
                 svg  = "image/svg+xml"
  )
  shiny::downloadHandler(
    filename = function() paste0(base_filename, "_", Sys.Date(), ".", format),
    contentType = mime,
    content  = function(file) {
      p <- plot_reactive()
      validate(need(!is.null(p), "No plot to save"))
      # SVG via svglite (no dpi)
      if (identical(format, "svg")) {
        svglite::svglite(file, width = as.numeric(width), height = as.numeric(height), bg = "white")
        on.exit(grDevices::dev.off(), add = TRUE)
        print(p)
        return(invisible())
      }
      # ggsave for the rest; omit dpi for PDF (vector), include for PNG/TIFF (raster)
      args <- list(
        filename = file,
        plot     = p,
        device   = format,
        width    = as.numeric(width),
        height   = as.numeric(height),
        units    = "in",
        bg       = "white"
      )
      if (!identical(format, "pdf")) args$dpi <- as.numeric(dpi)
      do.call(ggplot2::ggsave, args)
    }
  )
}


# Factory: table download (CSV or Excel)
create_table_download_handler <- function(
    table_reactive,
    base_filename = "table",
    format = c("csv","xlsx")
) {
  format <- match.arg(format)
  shiny::downloadHandler(
    filename = function() paste0(base_filename, "_", Sys.Date(), ".", format),
    content  = function(file) {
      df <- table_reactive()
      validate(need(!is.null(df) && nrow(df) > 0, "No data to save"))
      if (format == "csv") {
        readr::write_csv(df, file)
      } else {
        writexl::write_xlsx(df, file)
      }
    }
  )
}

# ----- Defensive emmeans/contrasts helpers (add once in server) -----

safe_emmeans <- function(model, specs, at = NULL, by = NULL) {
  tryCatch(
    emmeans::emmeans(model, specs = specs, at = at, by = by),
    error = function(e) {
      message("[safe_emmeans] ", conditionMessage(e))
      NULL
    }
  )
}

safe_contrasts <- function(emm, method, ..., adjust = "BH") {
  if (is.null(emm)) return(NULL)
  tryCatch(
    emmeans::contrast(emm, method = method, ..., adjust = adjust),
    error = function(e) {
      message("[safe_contrasts] ", conditionMessage(e))
      NULL
    }
  )
}

# --- Helpers: safe guards for observers (no silent aborts) ---

# --- Add once near DiD helpers ---
coalesce_first_existing <- function(df, cols, default = NA_character_) {
  present <- cols[cols %in% names(df)]
  if (length(present) == 0L) return(rep(default, nrow(df)))
  out <- as.character(df[[present[1]]])
  if (length(present) > 1L) {
    for (nm in present[-1]) out <- dplyr::coalesce(out, as.character(df[[nm]]))
  }
  out
}

# Helper: flag numeric-like labels (e.g., "1","02","15")
is_numeric_like <- function(x) {
  z <- trimws(as.character(x))
  nzchar(z) & grepl("^[0-9]+$", z)
}


# Helper to detect presence of any non-empty labels in a column
has_any_nonempty <- function(x) {
  if (is.null(x)) return(FALSE)
  any(!(is.na(x) | trimws(x) == ""))
}

safe_abort <- function(msg, notify = TRUE) {
  if (isTRUE(notify)) try(shiny::showNotification(msg, type = "error", duration = 6), silent = TRUE)
  message("[DiD] ABORT: ", msg)
  did_state(list(error = msg))
  return(invisible(NULL))
}

# Local utilities
did_canon_fiber <- function(x) {
  xl <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    xl %in% c("cot", "cotton") ~ "cotton",
    xl %in% c("pet", "polyester", "poly") ~ "pet",
    TRUE ~ xl
  )
}

did_canon_treatment <- function(x) {
  xl <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    xl %in% c("untreated", "control", "ctrl") ~ "untreated",
    xl %in% c("treated") ~ "treated",
    TRUE ~ NA_character_
  )
}

did_canon_conc_labels <- function(x) {
  v <- suppressWarnings(as.numeric(as.character(x)))
  allowed <- c(0, 100, 1000, 10000)
  snapped <- vapply(v, function(y) {
    if (!is.finite(y)) return(NA_real_)
    idx <- which.min(abs(allowed - y))
    tol <- max(1, 0.01 * allowed[idx])
    if (abs(allowed[idx] - y) <= tol) allowed[idx] else NA_real_
  }, numeric(1))
  as.character(as.integer(snapped))
}

nlev2 <- function(x) length(unique(stats::na.omit(as.character(x))))

# Scalar-safe coercion: always returns length-1 integer or NA
safe_int1 <- function(x) {
  v <- suppressWarnings(as.integer(x))
  if (length(v) == 0L) return(NA_integer_)
  v[1]
}

recovery_state <- shiny::reactiveVal(NULL)

# ============================================================================
#                                USER INTERFACE
# ============================================================================

ui <- fluidPage(
  useShinyjs(),
  
  # Enable Enter key to trigger Run buttons based on active tab
  tags$script(HTML("
    $(document).on('keypress', function(e) {
      if(e.which == 13) {  // Enter key pressed
        var activeTab = $('.nav-tabs .active').text().trim();
        if(activeTab === 'Regression Analysis') {
          $('#run_regression').click();
        } else if(activeTab === 'Mixed Effects Analysis') {
          $('#run_lmer').click();
        } else if(activeTab === 'DiD Analysis') {
          $('#run_did').click();
        }
      }
    });
  ")),
  
  titlePanel("The P-Value Prophet"),
  actionButton("toggleSidebar", "Toggle Sidebar"),
  
  fluidRow(
    column(
      width = 3, id = "sidebar",
      h4("Data Mode:"),
      radioButtons("active_dataset", NULL, 
                   choices = c("Assay Data" = "assay", "Physical Endpoint Data" = "physical"), 
                   selected = "assay"),
      
      # Assay controls
      conditionalPanel(
        condition = "input.active_dataset == 'assay'",
        h4("Processing mode:"),
        radioButtons("mode", NULL,
                     choices = c("Baseline (filter only)" = "baseline",
                                 "Well-level CV (wrangle_v3_avg)" = "wrangle", 
                                 "Inter-assay CV (inter_assay_CV)" = "inter"),
                     selected = "baseline"),
        br(),
        selectizeInput("assay", "Assay Type:", choices = endpoint_choices_assay, multiple = TRUE),
        selectizeInput("fiber", "Fiber Type:", choices = c("cotton", "pet"), multiple = TRUE),
        selectizeInput("sample_type", "Sample Type:",
                       choices = sample_types_assay,  # From global.R: gills, gland, hemolymph
                       multiple = TRUE),
        selectizeInput("week", "Week:", choices = week_choices_assay, multiple = TRUE),
        selectizeInput("treatment", "Treatment:", choices = c("Treated", "Untreated", "Control"), multiple = TRUE),
        selectizeInput("fiber_concentration", "Fiber Concentration:", choices = c("0", "100", "1000", "10000"), multiple = TRUE),
        conditionalPanel(condition = "input.mode != 'inter'",
                         selectizeInput("plate_replicate", "Plate Replicate:", choices = NULL, multiple = TRUE),
                         sliderInput("cv_threshold", "CV Threshold (%)", min = 0, max = 100, value = 15)),
        selectizeInput("x_axis_assay", "X-axis:", 
                       choices = c("Fiber Group" = "fiber_group", "Week" = "week", "Sample Type" = "sample_type", 
                                   "Treatment" = "treatment", "Fiber Concentration" = "fiber_concentration"), 
                       selected = "fiber_group"),
        actionButton("run_analysis", "Run Analysis"),
        br(), downloadButton("download_data", "Download Results")
      ),
      
      # Physical controls  
      conditionalPanel(
        condition = "input.active_dataset == 'physical'",
        radioButtons("mode_phys", "Mode:", choices = c("Baseline" = "baseline", "Inter CV" = "inter"), selected = "baseline"),
        selectizeInput("fiber_type_phys", "Fiber Type", choices = c("cotton", "pet"), multiple = TRUE),
        selectizeInput("week_phys", "Week", choices = sort(unique(physical_master$week)), multiple = TRUE),
        selectizeInput("endpoint", "Endpoint", choices = endpoint_choices_physical, multiple = TRUE),
        uiOutput("tissue_type_ui"),
        selectizeInput("fiber_concentration_phys", "Fiber Concentration", choices = c("0","100","1000","10000"), multiple = TRUE),
        selectizeInput("treatment_phys", "Treatment", choices = c("Control","Untreated","Treated"), multiple = TRUE),
        selectizeInput("x_axis_phys", "X Axis", 
                       choices = c("Fiber Group"="fiber_group", "Week"="week", "Fiber Type"="fiber_type", 
                                   "Endpoint"="endpoint", "Tissue Type"="tissue_type"), selected = "week"),
        actionButton("run_phys_analysis", "Run Physical Analysis"),
        downloadButton("download_phys_data", "Download Physical Data")
      )
    ),
    
    column(9,
           tabsetPanel(
             id = "analysis_tabs",
             
             # Visualization Tab
             tabPanel("Visualization",
                      h4("Results"), DT::dataTableOutput("table"), br(),
                      tags$div(style = "width: 100%; max-width: 1200px; aspect-ratio: 4/3;",
                               plotOutput("dist_plot", width = "100%", height = "100%")), br(), 
                      tags$div(style = "width: 100%; max-width: 1200px; aspect-ratio: 4/3;",
                               plotOutput("dot_plot", width = "100%", height = "100%")), br(),
                      verbatimTextOutput("info")),
             
             # Regression Analysis Tab 
             tabPanel("Regression Analysis",
                      sidebarLayout(
                        sidebarPanel(width = 4,
                                     h4("Multiple Linear Regression"),
                                     wellPanel(
                                       h5("Model Settings"),
                                       selectInput("regression_endpoint", "Select Endpoint:", choices = NULL),
                                       selectInput("regression_dataset", "Dataset:", 
                                                   choices = c("Assay Data" = "assay", "Physical Data" = "physical"), selected = "assay"),
                                       checkboxGroupInput("weeks_include", "Include Weeks:",
                                                          choices = c("1" = "1", "3" = "3", "5" = "5"),
                                                          selected = c("1","3","5"), inline = TRUE),
                                       checkboxInput("include_three_way", "Include 3-way: fiber × treatment × week", TRUE),
                                       checkboxInput("dose_by_fiber", "Include dose × fiber_type", FALSE),
                                       checkboxInput("dose_by_treat", "Include dose × chem_treatment", FALSE),
                                       actionButton("run_regression", "Run Regression Model", icon = icon("play"))
                                     ),
                                     
                                     wellPanel(
                                       h5("Data Filters"),
                                       # Assay: sample types
                                       conditionalPanel(
                                         condition = "input.regression_dataset == 'assay'",
                                         checkboxGroupInput(
                                           inputId = "reg_sample_types",
                                           label   = "Sample types:",
                                           choices = sample_types_assay,  # or sample_types_assay from global.R
                                           selected = c("gills"),
                                           inline  = TRUE
                                         ),
                                         helpText("Uncheck to exclude a sample type from modeling and EMMeans.")
                                       ),
                                       # Physical: tissues (mf_counts only)
                                       conditionalPanel(
                                         condition = "input.regression_dataset == 'physical'",
                                         uiOutput("reg_tissues_ui"),
                                         uiOutput("reg_tissues_hint")
                                       )
                                     ),
                                    
                                     wellPanel(
                                       h5("Post-hoc Comparisons"),
                                       selectInput("emmeans_by", "Calculate EM Means by:",
                                                   choices = list("Treatment Groups" = "treatment", "Fiber Type" = "fiber_type",
                                                                  "Week" = "week", "Treatment × Week" = "treatment_week",
                                                                  "Treatment × Fiber Type" = "treatment_fiber", 
                                                                  "Fiber Type × Week" = "fiber_week",
                                                                  "Treatment × Fiber × Week" = "all_interactions"),
                                                   selected = "treatment"),
                                       selectInput("emmeans_dose", "Evaluate at dose:",
                                                   choices = c("Controls (0 mf/L)" = "0", "Low dose (100 mf/L)" = "2",
                                                               "Medium dose (1000 mf/L)" = "3", "High dose (10000 mf/L)" = "4"),
                                                   selected = "0"),
                                       selectInput("comparison_type", "Comparison Type:",
                                                   choices = c("All Pairwise" = "pairwise", 
                                                               "Treatment vs Control" = "control"), 
                                                   selected = "pairwise"),
                                       selectInput("adjustment_method", "P-value Adjustment:",
                                                   choices = c("Tukey HSD" = "tukey", "Holm" = "holm", 
                                                               "BH (FDR)" = "fdr", "Bonferroni" = "bonferroni"), selected = "tukey"),
                                       actionButton("run_emmeans", "Calculate EM Means", icon = icon("chart-line"))
                                     )),
                        mainPanel(width = 8,
                                  h5("Model Summary"), 
                                  verbatimTextOutput("regression_model_summary"),
                                  
                                  h5("ANOVA Table"), 
                                  DT::dataTableOutput("regression_anova"),
                                  
                                  h5("Estimated Marginal Means"), 
                                  plotOutput("emmeans_plot", height = "420px"),
                                  DT::dataTableOutput("emmeans_table"),
                                  
                                  h5("Pairwise Comparisons"), 
                                  DT::dataTableOutput("pairwise_table"),
                                  
                                  # ============================================================================
                                  # TIME TRENDS ANALYSIS
                                  # ============================================================================
                                  hr(),
                                  h4("Optional: Time Trends Analysis"),
                                  
                                  # Sticky style for bottom controls (optional)
                                  tags$style(HTML("
  .trend-controls-bottom {
    position: sticky; top: 8px;
    background: #fff; padding: 8px 12px; margin-bottom: 8px;
    border: 1px solid #eee; border-radius: 6px; z-index: 2;
  }
")),
                                  
                                  wellPanel(
                                    style = "background-color: #f9f9f9;",
                                    
                                    checkboxInput(
                                      "show_trends",
                                      "Analyze linear trends over time",
                                      value = FALSE
                                    ),
                                    
                                    helpText(
                                      strong("What this does:"),
                                      "Tests whether outcomes change linearly across weeks for each fiber × treatment combination.",
                                      br(),
                                      "This analysis estimates the slope of change per week and compares slopes between groups.",
                                      br(), br(),
                                      strong("When to use:"),
                                      "Only useful when 'week' is included in your model and you want to know if trends differ between groups.",
                                      style = "color: #555;"
                                    ),
                                    
                                    conditionalPanel(
                                      condition = "input.show_trends",
                                      
                                      hr(),
                                      
                                      # 1) Linear trends output appears first
                                      h5("Linear Trends Output"),
                                      verbatimTextOutput("regression_trends", placeholder = TRUE),
                                      
                                      br(),
                                      
                                      # 2) MOVED HERE: controls below the summaries for easy access while viewing plot
                                      div(
                                        class = "trend-controls-bottom",
                                        hr(),
                                        h5("Adjust Trend Filters"),
                                        fluidRow(
                                          column(
                                            6,
                                            selectInput(
                                              "trend_filter_type", "Filter by:",
                                              choices = c(
                                                "All Groups" = "none",
                                                "Fiber Type" = "fiber",
                                                "Treatment"  = "treatment",
                                                "Both"       = "both"
                                              ),
                                              selected = "none"
                                            )
                                          ),
                                          column(
                                            6,
                                            selectInput(
                                              "trend_concentration", "Concentration for trends:",
                                              choices = c(
                                                "All doses (average)" = "all",
                                                "Controls (0 mf/L)"   = "0",
                                                "100 mf/L"            = "100",
                                                "1000 mf/L"           = "1000",
                                                "10000 mf/L"          = "10000"
                                              ),
                                              selected = "all"
                                            )
                                          )
                                        ),
                                        fluidRow(
                                          column(
                                            6,
                                            conditionalPanel(
                                              condition = "input.trend_filter_type == 'fiber' || input.trend_filter_type == 'both'",
                                              selectInput(
                                                "trend_fiber_select", "Select Fiber:",
                                                choices = c("cotton", "pet"),
                                                selected = "cotton"
                                              )
                                            )
                                          ),
                                          column(
                                            6,
                                            conditionalPanel(
                                              condition = "input.trend_filter_type == 'treatment' || input.trend_filter_type == 'both'",
                                              selectInput(
                                                "trend_treatment_select", "Select Treatment:",
                                                choices = c("treated", "untreated"),
                                                selected = "untreated"
                                              )
                                            )
                                          )
                                        )
                                      ),
                                      
                                      # 3) Plot comes after the moved controls
                                      h5("Trends Visualization"),
                                      plotOutput("regression_trend_plot", height = "500px"),
                                      
                                      downloadButton("download_regression_trend_png",  "PNG",  class = "btn-sm btn-primary", icon = icon("image")),
                                      downloadButton("download_regression_trend_pdf",  "PDF",  class = "btn-sm btn-primary", icon = icon("file-pdf")),
                                      downloadButton("download_regression_trend_tiff", "TIFF", class = "btn-sm btn-primary", icon = icon("image")),
                                      downloadButton("download_regression_trend_svg",  "SVG",  class = "btn-sm btn-primary", icon = icon("file-image"))
                                    )
                                  )
                        ))
                      
                      
             ),
             
             # ---------------------------------------------------------------------------
             # Combined Treatment Analysis tab
             # ---------------------------------------------------------------------------
             tabPanel(
               "Combined Treatment Analysis",
               icon = icon("layer-group"),
               br(),
               h3("14-Level Combined Treatment Model", style = "color:#2c3e50; font-weight: bold;"),
               
               # Data selection and setup
               wellPanel(
                 h5(icon("filter"), "Data Selection", style = "margin-top:0;"),
                 fluidRow(
                   column(
                     4,
                     # Endpoint selector
                     selectInput(
                       "combined_endpoint",
                       "Select Endpoint",
                       choices = c("Loading..." = "loading")  # updated in server
                     )
                   ),
                   column(
                     4,
                     # Weeks to include
                     checkboxGroupInput(
                       "combined_weeks",
                       "Select Weeks",
                       choices = c("1" = "1", "3" = "3", "5" = "5"),
                       selected = c("1", "3", "5"),
                       inline = TRUE
                     )
                   ),
                   column(
                     4,
                     # Quick data status
                     textOutput("combined_data_info"),
                     tags$style("#combined_data_info { font-weight: bold; margin-top: 25px; font-size: 14px; }")
                   )
                 ),
                 
                 fluidRow(
                   column(
                     6,
                     conditionalPanel(
                       condition = "input.active_dataset == 'assay'",
                       checkboxGroupInput(
                         inputId = "combined_sample_types",
                         label   = "Sample types:",
                         choices = sample_types_assay,  # or sample_types_assay
                         selected = c("gills","gland","hemolymph"),
                         inline  = TRUE
                       )
                     )
                   ),
                   column(
                     6,
                     conditionalPanel(
                       condition = "input.active_dataset == 'physical'",
                       uiOutput("combined_tissues_ui"),
                       uiOutput("combined_tissues_hint")
                     )
                   )
                 ),
                 hr(),
                 h4("Model Setup"),
                 fluidRow(
                   column(
                     4,
                     # Reference level for combined factor
                     selectInput(
                       "combined_reference",
                       "Reference Level",
                       choices = c(
                         "Control Cotton" = "Control Cotton",
                         "Control PET" = "Control PET",
                         "Untreated Cotton 100 mfL" = "Untreated Cotton 100 mfL",
                         "Treated Cotton 100 mfL"   = "Treated Cotton 100 mfL"
                       ),
                       selected = "Control Cotton"
                     )
                   ),
                   column(
                     4,
                     checkboxInput("combined_include_week", "Include Week (Time) Variable", value = FALSE)
                   ),
                   column(
                     4,
                     # Only show the weeks control once; the primary one above is used
                     # Left intentionally blank to keep layout aligned
                     div()
                   )
                 ),
                 fluidRow(
                   column(
                     6,
                     actionButton(
                       "run_combined_model",
                       "Run Combined Treatment Model",
                       icon = icon("play"),
                       class = "btn-primary btn-lg",
                       style = "width:100%; margin-top:10px;"
                     )
                   ),
                   column(
                     6,
                     # Enable download only after a model has been run
                     conditionalPanel(
                       condition = "input.run_combined_model > 0",
                       downloadButton(
                         "download_combined_results",
                         "Download Results",
                         class = "btn-success btn-lg",
                         style = "width:100%; margin-top:10px;"
                       )
                     )
                   )
                 )
               ),
               
               # Model Results gated by button click to avoid circular dependency on outputs
               h4("Model Results"),
               conditionalPanel(
                 condition = "input.run_combined_model > 0",
                 tabsetPanel(
                   id = "combined_results_tabs",
                   
                   tabPanel(
                     "Model Summary",
                     
                     # Collapsible card styling (scoped to this tab)
                     tags$style(HTML("
    details.fold { margin-bottom: 10px; border: 1px solid #e9ecef; border-radius: 6px; background:#fff; }
    details.fold > summary { cursor: pointer; padding: 8px 12px; font-weight: 600; background:#f7f7f7; border-bottom: 1px solid #e9ecef; list-style: none; }
    details.fold[open] > summary { background:#eef6ff; }
    details.fold .fold-body { padding: 10px 12px; }
  ")),
                     
                     br(),
                     
                     # --- Combined Treatment Variable Levels (collapsed by default) ---
                     tags$details(
                       class = "fold",
                       open = NA,  # set to TRUE to start expanded
                       
                       tags$summary("Combined Treatment Variable Levels"),
                       
                       div(
                         class = "fold-body",
                         # unchanged output id from server
                         verbatimTextOutput("combined_treatment_levels")
                       )
                     ),
                     
                     # --- Model Coefficients (collapsed by default) ---
                     tags$details(
                       class = "fold",
                       open = NA,
                       
                       tags$summary("Model Coefficients"),
                       
                       div(
                         class = "fold-body",
                         verbatimTextOutput("combined_model_summary")
                       )
                     ),
                     
                     # --- ANOVA Table (collapsed by default) ---
                     tags$details(
                       class = "fold",
                       open = NA,
                       
                       tags$summary("ANOVA Table"),
                       
                       div(
                         class = "fold-body",
                         verbatimTextOutput("combined_anova")
                       )
                     )
                   ),
                   
                   tabPanel(
                     "Estimated Means",
                     br(),
                     h5("Estimated Marginal Means by Treatment Level"),
                     DT::dataTableOutput("combined_emm"),
                     downloadButton("download_combined_emm_table_csv",  "CSV",  class = "btn-sm", icon = icon("file-csv")),
                     downloadButton("download_combined_emm_table_xlsx", "XLSX", class = "btn-sm", icon = icon("file-excel")),
                     br(),
                     
                     # New: controls to subset what the plot shows
                     wellPanel(
                       h5("EMMeans Plot Filters"),
                       fluidRow(
                         column(
                           4,
                           selectInput(
                             "combined_emm_filter_type", "Show levels by:",
                             choices = c(
                               "All Levels"   = "all",
                               "Fiber Type"   = "fiber",
                               "Treatment"    = "treat",
                               "Concentration"= "conc"
                             ),
                             selected = "all"
                           )
                         ),
                         column(
                           4,
                           conditionalPanel(
                             condition = "input.combined_emm_filter_type == 'fiber'",
                             checkboxGroupInput(
                               "combined_emm_fiber", "Fiber types:",
                               choices = c("Cotton", "PET"),
                               selected = c("Cotton","PET"), inline = TRUE
                             )
                           ),
                           conditionalPanel(
                             condition = "input.combined_emm_filter_type == 'treat'",
                             checkboxGroupInput(
                               "combined_emm_treat", "Treatments:",
                               choices = c("Control","Untreated","Treated"),
                               selected = c("Control","Untreated","Treated"), inline = TRUE
                             )
                           ),
                           conditionalPanel(
                             condition = "input.combined_emm_filter_type == 'conc'",
                             checkboxGroupInput(
                               "combined_emm_conc", "Concentrations (mf/L):",
                               choices = c("0","100","1000","10000"),
                               selected = c("0","100","1000","10000"), inline = TRUE
                             )
                           )
                         )
                       )
                     ),
                     
                     plotOutput("combined_emm_plot", height = "600px"),
                     downloadButton("download_combined_emm_png",  "PNG",  class = "btn-sm btn-success", icon = icon("image")),
                     downloadButton("download_combined_emm_pdf",  "PDF",  class = "btn-sm btn-success", icon = icon("file-pdf")),
                     downloadButton("download_combined_emm_tiff", "TIFF", class = "btn-sm btn-success", icon = icon("image")),
                     downloadButton("download_combined_emm_svg",  "SVG",  class = "btn-sm btn-success", icon = icon("file-image"))
                   ),
                   
                   # --- Pairwise Comparisons tab (reactive filter pickers added) ---
                   tabPanel(
                     "Pairwise Comparisons",
                     br(),
                     fluidRow(
                       column(
                         4,
                         selectInput(
                           "combined_comparison_filter",
                           "Filter Comparisons",
                           choices = c(
                             "All Comparisons"      = "all",
                             "vs. Reference Only"   = "ref",
                             "Within Fiber vs Control" = "fiber_control",
                             "Within Fiber Type"    = "fiber",
                             "Within Treatment"     = "treat",
                             "Within Concentration" = "conc"
                           ),
                           selected = "fiber_control"
                         )
                       ),
                       column(
                         4,
                         numericInput(
                           "combined_alpha",
                           "Significance Level",
                           value = 0.05, min = 0.001, max = 0.1, step = 0.01
                         )
                       ),
                       column(
                         4,
                         # Context pickers that appear only when relevant
                         conditionalPanel(
                           condition = "input.combined_comparison_filter == 'fiber' || input.combined_comparison_filter == 'fiber_control'",
                           selectInput(
                             "combined_filter_fiber", "Fiber:",
                             choices = c("Cotton","PET"), selected = "Cotton"
                           )
                         ),
                         conditionalPanel(
                           condition = "input.combined_comparison_filter == 'treat'",
                           selectInput(
                             "combined_filter_treat", "Treatment:",
                             choices = c("Control","Untreated","Treated"), selected = "Untreated"
                           )
                         ),
                         conditionalPanel(
                           condition = "input.combined_comparison_filter == 'conc'",
                           selectInput(
                             "combined_filter_conc", "Concentration (mf/L):",
                             choices = c("0","100","1000","10000"), selected = "0"
                           )
                         )
                       )
                     ),
                     
                     DT::dataTableOutput("combined_pairwise"),
                     downloadButton("download_combined_pairwise_csv",  "CSV",  class = "btn-sm", icon = icon("file-csv")),
                     downloadButton("download_combined_pairwise_xlsx", "XLSX", class = "btn-sm", icon = icon("file-excel")),
                     br(),
                     plotOutput("combined_pairwise_plot", height = "700px"),
                     downloadButton("download_combined_pairwise_png",  "PNG",  class = "btn-sm btn-warning", icon = icon("image")),
                     downloadButton("download_combined_pairwise_pdf",  "PDF",  class = "btn-sm btn-warning", icon = icon("file-pdf")),
                     downloadButton("download_combined_pairwise_tiff", "TIFF", class = "btn-sm btn-warning", icon = icon("image")),
                     downloadButton("download_combined_pairwise_svg",  "SVG",  class = "btn-sm btn-warning", icon = icon("file-image"))
                   ),
                   
                   tabPanel(
                     "Diagnostics",
                     br(),
                     plotOutput("combined_diagnostics", height = "800px")
                   )
                 )
               )
             ),
             
             # Mixed Effects Analysis Tab
             tabPanel("Mixed Effects Analysis",
                      sidebarLayout(
                        sidebarPanel(width = 4,
                                     h4("Linear Mixed Effects Regression"),
                                     # --- Mixed Effects: Data & filters (new) ---
                                     wellPanel(
                                       h5("Data & filters"),
                                       fluidRow(
                                         column(6,
                                                selectInput("lmer_dataset", "Dataset:",
                                                            choices = c("Assay Data" = "assay", "Physical Data" = "physical"),
                                                            selected = "assay")
                                         ),
                                         column(6,
                                                selectInput("lmer_endpoint", "Endpoint:", choices = c("Loading..."))
                                         )
                                       ),
                                       
                                       # Assay sample types OR Physical tissues (dynamic)
                                       conditionalPanel(
                                         condition = "input.lmer_dataset == 'assay'",
                                         checkboxGroupInput("lmer_sample_types", "Sample types:",
                                                            choices = sample_types_assay, selected = sample_types_assay, inline = TRUE)
                                       ),
                                       conditionalPanel(
                                         condition = "input.lmer_dataset == 'physical'",
                                         uiOutput("lmer_tissues_ui")  # shows choices for mf_counts; disabled “All” otherwise
                                       ),
                                       
                                       # Two short rows for core filters
                                       fluidRow(
                                         column(6,
                                                checkboxGroupInput("lmer_fiber_types", "Fiber types:",
                                                                   choices = c("cotton","pet"), selected = c("cotton","pet"), inline = TRUE)
                                         ),
                                         column(6,
                                                checkboxGroupInput("lmer_treatments", "Treatments:",
                                                                   choices = c("untreated","treated"), selected = c("untreated","treated"), inline = TRUE)
                                         )
                                       ),
                                       fluidRow(
                                         column(12,
                                                checkboxGroupInput("lmer_concentrations", "Fiber concentration (mf/L):",
                                                                   choices = c("0","100","1000","10000"),
                                                                   selected = c("0","100","1000","10000"), inline = TRUE),
                                                helpText("Weeks are controlled below in Include Weeks.")
                                         )
                                       )
                                     ),
                                     
                                     # --- Mixed Effects: Random-effects structure (single section) ---
                                     wellPanel(
                                       h5("Random-effects structure"),
                                       
                                       # Presets keep the UI small; users can flip to Advanced for full control
                                       radioButtons(
                                         "lmer_re_preset", "Preset:",
                                         choices = c(
                                           "Tank + batch (nested)" = "nested_tank_batch",
                                           "Tank only"             = "tank_only",
                                           "Batch only"            = "batch_only",
                                           "Custom (use Advanced)" = "custom"
                                         ),
                                         selected = "nested_tank_batch", inline = FALSE
                                       ),
                                       
                                       checkboxInput("lmer_re_show_advanced", "Show advanced options", value = FALSE),
                                       
                                       conditionalPanel(
                                         condition = "input.lmer_re_show_advanced",
                                         checkboxGroupInput(
                                           "lmer_random_intercepts", "Random intercepts:",
                                           choices  = c("Tank" = "tank", "Experiment batch" = "experiment_batch", "Individual" = "individual_id"),
                                           selected = c("tank","experiment_batch")
                                         ),
                                         checkboxInput("lmer_nested_batch_tank", "Nest tank within experiment_batch", value = TRUE),
                                         selectInput("lmer_random_slope", "Random slope:",
                                                     choices = c("None"="none","Week within tank"="week|tank",
                                                                 "Dose within tank"="dose_log10|tank","Week within batch"="week|experiment_batch"),
                                                     selected = "none")
                                       )
                                     ),
                                     
                                     # --- Mixed Effects: Model settings (remove the duplicate random-intercepts here) ---
                                     wellPanel(
                                       h5("Model Settings"),
                                       checkboxInput("lmer_include_three_way", "Include 3‑way: fiber × treatment × week", value = TRUE),
                                       checkboxInput("lmer_dose_by_fiber", "Include dose × fiber_type", value = FALSE),
                                       checkboxInput("lmer_dose_by_treat", "Include dose × chem_treatment", value = FALSE),
                                       checkboxGroupInput("lmer_weeks_include", "Include Weeks:", choices = c("1","3","5"),
                                                          selected = c("1","3","5"), inline = TRUE),
                                       actionButton("run_lmer", "Run Mixed Model", icon = icon("play"))
                                     ),
                                     wellPanel(
                                       h5("Post-hoc Comparisons"),
                                       selectInput("lmer_emmeans_by", "Calculate EM Means by:",
                                                   choices = list("Treatment Groups" = "treatment", "Fiber Type" = "fiber_type", 
                                                                  "Week" = "week", "Treatment × Week" = "treatment_week",
                                                                  "Treatment × Fiber Type" = "treatment_fiber",
                                                                  "Fiber Type × Week" = "fiber_week", 
                                                                  "Treatment × Fiber × Week" = "all_interactions"), selected = "treatment"),
                                       selectInput("lmer_emmeans_dose", "Evaluate at dose:",
                                                   choices = c("Controls (0 mf/L)" = "0", "Low dose (100 mf/L)" = "2",
                                                               "Medium dose (1000 mf/L)" = "3", "High dose (10000 mf/L)" = "4"), selected = "0"),
                                       selectInput("lmer_comparison_type", "Comparison Type:",
                                                   choices = c("All Pairwise" = "pairwise", 
                                                               "Treatment vs Control" = "control"),
                                                   selected = "pairwise"),
                                       selectInput("lmer_adjustment_method", "P-value Adjustment:",
                                                   choices = c("Tukey HSD" = "tukey", "Holm" = "holm",
                                                               "BH (FDR)" = "fdr", "Bonferroni" = "bonferroni"), selected = "tukey"),
                                       actionButton("lmer_run_emmeans", "Calculate EM Means", icon = icon("chart-line"))
                                     )),
                        mainPanel(width = 8,
                                  h5("Mixed Model Summary"), verbatimTextOutput("lmer_model_summary"),
                                  h5("ANOVA Table"), DT::dataTableOutput("lmer_anova"), 
                                  h5("Estimated Marginal Means"), plotOutput("lmer_emmeans_plot", height = "420px"),
                                  DT::dataTableOutput("lmer_emmeans_table"),
                                  h5("Pairwise / Custom Contrasts"), DT::dataTableOutput("lmer_pairwise_table"),
                                  h5("Model Comparison (lm vs lmer)"), DT::dataTableOutput("model_compare")))
             ),
             
             # ============================================================================
             # RECOVERY & TRANSLOCATION TAB
             # ============================================================================
             tabPanel(
               "Recovery & Translocation Analysis",
               fluidRow(
                 column(
                   4,
                   wellPanel(
                     h4("Recovery Analysis (Week vs Week Comparison)"),
                     
                     # Dataset and endpoint
                     fluidRow(
                       column(
                         6,
                         selectInput(
                           "recovery_dataset", "Dataset:",
                           choices = c("Assay Data" = "assay", "Physical Data" = "physical"),
                           selected = "assay"
                         )
                       ),
                       column(
                         6,
                         selectInput("recovery_endpoint", "Endpoint:", choices = NULL)
                       )
                     ),
                     
                     # Fixed weeks
                     selectInput(
                       inputId = "recovery_baseline_week",
                       label   = "Baseline Week:",
                       choices = c(),   # server will update
                       selected = NULL
                     ),
                     
                     selectInput(
                       inputId = "recovery_recovery_week",
                       label   = "Recovery Week:",
                       choices = c(),   # server will update
                       selected = NULL
                     ),
                     h5("Recovery-Specific Filters:"),
                     
                     # Row 1: Fiber + Treatment (always visible; fiber selector restored)
                     fluidRow(
                       column(
                         6,
                         checkboxGroupInput(
                           inputId  = "recovery_fiber_types",
                           label    = "Fiber Types:",
                           choices  = c("Cotton" = "cotton", "PET" = "pet"),
                           selected = c("cotton", "pet"),
                           inline   = TRUE
                         )
                       ),
                       column(
                         6,
                         checkboxGroupInput(
                           inputId  = "recovery_treatments",
                           label    = "Treatments:",
                           choices  = c("Untreated" = "untreated", "Treated" = "treated"),
                           selected = c("untreated", "treated")
                         )
                       )
                     ),
                     
                     # Row 2: Assay-only sample types and dose
                     conditionalPanel(
                       condition = "input.recovery_dataset == 'assay'",
                       checkboxGroupInput(
                         "recovery_samples", "Sample Types:",
                         choices  = c("Hemolymph", "Gills", "Gland"),
                         selected = c("Hemolymph", "Gills", "Gland")
                       ),
                       checkboxGroupInput(
                         "recovery_concentration", "Fiber concentration (mf/L):",
                         choices  = c("0"="0","100"="100","1000"="1000","10000"="10000"),
                         selected = c("0","100","1000","10000"),
                         inline   = TRUE
                       )
                     ),
                     
                     # Physical: tissues + dose for mf_counts
                     conditionalPanel(
                       condition = "input.recovery_dataset == 'physical'",
                       uiOutput("recovery_tissues_ui"),
                       uiOutput("recovery_tissues_hint"),
                       checkboxGroupInput(
                         "recovery_concentration_physical", "Fiber concentration (mf/L):",
                         choices  = c("0"="0","100"="100","1000"="1000","10000"="10000"),
                         selected = c("0","100","1000","10000"),
                         inline   = TRUE
                       )
                     ),
                     
                     # Run controls
                     fluidRow(
                       column(6, actionButton("run_recovery", "Run Recovery Analysis", class = "btn-primary")),
                       column(6, checkboxInput("recovery_include_emmeans", "Include EM Means Analysis", value = TRUE))
                     ),
                     
                     helpText("Recovery compares Week 6 vs the selected baseline week (1 or 5) using filters independent of the main sidebar.")
                   ),
                   
                   # Translocation stub
                   wellPanel(
                     h4("Translocation Analysis"),
                     # Week selector (single)
                     selectizeInput(
                       "transloc_week", "Week:",
                       choices = sort(unique(physical_master$week)),     # assumes weeks present in physical
                       selected = max(sort(unique(physical_master$week))), 
                       multiple = FALSE
                     ),
                     # Two-tissue selector (exactly 2)
                     selectizeInput(
                       "transloc_tissues", "Tissues to compare:",
                       choices = tissue_types_physical,                   # c("gills","gland","tissue") from global.R
                       multiple = TRUE, options = list(maxItems = 2)
                     ),
                     # Fiber concentration (multi)
                     selectizeInput(
                       "transloc_conc", "Fiber concentration (mF/L):",
                       choices = c("0","100","1000","10000"),
                       selected = c("0","100","1000","10000"),
                       multiple = TRUE
                     ),
                     
                     # --- Extra controls for the new multi-week translocation plot ---
                     h5("Multi‑week plot controls"),
                     
                     # Weeks: allow multiple
                     selectizeInput(
                       inputId = "transloc_weeks_multi",
                       label   = "Weeks (multi):",
                       choices = sort(unique(physical_master$week)), # uses your loaded physical data
                       selected = sort(unique(physical_master$week)),  # preselect all available
                       multiple = TRUE
                     ),
                     
                     # Tissue: enforce a single tissue
                     selectizeInput(
                       inputId = "transloc_tissue_single",
                       label   = "Tissue (single):",
                       choices = tissue_types_physical,  # from global.R
                       selected = head(tissue_types_physical, 1),
                       multiple = FALSE
                     ),
                     
                     # Concentrations: default to 100/1000/10000
                     checkboxGroupInput(
                       inputId  = "transloc_conc_multi",
                       label    = "Concentrations (mf/L):",
                       choices  = c("100","1000","10000"),
                       selected = c("100","1000","10000"),
                       inline   = TRUE
                     ),
                     actionButton("run_translocation", "Run Translocation Analysis")
                   )
                 ),
                 
                 # Results
                 column(
                   8,
                   tabsetPanel(
                     tabPanel(
                       "Recovery Model",
                       h4("Model Summary"),
                       verbatimTextOutput("recovery_model_summary"),
                       br(),
                       h4("Pairwise Comparisons"),
                       DT::dataTableOutput("recovery_comparisons")
                     ),
                     tabPanel("Recovery Plot", plotOutput("recovery_plot", height = "500px")),
                     tabPanel(
                       "Recovery EM Means",
                       conditionalPanel(
                         condition = "input.recovery_include_emmeans",
                         
                         h4("Estimated Marginal Means"),
                         DT::dataTableOutput("recovery_emmeans_table"),
                         br(),
                         fluidRow(
                           column(
                             width = 4,
                             selectInput(
                               inputId = "recovery_p_adjust",
                               label   = "P-value adjust:",
                               choices  = c("BH (FDR)" = "BH",
                                            "Holm"     = "holm",
                                            "Bonferroni" = "bonferroni",
                                            "None"     = "none"),
                               selected = "BH"
                             )
                           )
                         ),
                         
                         h4("p-values"),
                         DT::dataTableOutput("recovery_significance_table"),
                         
                         br(),
                         h4("EM Means Plot"),
                         plotOutput("recovery_emmeans_plot", height = "500px")
                       ),
                       conditionalPanel(
                         condition = "!input.recovery_include_emmeans",
                         div(class = "alert alert-info", "Enable 'Include EM Means Analysis' to view emmeans results.")
                       )
                     ),
                     tabPanel("Translocation Summary", DT::dataTableOutput("transloc_summary_table")),
                     tabPanel("Translocation Plots",
                              plotOutput("transloc_counts_plot", height = 420),
                              plotOutput("transloc_diff_plot",   height = 360),
                              plotOutput("transloc_multiweek_plot", height = 420)
                              )
                     
                   )
                 )
               )
             ),
             
             # ============================================================================
             #                                  DID TAB 
             # ============================================================================
             tabPanel(
               "Difference-in-Differences",
               fluidRow(
                 column(
                   4,
                   wellPanel(
                     h4("DiD Setup"),
                     
                     # Dataset + endpoint (unchanged)
                     fluidRow(
                       column(6,
                              selectInput(
                                "did_dataset", "Dataset:",
                                choices = c("Assay Data" = "assay", "Physical Data" = "physical"),
                                selected = "assay"
                              )
                       ),
                       column(6,
                              selectInput("did_endpoint", "Endpoint:", choices = NULL)
                       )
                     ),
                     
                     # DiD mode selection (unchanged semantics, label text harmonized)
                     radioButtons(
                       "did_mode", "DiD Mode:",
                       choices = c(
                         "Across Fiber (at week)" = "across_fiber",
                         "Across Week (within fiber)" = "across_week",
                         "Across Treatment (within week)" = "across_treatment"
                       ),
                       selected = "across_fiber"
                     ),
                     
                     # Sample/Tissue selection - only for across_fiber and across_week
                     conditionalPanel(
                       condition = "input.did_mode == 'across_fiber' || input.did_mode == 'across_week'",
                       conditionalPanel(
                         condition = "input.did_dataset == 'assay'",
                         div(
                           id = "did_samples_box",
                           checkboxGroupInput(
                             inputId = "did_samples",
                             label   = "Sample types:",
                             choices = character(0),     # filled by server
                             selected = NULL,
                             inline  = TRUE
                           ),
                           uiOutput("did_samples_hint")  # hint text managed by server
                         )
                       ),
                       
                       # Physical selectors
                       conditionalPanel(
                         condition = "input.did_dataset == 'physical'",
                         div(
                           id = "did_tissues_box",
                           checkboxGroupInput(
                             inputId = "did_tissues",
                             label   = "Tissues:",
                             choices = character(0),     # filled by server
                             selected = NULL,
                             inline  = TRUE
                           ),
                           uiOutput("did_tissues_hint")  # hint text managed by server
                         )
                       )
                     ),
                     
                     # Week selection - for across_fiber and across_treatment
                     conditionalPanel(
                       condition = "input.did_mode == 'across_fiber' || input.did_mode == 'across_treatment'",
                       selectInput("did_week_select", "Select Week:",
                                   choices = c("1","3","5"), selected = "5")
                     ),
                     
                     # Fiber selection - for across_week and across_treatment
                     conditionalPanel(
                       condition = "input.did_mode == 'across_week' || input.did_mode == 'across_treatment'",
                       selectInput("did_fiber_select", "Select Fiber:",
                                   choices = c("Cotton" = "cotton", "PET" = "pet"),
                                   selected = "cotton")
                     ),
                     
                     # Contrast type - only for across_fiber
                     conditionalPanel(
                       condition = "input.did_mode == 'across_fiber'",
                       radioButtons(
                         "did_across_fiber_contrast", "Contrast Type:",
                         choices = c(
                           "DiD (fiber effect on treatment)" = "did",
                           "Single-treatment fiber contrast" = "single_trt"
                         ),
                         selected = "did"
                       )
                     ),
                     
                     # Single treatment selector - only when single_trt is selected
                     conditionalPanel(
                       condition = "input.did_mode == 'across_fiber' && input.did_across_fiber_contrast == 'single_trt'",
                       selectInput(
                         "did_single_treatment", "Select Treatment:",
                         choices = c("Treated" = "treated", "Untreated" = "untreated"),
                         selected = "treated"
                       )
                     ),
                     
                     # Fiber A and B - only for across_fiber DiD mode
                     conditionalPanel(
                       condition = "input.did_mode == 'across_fiber' && input.did_across_fiber_contrast == 'did'",
                       selectInput("did_fiber_a", "Fiber A:",
                                   choices = c("Cotton" = "cotton", "PET" = "pet"),
                                   selected = "pet"),
                       selectInput("did_fiber_b", "Fiber B:",
                                   choices = c("Cotton" = "cotton", "PET" = "pet"),
                                   selected = "cotton")
                     ),
                     
                     # Week A and B - only for across_week
                     conditionalPanel(
                       condition = "input.did_mode == 'across_week'",
                       selectInput("did_week_a", "Week A:",
                                   choices = c("1", "3", "5"), selected = "3"),
                       selectInput("did_week_b", "Week B:",
                                   choices = c("1", "3", "5"), selected = "5")
                     ),
                     
                     # Baseline selector (now includes untreated_ctrl)
                     conditionalPanel(
                       condition = "input.did_mode == 'across_fiber'",
                       radioButtons(
                         inputId = "did_baseline_type",
                         label   = "Baseline for DiD:",
                         choices = c(
                           "Untreated vs Treated (current)"        = "treatment",
                           "Dose vs 0 control (new)"               = "dose",
                           "Untreated dose vs 0 (control-offset)"  = "untreated_ctrl"   # NEW
                         ),
                         selected = "treatment",
                         inline = FALSE
                       ),
                       # Dose picker shown for both 'dose' and 'untreated_ctrl'
                       conditionalPanel(
                         condition = "input.did_baseline_type == 'dose' || input.did_baseline_type == 'untreated_ctrl'",
                         selectInput(
                           inputId = "did_treated_dose",
                           label   = "Dose for baseline (mf/L):",
                           choices = c("100","1000","10000"),
                           selected = "100"
                         )
                       )
                     ),
                     
                     # Treatment A and B - only for across_treatment
                     conditionalPanel(
                       condition = "input.did_mode == 'across_treatment'",
                       selectInput("did_trt_a", "Treatment A:",
                                   choices = c("Treated" = "treated", "Untreated" = "untreated"),
                                   selected = "treated"),
                       selectInput("did_trt_b", "Treatment B:",
                                   choices = c("Treated" = "treated", "Untreated" = "untreated"),
                                   selected = "untreated")
                     ),
                     
                     # Run button
                     hr(),
                     actionButton("run_did", "Run DiD Analysis", class = "btn-primary")
                   )
                 ),
                 
                 column(
                   8,
                   tabsetPanel(
                     tabPanel(
                       "DiD Result",
                       h4("Custom DiD Contrast"),
                       verbatimTextOutput("did_summary"),
                       br(),
                       DT::dataTableOutput("did_table")
                     ),
                     tabPanel(
                       "EMMeans Grid",
                       DT::dataTableOutput("did_emm_table")
                     ),
                     tabPanel(
                       "Visualization",
                       plotOutput("did_plot", height = "480px")
                     )
                   )
                 )
               )
             )
           )
    )
  )
)

# ============================================================================
#                             SERVER FUNCTION 
# ============================================================================

server <- function(input, output, session) {
  
  # Toggle sidebar
  observeEvent(input$toggleSidebar, { toggle("sidebar") })
  
  # Update plate replicate choices from pre-loaded data
  observe({
    updateSelectizeInput(session, "plate_replicate", 
                         choices = sort(unique(final_data$plate_replicate)))
  })
  
  # Dynamic tissue filter UI for mf_counts endpoint
  output$tissue_type_ui <- renderUI({
    # Only render in the Physical Endpoint Data mode
    req(input$active_dataset == "physical")  # FIX: correct id
    
    # Only show when mf_counts is among the selected endpoints
    if (!is.null(input$endpoint) && "mf_counts" %in% input$endpoint) {
      selectizeInput(
        inputId = "tissue_type_phys",
        label   = "Tissue Type",
        choices = tissue_types_physical,   # from global.R
        multiple = TRUE,
        selected = NULL
      )
    } else {
      NULL  # Hide the input for non-mf_counts endpoints
    }
  })
  
  # ---------------------------
  # Regression: tissues UI (physical)
  # ---------------------------
  output$reg_tissues_ui <- renderUI({
    req(input$regression_dataset == "physical")
    ep <- input$regression_endpoint %||% ""
    if (identical(ep, "mf_counts")) {
      checkboxGroupInput(
        inputId = "reg_tissues",
        label   = "Tissues:",
        choices = tissue_types_physical,   # from global.R, e.g., c("gills","gland","tissue")
        selected = tissue_types_physical,
        inline   = TRUE
      )
    } else {
      # Show a disabled “All” when tissues don’t apply
      checkboxGroupInput("reg_tissues", "Tissues:", choices = c("All"), selected = "All", inline = TRUE) %>%
        shinyjs::disable()
    }
  })
  output$reg_tissues_hint <- renderUI({
    req(input$regression_dataset == "physical")
    ep <- input$regression_endpoint %||% ""
    if (!identical(ep, "mf_counts")) {
      div(class = "text-muted small", "Tissue filters apply only to mf_counts; other physical endpoints have no tissue labels.")
    }
  })
  
  # ---------------------------
  # Combined: tissues UI (physical)
  # ---------------------------
  output$combined_tissues_ui <- renderUI({
    req(input$active_dataset == "physical")
    ep <- input$combined_endpoint %||% ""
    if (identical(ep, "mf_counts")) {
      checkboxGroupInput(
        inputId = "combined_tissues",
        label   = "Tissues:",
        choices = tissue_types_physical,
        selected = tissue_types_physical,
        inline   = TRUE
      )
    } else {
      checkboxGroupInput("combined_tissues", "Tissues:", choices = c("All"), selected = "All", inline = TRUE) %>%
        shinyjs::disable()
    }
  })
  output$combined_tissues_hint <- renderUI({
    req(input$active_dataset == "physical")
    ep <- input$combined_endpoint %||% ""
    if (!identical(ep, "mf_counts")) {
      div(class = "text-muted small", "Tissue filters apply only to mf_counts; other physical endpoints have no tissue labels.")
    }
  })
  
  # ============================================================================
  # Dynamic sample type choices for assay data
  # ============================================================================
  observe({
    # Only run for assay data
    req(input$active_dataset == "assay")
    
    # Get available sample types from actual data
    if (exists("final_data") && nrow(final_data) > 0) {
      actual_samples <- sort(unique(tolower(final_data$sample_type)))
      
      # Filter to allowed sample types only
      available_samples <- intersect(sample_types_assay, actual_samples)
      if (length(available_samples) == 0) {
        available_samples <- sample_types_assay  # Fallback to all allowed types
      }
      
      # Update the choices dynamically
      updateSelectizeInput(session, "sample_type", 
                           choices = available_samples,
                           selected = NULL)  # Select all by default
    }
  })
  
  # Apply presets to advanced widgets; users can still tweak after toggling Advanced
  observeEvent(input$lmer_re_preset, {
    preset <- input$lmer_re_preset %||% "nested_tank_batch"
    
    if (preset == "nested_tank_batch") {
      updateCheckboxGroupInput(session, "lmer_random_intercepts",
                               selected = c("tank","experiment_batch"))
      updateCheckboxInput(session, "lmer_nested_batch_tank", value = TRUE)
      updateSelectInput(session, "lmer_random_slope", selected = "none")
    } else if (preset == "tank_only") {
      updateCheckboxGroupInput(session, "lmer_random_intercepts",
                               selected = c("tank"))
      updateCheckboxInput(session, "lmer_nested_batch_tank", value = FALSE)
      updateSelectInput(session, "lmer_random_slope", selected = "none")
    } else if (preset == "batch_only") {
      updateCheckboxGroupInput(session, "lmer_random_intercepts",
                               selected = c("experiment_batch"))
      updateCheckboxInput(session, "lmer_nested_batch_tank", value = FALSE)
      updateSelectInput(session, "lmer_random_slope", selected = "none")
    }
    # Custom: leave existing selections unchanged
  })
  
  # ---- Helpers for normalization ----
  normalize_fiber <- function(x) {
    x_low <- tolower(x)
    dplyr::case_when(
      x_low %in% c("cotton", "cot") ~ "cotton",
      x_low %in% c("pet", "polyester") ~ "pet",
      TRUE ~ x_low
    )
  }
  
  to_chem_treatment <- function(x) {
    x_low <- tolower(x)
    dplyr::case_when(
      x_low %in% c("untreated", "control") ~ "untreated",
      x_low %in% c("treated") ~ "treated",
      TRUE ~ NA_character_
    )
  }
  
  
  # ============================================================================
  # DATA FILTERING - FIXED: Use reactive() with bindCache(), then bindEvent()
  # ============================================================================
  
  assay_filtered_cached <- reactive({
    df <- final_data  # Pre-processed data
    
    if (length(input$assay)) {
      df <- df %>% filter(assay_type %in% input$assay)
      message("After assay filter: ", nrow(df), " rows")
    }
    
    # FIXED: Case-insensitive fiber type filtering
    if (length(input$fiber)) {
      df <- df %>% filter(tolower(fiber_type) %in% tolower(input$fiber))
      message("After fiber filter: ", nrow(df), " rows")
    }
    
    # FIXED: Case-insensitive sample type filtering
    if (length(input$sample_type)) {
      df <- df %>% filter(tolower(sample_type) %in% tolower(input$sample_type))
      message("After sample_type filter: ", nrow(df), " rows")
    }
    
    if (length(input$week)) {
      df <- df %>% filter(week %in% as.integer(input$week))
      message("After week filter: ", nrow(df), " rows")
    }
    
    if (length(input$fiber_concentration)) {
      df <- df %>% filter(fiber_concentration %in% input$fiber_concentration)
      message("After fiber_concentration filter: ", nrow(df), " rows")
    }
    
    if (length(input$treatment)) {
      if ("Control" %in% input$treatment) {
        df <- df %>% filter(fiber_concentration == "0" | treatment %in% setdiff(input$treatment, "Control"))
      } else {
        df <- df %>% filter(treatment %in% input$treatment)
      }
      message("After treatment filter: ", nrow(df), " rows")
    }
    
    if (input$mode != "inter" && length(input$plate_replicate)) {
      df <- df %>% filter(plate_replicate %in% input$plate_replicate)
      message("After plate_replicate filter: ", nrow(df), " rows")
    }
    
    message("Final filtered data rows: ", nrow(df))
    df
  }) %>% bindCache(input$assay, input$fiber, input$sample_type, input$week, 
                   input$fiber_concentration, input$treatment, input$plate_replicate, input$mode)
  
  assay_filtered <- assay_filtered_cached %>% bindEvent(input$run_analysis)
  
  #Physical data filter - Use reactive with proper caching  
  phys_filtered_cached <- reactive({
    df <- physical_master  # Pre-processed data
    
    if (length(input$fiber_type_phys)) {
      df <- df %>% filter(fiber_type %in% input$fiber_type_phys)
    }
    if (length(input$week_phys)) {
      df <- df %>% filter(week %in% input$week_phys)
    }
    if (length(input$endpoint)) {
      df <- df %>% filter(endpoint %in% input$endpoint)
    }
    
    # Apply tissue_type filter for mf_counts only
    if (length(input$tissue_type_phys) > 0) {
      if ("mf_counts" %in% df$endpoint) {
        # For mf_counts: filter by tissue_type
        # For other endpoints: keep all (tissue_type is NA for them)
        df <- df %>%
          filter(
            (endpoint == "mf_counts" & tissue_type %in% input$tissue_type_phys) |
              (endpoint != "mf_counts")
          )
      }
    }
    if (length(input$fiber_concentration_phys)) df <- df %>% filter(fiber_concentration %in% input$fiber_concentration_phys)
    if (length(input$treatment_phys)) df <- df %>% filter(treatment %in% input$treatment_phys)
    df
  }) %>% bindCache(input$fiber_type_phys, input$week_phys, input$endpoint, 
                   input$tissue_type_phys, input$fiber_concentration_phys, input$treatment_phys)
  
  # Now use bindEvent
  phys_filtered <- phys_filtered_cached %>% bindEvent(input$run_phys_analysis)
  
  # ============================================================================  
  # DATA SUMMARIZATION - Using helper functions from helpers.R
  # ============================================================================
  
  summarized <- reactive({
    req(input$active_dataset == "assay")
    df <- assay_filtered()
    if (nrow(df) == 0) return(NULL)
    
    threshold <- input$cv_threshold
    if (input$mode == "baseline") {
      summarize_baseline(df, threshold)
    } else if (input$mode == "wrangle") {
      summarize_wrangle(df, threshold)  
    } else {
      summarize_inter(df, threshold)
    }
  })
  
  phys_summarized <- reactive({
    df <- phys_filtered()
    req(nrow(df) > 0)
    
    # Control correction is already applied in convert_data_debugged.R
    # Just select relevant columns and clean up
    display_cols <- c("fiber_type", "week", "endpoint", "tank", "sample", 
                      "treatment", "fiber_concentration", "fiber_group", "value")
    
    # Add tissue_type for mf_counts
    if ("tissue_type" %in% names(df) && "mf_counts" %in% df$endpoint) {
      display_cols <- c(display_cols, "tissue_type")
    }
    
    # Add value_corr for mf_counts (control-corrected values)
    if ("value_corr" %in% names(df) && "mf_counts" %in% df$endpoint) {
      display_cols <- c(display_cols, "value_corr")
    }
    
    # Select only relevant columns
    df <- df %>%
      dplyr::select(dplyr::any_of(display_cols)) %>%
      mutate(
        value = if_else(endpoint == "mf_counts", round(value, 1), value)
      )
    
    if (input$mode_phys == "baseline") {
      # Baseline mode: use raw values
      df
    } else {
      # Summary modes: use corrected values for mf_counts
      df <- df %>%
        mutate(
          use_val = if_else(endpoint == "mf_counts" & "value_corr" %in% names(.), value_corr, value)
        ) %>%
        group_by(fiber_group, week, endpoint, tissue_type, fiber_concentration) %>%
        summarise(
          mean_value = mean(use_val, na.rm = TRUE),
          sd_value = sd(use_val, na.rm = TRUE),
          n = n(),
          cv = if_else(mean_value != 0, 100 * sd_value / mean_value, NA_real_),
          .groups = "drop"
        ) %>%
        mutate(
          mean_value = if_else(abs(mean_value) < 1, round(mean_value, 4), round(mean_value, 2)),
          sd_value = signif(sd_value, 5),
          cv = round(cv, 2)
        )
      
      df
    }
  })
  
  # ============================================================================
  # OUTPUTS - Optimized rendering  
  # ============================================================================
  
  active_table <- reactive({
    if (input$active_dataset == "assay") summarized() else phys_summarized()
  })
  
  output$table <- DT::renderDataTable({
    df <- active_table()
    if (is.null(df) || nrow(df) == 0) {
      return(datatable(data.frame(Message = "No data for selected filters / mode")))
    }
    
    if("endpoint" %in% names(df) && "value" %in% names(df)) {
      df <- df %>% mutate(value = if_else(endpoint == "mf_counts", round(value), value))
    }
    
    datatable(df, filter = "top", options = list(pageLength = 15))
  })
  
  # Distribution plot - optimized
  output$dist_plot <- renderPlot({
    df <- active_table()
    req(df, nrow(df) > 0)
    
    if (input$active_dataset == "assay") {
      x_axis_var <- input$x_axis_assay %||% "fiber_group"
      if (x_axis_var == "week") df$week <- factor(df$week)
      
      ggplot(df, aes_string(x = x_axis_var, y = "mean_activity", fill = "fiber_concentration")) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
        facet_wrap(~ assay_type) + theme_minimal() +
        labs(title = glue("Assay Activity by {str_to_title(gsub('_', ' ', x_axis_var))}"),
             x = str_to_title(gsub("_", " ", x_axis_var)), y = "Mean Activity", fill = "Fiber Concentration") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } else {
      x_axis_var <- input$x_axis_phys %||% "week"
      x_axis_factor <- paste0("factor(", x_axis_var, ")")
      
      if (input$mode_phys == "baseline") {
        ggplot(df, aes_string(x = x_axis_factor, y = "value", fill = "fiber_concentration")) +
          geom_boxplot(outlier.shape = NA, alpha = 0.7) + geom_jitter(width = 0.18, alpha = 0.6, size = 2) +
          facet_wrap(~ endpoint, scales = "free_y") + theme_minimal() +
          labs(title = glue("Physical Baseline by {str_to_title(gsub('_', ' ', x_axis_var))}"),
               x = str_to_title(gsub("_", " ", x_axis_var)), y = "Value", fill = "Fiber Concentration") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      } else {
        ggplot(df, aes_string(x = x_axis_factor, y = "mean_value", fill = "fiber_concentration")) +
          geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.7) +
          geom_point(aes_string(group = "fiber_concentration"), position = position_dodge(width = 0.8), 
                     size = 2, color = "black", show.legend = FALSE) +
          facet_wrap(~ endpoint, scales = "free_y") + theme_minimal() +
          labs(title = glue("Physical Inter CV by {str_to_title(gsub('_', ' ', x_axis_var))}"),
               x = str_to_title(gsub("_", " ", x_axis_var)), y = "Mean Value", fill = "Fiber Concentration") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
    }
  })
  
  # Dot plot - assay only
  output$dot_plot <- renderPlot({
    req(input$active_dataset == "assay")
    df <- assay_filtered()
    req(!is.null(df), nrow(df) > 0)
    
    required_cols <- c("fiber_group", "calculated_concentration", "fiber_concentration", "assay_type")
    if (!all(required_cols %in% colnames(df))) return(NULL)
    
    ggplot(df, aes(x = fiber_group, y = calculated_concentration, color = fiber_concentration)) +
      geom_jitter(width = 0.3, alpha = 0.7, size = 2) + facet_wrap(~ assay_type) + theme_minimal() +
      labs(title = "Sample-Level Distribution by Group and Fiber Concentration",
           x = "Fiber Group", y = "Calculated Concentration", color = "Fiber Concentration") +
      theme(plot.title = element_text(size = 14, face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Info text
  output$info <- renderText({
    df <- active_table()
    if (is.null(df)) return("No data for current filters / mode.")
    
    if (input$active_dataset == "assay" && "cv_flag" %in% colnames(df)) {
      flagged <- sum(coalesce(df$cv_flag, FALSE), na.rm = TRUE)
      total <- nrow(df)
      glue("Assay mode: {input$mode}. Showing {total} summary rows. {flagged} rows exceed CV threshold ({input$cv_threshold}%).")
    } else {
      glue("Physical mode: {input$mode_phys}. Showing {nrow(df)} rows. Filtered endpoint: {paste(unique(df$endpoint), collapse = ', ')}")
    }
  })
  
  # Download handlers
  output$download_data <- downloadHandler(
    filename = function() { paste0("ECOTOX_", input$mode, "_", Sys.Date(), ".xlsx") },
    content = function(file) {
      df <- summarized()
      writexl::write_xlsx(if (is.null(df)) data.frame(Message = "No data") else df, file)
    }
  )
  
  output$download_phys_data <- downloadHandler(
    filename = function() { paste0("ECOTOX_physical_", Sys.Date(), ".xlsx") },
    content = function(file) {
      df <- phys_summarized()
      writexl::write_xlsx(if (is.null(df)) data.frame(Message = "No data") else df, file)
    }
  )
  
  # ============================================================================
  # REGRESSION ANALYSIS MODULE - Using modeling.R functions
  # ============================================================================
  
  # Dynamic endpoint selection for regression
  observe({
    req(input$regression_dataset)
    if (identical(input$regression_dataset, "assay")) {
      if (exists("final_data") && nrow(final_data) > 0) {
        choices <- sort(unique(final_data$assay_type))
      } else {
        choices <- ""
      }
    } else {
      if (exists("physical_master") && nrow(physical_master) > 0) {
        choices <- sort(unique(physical_master$endpoint))
      } else {
        choices <- ""
      }
    }
    
    if (length(choices) == 0) choices <- ""
    updateSelectInput(session, "regression_endpoint", choices = choices, selected = choices[1])
  })
  
  # Dynamic endpoint selection for recovery
  observe({
    req(input$recovery_dataset)
    if (identical(input$recovery_dataset, "assay")) {
      if (exists("final_data") && nrow(final_data) > 0) {
        choices <- sort(unique(final_data$assay_type))
      } else {
        choices <- ""
      }
    } else {
      if (exists("physical_master") && nrow(physical_master) > 0) {
        choices <- sort(unique(physical_master$endpoint))
      } else {
        choices <- ""
      }
    }
    
    if (length(choices) == 0) choices <- ""
    updateSelectInput(session, "recovery_endpoint", choices = choices, selected = choices[1])
  })
  
  # --- Build random-effects formula string safely ---
  build_random_effects <- function(df, intercepts, nested_batch_tank = TRUE, slope_spec = "none") {
    parts <- character(0)
    
    # Intercepts
    add_intercept_if <- function(var, min_levels = 2L) {
      if (var %in% names(df)) {
        nlev <- dplyr::n_distinct(df[[var]])
        if (nlev >= min_levels) parts <<- c(parts, sprintf("(1|%s)", var))
      }
    }
    
    if ("tank" %in% intercepts)             add_intercept_if("tank", 3L)
    if ("experiment_batch" %in% intercepts) add_intercept_if("experiment_batch", 2L)
    if ("individual_id" %in% intercepts)    add_intercept_if("individual_id", 5L)
    
    # Nesting: experiment_batch/tank ≡ (1|experiment_batch) + (1|experiment_batch:tank)
    if (isTRUE(nested_batch_tank) &&
        all(c("experiment_batch","tank") %in% names(df))) {
      nb <- dplyr::n_distinct(df$experiment_batch)
      nt <- dplyr::n_distinct(df$tank)
      if (nb >= 2L && nt >= 3L) {
        parts <- unique(c(parts, "(1|experiment_batch)", "(1|experiment_batch:tank)"))
        # If we added nesting, drop a lone (1|tank) to avoid redundancy
        parts <- setdiff(parts, "(1|tank)")
      }
    }
    
    # Random slopes: pattern "var|group"
    if (!identical(slope_spec, "none") && grepl("\\|", slope_spec, fixed = TRUE)) {
      lhs <- sub("\\|.*$", "", slope_spec)
      grp <- sub("^.*\\|", "", slope_spec)
      lhs <- trimws(lhs); grp <- trimws(grp)
      if (all(c(lhs, grp) %in% names(df))) {
        if (dplyr::n_distinct(df[[grp]]) >= 2L) {
          parts <- c(parts, sprintf("(%s|%s)", lhs, grp))
        }
      }
    }
    
    if (length(parts) == 0L) return("")
    paste(parts, collapse = " + ")
  }
  
  # ------------------------------------------------------------------
  # Regression: model-ready data with sample/tissue filters (no is_recovery)
  # ------------------------------------------------------------------
  regression_data_cached <- reactive({
    req(input$regression_dataset, input$regression_endpoint)
    
    if (identical(input$regression_dataset, "assay")) {
      df <- final_data %>%
        create_enhanced_treatment_categories() %>%
        dplyr::mutate(outcome = calculated_concentration) %>%
        dplyr::filter(assay_type == input$regression_endpoint)
      
      # Assay sample-type filter (case-insensitive)
      if (!is.null(input$reg_sample_types) && length(input$reg_sample_types)) {
        df <- df %>%
          dplyr::filter(tolower(sample_type) %in% tolower(input$reg_sample_types))
      }
      
      # Average replicates after filtering
      df <- df %>% average_assay_replicates(outcome)
      
    } else {
      df <- physical_master %>%
        create_enhanced_treatment_categories() %>%
        dplyr::mutate(outcome = value) %>%
        dplyr::filter(endpoint == input$regression_endpoint)
      
      # Tissue filter only for mf_counts
      if (identical(input$regression_endpoint, "mf_counts") &&
          !is.null(input$reg_tissues) && length(input$reg_tissues)) {
        df <- df %>% dplyr::filter(tissue_type %in% input$reg_tissues)
      }
    }
    
    # Week filter and normalization
    wks <- input$weeks_include %||% c("1","3","5")
    df <- df %>%
      dplyr::filter(as.character(week) %in% wks) %>%
      normalize_controls_and_dose() %>%
      droplevels()
    
    df
  }) %>%
    bindCache(input$regression_dataset,
              input$regression_endpoint,
              input$reg_sample_types,
              input$reg_tissues,
              input$weeks_include)
  
  # Keep existing alias if used downstream
  regression_data <- regression_data_cached
  
  # --- Mixed Effects: endpoint choices ---
  observe({
    req(input$lmer_dataset)
    if (identical(input$lmer_dataset, "assay")) {
      choices <- if (exists("final_data") && nrow(final_data) > 0) {
        sort(unique(final_data$assay_type))
      } else character(0)
    } else {
      choices <- if (exists("physical_master") && nrow(physical_master) > 0) {
        sort(unique(physical_master$endpoint))
      } else character(0)
    }
    if (length(choices) == 0) choices <- "loading"
    updateSelectInput(session, "lmer_endpoint", choices = choices, selected = choices[1])
  })
  
  # --- Mixed Effects: tissues UI (physical mf_counts only) ---
  output$lmer_tissues_ui <- renderUI({
    req(input$lmer_dataset == "physical")
    ep <- input$lmer_endpoint %||% ""
    if (identical(ep, "mf_counts")) {
      checkboxGroupInput(
        inputId = "lmer_tissue_types",
        label   = "Tissue types:",
        choices = tissue_types_physical,   # from global.R
        selected = tissue_types_physical,
        inline  = TRUE
      )
    } else {
      # Show a disabled "All" when not applicable
      checkboxGroupInput(
        inputId = "lmer_tissue_types",
        label   = "Tissue types:",
        choices = c("All"),
        selected = "All",
        inline = TRUE
      ) %>%
        shinyjs::disable()
    }
  })
  
  # --- Mixed Effects: model-ready data (assay vs physical) ---
  lmer_data_cached <- reactive({
    req(input$lmer_dataset, input$lmer_endpoint)
    
    if (identical(input$lmer_dataset, "assay")) {
      # Assay branch
      df <- final_data %>%
        create_experiment_batch() %>%
        create_enhanced_treatment_categories() %>%
        filter(assay_type == input$lmer_endpoint) %>%
        # Filters
        { if (length(input$lmer_fiber_types)) filter(., tolower(fiber_type) %in% tolower(input$lmer_fiber_types)) else . } %>%
        { if (length(input$lmer_sample_types)) filter(., tolower(sample_type) %in% tolower(input$lmer_sample_types)) else . } %>%
        { if (length(input$lmer_treatments)) filter(., tolower(chem_treatment) %in% tolower(input$lmer_treatments)) else . } %>%
        { if (length(input$lmer_concentrations)) filter(., fiber_concentration %in% input$lmer_concentrations) else . } %>%
        { 
          wks <- input$lmer_weeks_include %||% c("1", "3", "5")
          filter(., as.character(week) %in% wks)
        } %>%
        mutate(
          outcome = calculated_concentration,
          tank = as.factor(tank),
          week = factor(as.character(week), levels = sort(unique(as.character(week)))),
          sample_type = as.factor(sample_type)
        ) %>%
        average_assay_replicates(outcome) %>%
        normalize_controls_and_dose() %>%
        mutate(grouping_var = tank) %>%
        droplevels()
      
    } else {
      # Physical branch
      df <- physical_master %>%
        create_experiment_batch() %>%
        create_enhanced_treatment_categories() %>%
        filter(endpoint == input$lmer_endpoint) %>%
        # Filters
        { if (length(input$lmer_fiber_types)) filter(., tolower(fiber_type) %in% tolower(input$lmer_fiber_types)) else . } %>%
        { if (length(input$lmer_treatments)) filter(., tolower(chem_treatment) %in% tolower(input$lmer_treatments)) else . } %>%
        { if (length(input$lmer_concentrations)) filter(., fiber_concentration %in% input$lmer_concentrations) else . } %>%
        {
          wks <- input$lmer_weeks_include %||% c("1", "3", "5")
          filter(., as.character(week) %in% wks)
        } %>%
        {
          # mf_counts gets tissue_type filtering; others ignore tissue labels
          if (identical(input$lmer_endpoint, "mf_counts") && length(input$lmer_tissue_types)) {
            filter(., tissue_type %in% input$lmer_tissue_types)
          } else .
        } %>%
        mutate(
          outcome = value,
          week = factor(as.character(week), levels = sort(unique(as.character(week)))),
          tissue_type = as.factor(tissue_type),
          # Per-individual grouping (robust to naming differences)
          individual_id = dplyr::coalesce(
            as.character(sample),
            stringr::str_extract(sample, "^[0-9]+")
          ),
          grouping_var = as.factor(individual_id)
        ) %>%
        normalize_controls_and_dose() %>%
        droplevels()
    }
    
    validate(need(nrow(df) > 0, "No rows after filters"))
    df
  }) %>%
    bindCache(input$lmer_dataset, input$lmer_endpoint, input$lmer_fiber_types,
              input$lmer_sample_types, input$lmer_tissue_types, input$lmer_treatments,
              input$lmer_concentrations, input$lmer_weeks_include)
  
  lmer_data <- lmer_data_cached
  
  # Regression model fitting
  regression_model <- eventReactive(input$run_regression, {
    tryCatch({
      # regression_data() already applies weeks_include and (after your last changes)
      # may add is_recovery only when Week6 and a baseline are both present
      df <- regression_data()
      req(nrow(df) > 0)
      
      # Does this run have 2+ distinct weeks?
      has_week <- dplyr::n_distinct(df$week) >= 2
      
      # Base formula without week terms
      form <- outcome ~ fiber_type * chem_treatment +
        is_control * fiber_type +
        dose_log10
      
      # Add week main effect and interactions only when estimable
      if (has_week) {
        form <- update(form, . ~ . + week + fiber_type:week + chem_treatment:week)
        if (isTRUE(input$include_three_way)) {
          form <- update(form, . ~ . + fiber_type:chem_treatment:week)
        }
      }
      
      # Optional dose interactions
      if (isTRUE(input$dose_by_fiber)) form <- update(form, . ~ . + dose_log10:fiber_type)
      if (isTRUE(input$dose_by_treat)) form <- update(form, . ~ . + dose_log10:chem_treatment)
      
      # Recovery terms only if present and has two levels
      has_recovery <- ("is_recovery" %in% names(df)) &&
        dplyr::n_distinct(stats::na.omit(df$is_recovery)) >= 2
      if (isTRUE(input$include_recovery_interaction) && has_recovery) {
        form <- update(form, . ~ . + is_recovery +
                         is_recovery:chem_treatment +
                         is_recovery:dose_log10)
      }
      
      stats::lm(form, data = df)
    }, error = function(e) {
      list(error = e$message)
    })
  })
  
  # --- Mixed effects model fitting (updated) ---
  lmer_model <- eventReactive(input$run_lmer, {
    tryCatch({
      df <- lmer_data_cached()
      req(nrow(df) > 0)
      
      # Fixed effects: include week terms only if estimable (≥2 levels)
      has_week <- dplyr::n_distinct(df$week) >= 2
      fixed_terms <- c("outcome ~ fiber_type + chem_treatment", "dose_log10", "is_control")
      if (has_week) {
        fixed_terms <- c(
          fixed_terms,
          "week", "fiber_type:chem_treatment", "fiber_type:week", "chem_treatment:week"
        )
        if (isTRUE(input$lmer_include_three_way) && nrow(df) > 100) {
          fixed_terms <- c(fixed_terms, "fiber_type:chem_treatment:week")
        }
      } else {
        fixed_terms <- c(fixed_terms, "fiber_type:chem_treatment")
      }
      if (isTRUE(input$lmer_dose_by_fiber)) fixed_terms <- c(fixed_terms, "dose_log10:fiber_type")
      if (isTRUE(input$lmer_dose_by_treat)) fixed_terms <- c(fixed_terms, "dose_log10:chem_treatment")
      
      fixed_formula <- paste(fixed_terms, collapse = " + ")
      
      # Random-effects builder
      rand_str <- build_random_effects(
        df = df,
        intercepts = input$lmer_random_intercepts %||% character(0),
        nested_batch_tank = isTRUE(input$lmer_nested_batch_tank),
        slope_spec = input$lmer_random_slope %||% "none"
      )
      
      # Assemble and fit
      full_formula <- if (nzchar(rand_str)) {
        as.formula(paste(fixed_formula, "+", rand_str))
      } else {
        as.formula(fixed_formula)
      }
      
      lme4::lmer(
        full_formula,
        data = df,
        control = lme4::lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 20000),
          check.conv.singular = "warning",
          calc.derivs = FALSE
        )
      )
    }, error = function(e) {
      list(error = e$message)
    })
  })
  
  
  # ============================================================================
  # ENHANCED EMMEANS WITH COLORED SIGNIFICANCE & DID ANALYSIS
  # ============================================================================
  
  # EMMEANS calculation for regular regression 
  emmeans_results <- eventReactive(input$run_emmeans, {
    tryCatch({
      mdl <- regression_model()
      validate(need(!is.null(mdl) && !inherits(mdl, "list"), "Run the regression model first"))
      
      # Dose settings for 'at'
      dose_val <- switch(input$emmeans_dose,
                         "0" = 0, "2" = log10(100), "3" = log10(1000), "4" = log10(10000))
      is_ctrl <- ifelse(dose_val == 0, 1, 0)
      at_list <- list(dose_log10 = dose_val, is_control = is_ctrl)
      
      # Is 'week' in the fitted model?
      model_terms <- attr(terms(mdl), "term.labels")
      has_week <- any(grepl("\\bweek\\b", model_terms))
      
      # Map UI to a valid emmeans spec, avoiding 'week' if absent
      emmeans_spec <- switch(input$emmeans_by,
                             "treatment"        = ~ chem_treatment | fiber_type,
                             "fiber_type"       = ~ fiber_type,
                             "week"             = if (has_week) ~ week else ~ chem_treatment | fiber_type,
                             "treatment_week"   = if (has_week) ~ chem_treatment | week + fiber_type else ~ chem_treatment | fiber_type,
                             "treatment_fiber"  = ~ chem_treatment | fiber_type,
                             "fiber_week"       = if (has_week) ~ fiber_type | week else ~ fiber_type,
                             "all_interactions" = if (has_week) ~ chem_treatment | fiber_type + week else ~ chem_treatment | fiber_type
      )
      
      emm <- emmeans::emmeans(mdl, specs = emmeans_spec, at = at_list)
      list(emmeans = emm, emmeans_df = as.data.frame(emm))
    }, error = function(e) {
      message("Emmeans error: ", e$message)
      list(error = e$message)
    })
  })
  
  # PAIRWISE COMPARISONS 
  output$pairwise_table <- DT::renderDataTable({
    # Require a fitted regression model
    m <- regression_model()
    validate(need(!is.null(m), "Run the regression model first."))
    
    # Multiple-comparison adjustment (fallback to Tukey)
    adj_method <- if (!is.null(input$adjustment_method) && nzchar(input$adjustment_method)) {
      input$adjustment_method
    } else {
      "tukey"
    }
    
    # Pull emmeans object prepared elsewhere in the regression workflow
    emm_res <- emmeans_results()
    validate(need(!is.null(emm_res$emmeans), "Calculate EM Means first"))
    emm <- emm_res$emmeans
    
    # Contrasts: control (trt.vs.ctrl with untreated as reference) or all pairwise
    contrasts <- if (identical(input$comparison_type, "control")) {
      emmeans::contrast(emm, method = "trt.vs.ctrl", ref = "untreated", adjust = adj_method)
    } else {
      emmeans::contrast(emm, method = "pairwise", adjust = adj_method)
    }
    
    # Coerce to data.frame and ensure p.value exists for consistent formatting
    dfp <- as.data.frame(contrasts)
    if (!"p.value" %in% names(dfp)) dfp$p.value <- NA_real_
    
    # Format numbers and add a significance label
    numeric_cols <- intersect(
      c("estimate", "SE", "df", "t.ratio", "p.value", "lower.CL", "upper.CL"),
      names(dfp)
    )
    
    dfp <- dfp %>%
      dplyr::mutate(
        dplyr::across(dplyr::all_of(numeric_cols), ~ signif(.x, 5)),
        Significant = dplyr::case_when(
          is.na(p.value)        ~ "Unable to compute",
          p.value < 0.001       ~ "Yes (***)",
          p.value < 0.01        ~ "Yes (**) ",
          p.value < 0.05        ~ "Yes (*)",
          TRUE                  ~ "No"
        ),
        Adjustment = adj_method
      )
    
    # Render styled table
    DT::datatable(
      dfp,
      options = list(pageLength = 20, scrollX = TRUE, dom = "Blfrtip"),
      filter = "top",
      rownames = FALSE
    ) %>%
      DT::formatStyle(
        "Significant",
        backgroundColor = DT::styleEqual(
          c("Yes (***)", "Yes (**) ", "Yes (*)", "No", "Unable to compute"),
          c("#28a745",  "#28a745",    "#28a745", "#ffffff", "#ffe6cc")
        ),
        color = DT::styleEqual(
          c("Yes (***)", "Yes (**) ", "Yes (*)", "No", "Unable to compute"),
          c("#ffffff",  "#ffffff",    "#ffffff", "#000000", "#000000")
        )
      ) %>%
      { if (all(c("estimate","SE","t.ratio") %in% names(dfp)))
        DT::formatRound(., columns = c("estimate","SE","t.ratio"), digits = 4) else . } %>%
      { if ("p.value" %in% names(dfp))
        DT::formatRound(., columns = "p.value", digits = 6) else . }
  })
  
  # ENHANCED EMMEANS TABLE WITH SIGNIFICANCE COLORS  
  output$emmeans_table <- DT::renderDataTable({
    results <- emmeans_results()
    validate(need(!is.null(results$emmeans_df), "Click 'Calculate EM Means' to generate results"))
    
    df_emm <- results$emmeans_df %>%
      mutate(across(where(is.numeric), ~ signif(.x, 4)))
    
    DT::datatable(
      df_emm, 
      options = list(pageLength = 15, scrollX = TRUE, dom = 'Blfrtip'), 
      filter = "top", 
      rownames = FALSE
    ) %>%
      DT::formatRound(columns = c("emmean", "SE", "lower.CL", "upper.CL"), digits = 4)
  })
  
  # ENHANCED EMMEANS PLOT
  output$emmeans_plot <- renderPlot({
    results <- emmeans_results()
    validate(need(!is.null(results$emmeans), "Calculate EM Means first"))
    
    tryCatch({
      # USE THE SAME APPROACH AS YOUR OLD APP - built-in emmeans plot
      plot(results$emmeans) + 
        theme_minimal() + 
        labs(
          title = "Estimated Marginal Means with 95% Confidence Intervals",
          subtitle = paste("Dose level:", input$emmeans_dose, "| Adjustment:", input$adjustment_method),
          x = "Estimated Marginal Mean"
        ) +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          panel.grid.minor = element_blank()
        )
    }, error = function(e) {
      # Fallback if emmeans plot fails - use FIXED geom_errorbar (not geom_errorbarh)
      df <- results$emmeans_df
      ggplot(df, aes(x = emmean, y = rownames(df))) +
        geom_point(size = 3, color = "#2c3e50") +
        # FIXED: Use geom_errorbar with orientation instead of deprecated geom_errorbarh
        geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), 
                      width = 0.2, orientation = "y") +
        theme_minimal() +
        labs(
          title = "Estimated Marginal Means", 
          x = "Estimated Mean", 
          y = "Groups"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
  })
  
  # ============================================================================
  # LINEAR TRENDS ANALYSIS - SERVER CODE (OPTION 2: REFIT WITH NUMERIC WEEK)
  # ============================================================================
  
  # Helper: Get filtered trend specification
  get_filtered_trend_spec <- reactive({
    req(regression_model())
    
    filter_type <- input$trend_filter_type
    
    if (filter_type == "none") {
      return(~ fiber_type + chem_treatment)
    } else if (filter_type == "fiber") {
      return(list(~ chem_treatment, fiber_type = input$trend_fiber_select))
    } else if (filter_type == "treatment") {
      return(list(~ fiber_type, chem_treatment = input$trend_treatment_select))
    } else if (filter_type == "both") {
      return(list(fiber_type = input$trend_fiber_select, 
                  chem_treatment = input$trend_treatment_select))
    }
  })
  
  # Helper: map trend_concentration to at-list for emmeans/emtrends
  trend_at_list <- reactive({
    conc <- input$trend_concentration
    if (is.null(conc) || conc == "all") return(NULL)
    
    if (conc == "0") {
      # control: dose_log10 = 0, is_control = 1
      list(dose_log10 = 0, is_control = 1)
    } else {
      # non-control dose: log10(conc), is_control = 0
      conc_num <- suppressWarnings(as.numeric(conc))
      if (isTRUE(is.na(conc_num)) || conc_num <= 0) return(NULL)
      list(dose_log10 = log10(conc_num), is_control = 0)
    }
  })
  
  # ===========================
  # LINEAR TRENDS TEXT OUTPUT
  # ===========================
  output$regression_trends <- renderPrint({
    req(input$show_trends == TRUE)
    req(regression_model())
    
    tryCatch({
      mdl <- regression_model()
      
      # Check if week is in the model
      model_terms <- attr(terms(mdl), "term.labels")
      has_week <- any(grepl("week", model_terms, ignore.case = TRUE))
      
      if (!has_week) {
        cat("Week is not included in the current model.\n")
        cat("To use time trends analysis, select week checkboxes in 'Include Weeks'.\n")
        return(invisible())
      }
      
      # Get the original data and formula from the model
      model_data <- mdl$model
      original_formula <- formula(mdl)
      
      # Check if week variable exists
      if (!"week" %in% names(model_data)) {
        cat("Week variable not found in model data.\n")
        return(invisible())
      }
      
      # Create a CLEAN dataset with numeric week
      trend_data <- model_data
      
      if (is.factor(trend_data$week)) {
        trend_data$week <- as.numeric(as.character(trend_data$week))
        cat("Note: Week converted from factor to numeric for trend analysis\n")
        cat("Week values: ", paste(sort(unique(trend_data$week)), collapse = ", "), "\n\n")
      }
      
      cat("=== LINEAR TRENDS OVER TIME ===\n\n")
      
      # Display what's being analyzed
      if (input$trend_filter_type == "none") {
        cat("Analyzing trends for: ALL fiber types and treatments\n\n")
      } else if (input$trend_filter_type == "fiber") {
        cat("Analyzing trends for: ", toupper(input$trend_fiber_select), 
            " fiber (both treatments)\n\n")
      } else if (input$trend_filter_type == "treatment") {
        cat("Analyzing trends for: ", toupper(input$trend_treatment_select), 
            " (both fiber types)\n\n")
      } else {
        cat("Analyzing trends for: ", toupper(input$trend_fiber_select), 
            " fiber, ", toupper(input$trend_treatment_select), " treatment\n\n")
      }
      
      # Report concentration choice
      if (!is.null(input$trend_concentration) && input$trend_concentration != "all") {
        cat("Evaluated at concentration:", input$trend_concentration, "mf/L\n\n")
      } else {
        cat("Evaluated across all doses (averaged)\n\n")
      }
      
      # Fit a NEW model from scratch with numeric week
      cat("Refitting model with numeric week...\n")
      mdl_numeric <- lm(original_formula, data = trend_data)
      cat("Model refit successful.\n\n")
      
      # Build 'at' list for emtrends (NULL means average over doses)
      at_list <- trend_at_list()
      
      # emtrends calls with optional 'at'
      if (input$trend_filter_type == "none") {
        trends <- emmeans::emtrends(mdl_numeric, ~ fiber_type + chem_treatment,
                                    var = "week", at = at_list)
      } else if (input$trend_filter_type == "fiber") {
        trends <- emmeans::emtrends(mdl_numeric, ~ chem_treatment,
                                    var = "week",
                                    at = c(list(fiber_type = input$trend_fiber_select), at_list))
      } else if (input$trend_filter_type == "treatment") {
        trends <- emmeans::emtrends(mdl_numeric, ~ fiber_type,
                                    var = "week",
                                    at = c(list(chem_treatment = input$trend_treatment_select), at_list))
      } else {
        trends <- emmeans::emtrends(mdl_numeric, ~ 1,
                                    var = "week",
                                    at = c(list(fiber_type = input$trend_fiber_select,
                                                chem_treatment = input$trend_treatment_select),
                                           at_list))
      }
      
      cat("--- Estimated Linear Trends (Slopes) ---\n")
      cat("How much does the outcome change per 1-week increase?\n\n")
      print(summary(trends, infer = TRUE))
      
      if (input$trend_filter_type == "none" || 
          (input$trend_filter_type != "both")) {
        cat("\n\n--- Pairwise Comparisons of Trends ---\n")
        cat("Are slopes significantly different between groups?\n\n")
        trend_pairs <- pairs(trends)
        print(summary(trend_pairs))
      }
      
      cat("\n\n=== INTERPRETATION GUIDE ===\n")
      cat("- week.trend: Change in outcome per 1-week increase\n")
      cat("- SE: Standard error of the trend estimate\n")
      cat("- t.ratio & p.value: Test if trend differs from zero\n")
      cat("- lower.CL & upper.CL: 95% confidence interval\n")
      cat("\nPositive trend = outcome increases over time\n")
      cat("Negative trend = outcome decreases over time\n")
      cat("\nNote: Trends reflect the selected concentration when a dose is chosen; otherwise they are averaged across all dose levels.\n")
      
    }, error = function(e) {
      cat("Error computing trends:\n")
      cat(conditionMessage(e), "\n\n")
      cat("This may occur if:\n")
      cat("- Model cannot be refit with numeric week\n")
      cat("- Insufficient data for selected filter\n")
      cat("- Model lacks necessary interaction terms\n")
    })
  })
  
  # ===========================
  # TREND VISUALIZATION PLOT
  # ===========================
  regression_trend_plot_reactive <- reactive({
    req(input$show_trends == TRUE)
    req(regression_model())
    mdl <- regression_model()
    
    # Guard: 'week' must be in the model
    model_terms <- attr(terms(mdl), "term.labels")
    if (!any(grepl("week", model_terms, ignore.case = TRUE))) return(NULL)
    
    at_list <- trend_at_list()
    
    if (input$trend_filter_type == "none") {
      emm_week <- emmeans::emmeans(mdl, ~ week + fiber_type + chem_treatment, at = at_list)
    } else if (input$trend_filter_type == "fiber") {
      emm_week <- emmeans::emmeans(mdl, ~ week + chem_treatment,
                                   at = c(list(fiber_type = input$trend_fiber_select), at_list))
    } else if (input$trend_filter_type == "treatment") {
      emm_week <- emmeans::emmeans(mdl, ~ week + fiber_type,
                                   at = c(list(chem_treatment = input$trend_treatment_select), at_list))
    } else {
      emm_week <- emmeans::emmeans(mdl, ~ week,
                                   at = c(list(fiber_type = input$trend_fiber_select,
                                               chem_treatment = input$trend_treatment_select), at_list))
    }
    
    emm_df <- as.data.frame(emm_week)
    emm_df$week_num <- if (is.factor(emm_df$week)) as.numeric(as.character(emm_df$week)) else emm_df$week
    
    dose_subtitle <- if (!is.null(input$trend_concentration) && input$trend_concentration != "all") {
      paste0("Evaluated at ", input$trend_concentration, " mf/L")
    } else "Averaged across doses"
    
    library(ggplot2)
    if (input$trend_filter_type == "none") {
      p <- ggplot(emm_df, aes(x = week_num, y = emmean, color = chem_treatment,
                              linetype = chem_treatment, shape = chem_treatment)) +
        geom_line(size = 1.2) +
        geom_point(size = 3.5) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, size = 0.8, alpha = 0.7) +
        facet_wrap(~ fiber_type, nrow = 1, scales = "free_y") +
        scale_color_manual(values = get_treatment_palette()) +
        theme_publication() +
        labs(title = "Estimated Marginal Means Trends Over Time",
             subtitle = paste("Linear trends; error bars = 95% CI •", dose_subtitle),
             x = "Week", y = "Estimated Mean Outcome",
             color = "Treatment", linetype = "Treatment", shape = "Treatment")
    } else if (input$trend_filter_type == "fiber") {
      p <- ggplot(emm_df, aes(x = week_num, y = emmean, color = chem_treatment,
                              linetype = chem_treatment, shape = chem_treatment)) +
        geom_line(size = 1.3) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, size = 0.9, alpha = 0.7) +
        scale_color_manual(values = get_treatment_palette()) +
        theme_publication() +
        labs(title = paste0("Trends Over Time (", toupper(input$trend_fiber_select), " Fiber)"),
             subtitle = paste("Linear trends; error bars = 95% CI •", dose_subtitle),
             x = "Week", y = "Estimated Mean Outcome",
             color = "Treatment", linetype = "Treatment", shape = "Treatment")
    } else if (input$trend_filter_type == "treatment") {
      p <- ggplot(emm_df, aes(x = week_num, y = emmean, color = fiber_type,
                              linetype = fiber_type, shape = fiber_type)) +
        geom_line(size = 1.3) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, size = 0.9, alpha = 0.7) +
        scale_color_manual(values = get_fiber_palette()) +
        theme_publication() +
        labs(title = paste0("Trends Over Time (", toupper(input$trend_treatment_select), " Treatment)"),
             subtitle = paste("Linear trends; error bars = 95% CI •", dose_subtitle),
             x = "Week", y = "Estimated Mean Outcome",
             color = "Fiber Type", linetype = "Fiber Type", shape = "Fiber Type")
    } else {
      p <- ggplot(emm_df, aes(x = week_num, y = emmean)) +
        geom_line(size = 1.5, color = "#0073C2") +
        geom_point(size = 5, color = "#0073C2") +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, size = 1, color = "#0073C2", alpha = 0.7) +
        theme_publication() +
        labs(title = paste0("Trend (", toupper(input$trend_fiber_select), ", ", toupper(input$trend_treatment_select), ")"),
             subtitle = paste("Linear trend; error bars = 95% CI •", dose_subtitle),
             x = "Week", y = "Estimated Mean Outcome")
    }
    
    p + scale_x_continuous(breaks = unique(emm_df$week_num))
  })
  
  output$regression_trend_plot <- renderPlot({
    p <- regression_trend_plot_reactive()
    if (is.null(p)) {
      plot.new()
      text(0.5, 0.5, "Week not included in model (check week checkboxes in 'Include Weeks')",
           cex = 1.3, col = "#d9534f", font = 2)
    } else p
  })
  
  output$download_regression_trend_png  <- create_plot_download_handler(regression_trend_plot_reactive, "regression_trends", "png")
  output$download_regression_trend_pdf  <- create_plot_download_handler(regression_trend_plot_reactive, "regression_trends", "pdf")
  output$download_regression_trend_tiff <- create_plot_download_handler(regression_trend_plot_reactive, "regression_trends", "tiff")
  output$download_regression_trend_svg  <- create_plot_download_handler(regression_trend_plot_reactive, "regression_trends", "svg")
  
  # =============================================================================
  # Combined Treatment Analysis (server)
  # - Uses snake_case IDs
  # - Harmonizes legacy/snake_case columns
  # - Avoids circular gating by rendering even when hidden and gating by button
  # =============================================================================
  
  # Small helper to format a quick data status line in the UI
  output$combined_data_info <- renderText({
    if (is.null(input$combined_endpoint) || identical(input$combined_endpoint, "loading")) {
      return("No data loaded")  # initial state
    }
    # Try to count rows safely
    n <- tryCatch({
      df <- combined_base_data()
      if (is.null(df)) 0L else nrow(df)
    }, error = function(e) 0L)
    if (n > 0L) paste0(format(n, big.mark = ","), " observations loaded") else "No data for selected endpoint"
  })
  
  # -----------------------------------------------------------------------------
  # REACTIVE: Base data for combined analysis (with sample/tissue filters)
  # -----------------------------------------------------------------------------
  combined_base_data <- reactive({
    req(input$combined_endpoint)
    req(input$active_dataset)
    
    cat("\n=== combined_base_data reactive triggered ===\n")
    cat("Dataset:", input$active_dataset, "\n")
    cat("Endpoint:", input$combined_endpoint, "\n")
    
    if (identical(input$active_dataset, "assay")) {
      req(final_data)
      df <- final_data %>%
        create_enhanced_treatment_categories() %>%           # normalize fiber/treatment fields
        dplyr::mutate(outcome = calculated_concentration)
      
      # Coerce outcome to numeric and drop invalid rows BEFORE averaging
      df <- df %>%
        dplyr::mutate(outcome = suppressWarnings(as.numeric(outcome))) %>%
        dplyr::filter(!is.na(outcome))
      
      # Endpoint filter then replicate averaging
      df <- df %>%
        dplyr::filter(assay_type == input$combined_endpoint)
      
      # Assay sample-type filter (case-insensitive), if provided
      if (!is.null(input$combined_sample_types) && length(input$combined_sample_types)) {
        df <- df %>%
          dplyr::filter(tolower(sample_type) %in% tolower(input$combined_sample_types))
      }
      
      # Average assay replicates after filtering
      df <- df %>% average_assay_replicates(outcome)
      
      cat("Assay data filtered. Rows:", nrow(df), "\n")
      
    } else if (identical(input$active_dataset, "physical")) {
      req(physical_master)
      df <- physical_master %>%
        create_enhanced_treatment_categories() %>%
        dplyr::mutate(outcome = value) %>%
        dplyr::filter(endpoint == input$combined_endpoint)
      
      # Tissue filter applies only to mf_counts
      if (identical(input$combined_endpoint, "mf_counts") &&
          !is.null(input$combined_tissues) && length(input$combined_tissues)) {
        df <- df %>% dplyr::filter(tissue_type %in% input$combined_tissues)
      }
      
      cat("Physical data filtered. Rows:", nrow(df), "\n")
      
    } else {
      stop("Invalid dataset selection")
    }
    
    # Weeks filter and normalization
    wks <- input$combined_weeks
    if (is.null(wks) || length(wks) == 0) {
      wks <- c("1","3","5")
      cat("No weeks selected, using default: 1, 3, 5\n")
    }
    
    df <- df %>%
      dplyr::filter(as.character(week) %in% wks) %>%
      normalize_controls_and_dose() %>%
      droplevels()
    
    cat("After week filter and normalization. Final rows:", nrow(df), "\n")
    cat("Columns:", paste(names(df), collapse = ", "), "\n\n")
    df
  }) %>%
    bindCache(input$active_dataset,
              input$combined_endpoint,
              input$combined_sample_types,
              input$combined_tissues,
              input$combined_weeks)
  
  # =============================================================================
  # Combined Treatment Analysis: Harmonize and build 14-level factor
  # =============================================================================
  combined_data <- reactive({
    df <- combined_base_data()
    req(df, nrow(df) > 0)
    
    # 1) Harmonize column names
    if ("fibertype" %in% names(df) && !"fiber_type" %in% names(df)) {
      df <- dplyr::rename(df, fiber_type = fibertype)
    }
    if ("chemtreatment" %in% names(df) && !"chem_treatment" %in% names(df)) {
      df <- dplyr::rename(df, chem_treatment = chemtreatment)
    }
    if ("iscontrol" %in% names(df) && !"is_control" %in% names(df)) {
      df <- dplyr::rename(df, is_control = iscontrol)
    }
    
    # 2) Build fiber_concentration_numeric ONLY from columns that exist
    if (!"fiber_concentration_numeric" %in% names(df)) {
      # Start with a fallback column of zeros
      fcn <- rep(0, nrow(df))
      
      # Overwrite with first available numeric source
      if ("fiber_concentration" %in% names(df)) {
        tmp <- suppressWarnings(as.numeric(df[["fiber_concentration"]]))
        if (length(tmp) == nrow(df)) fcn <- dplyr::coalesce(tmp, fcn)
      }
      # Do NOT reference dose_log10/doselog10 unless they exist
      # (combined_base_data likely doesn't have them)
      
      df[["fiber_concentration_numeric"]] <- fcn
    }
    
    # 3) Call the wrapper to create the combined factor
    df_combined <- create_combined_treatment(
      data = df,
      reference_level = input$combined_reference
    )
    
    # 4) Normalize final names
    if ("treatmentcombined" %in% names(df_combined) && !"treatment_combined" %in% names(df_combined)) {
      df_combined <- dplyr::rename(df_combined, treatment_combined = treatmentcombined)
    }
    
    # 5) *** FIXED: TRIM WHITESPACE FROM treatment_combined FACTOR ***
    # This fixes the "Control Cotton " (with trailing space) issue
    if ("treatment_combined" %in% names(df_combined)) {
      df_combined <- df_combined %>%
        dplyr::mutate(
          # Convert to character, trim whitespace, then back to factor
          treatment_combined = stringr::str_trim(as.character(treatment_combined)),
          # Re-factorize with cleaned levels
          treatment_combined = factor(
            treatment_combined,
            levels = unique(stringr::str_trim(as.character(treatment_combined)))
          )
        )
      
      # Remove "Other" or NA levels if present
      df_combined <- df_combined %>%
        dplyr::filter(!is.na(treatment_combined), treatment_combined != "Other")
      
      # Drop unused factor levels after filtering
      df_combined$treatment_combined <- droplevels(df_combined$treatment_combined)
    }
    
    # 6) Final validation
    req("outcome" %in% names(df_combined), "treatment_combined" %in% names(df_combined))
    req(nrow(df_combined) > 0)  # Also check we have data left after filtering
    
    df_combined
  })
  
  # =============================================================================
  #         OBSERVER TO DYNAMICALLY UPDATE THE REFERENCE LEVEL CHOICES
  # =============================================================================
  
  observe({
    # Trigger when combined_base_data changes (endpoint/week selection)
    req(input$combined_endpoint, length(input$combined_weeks) > 0)
    
    tryCatch({
      # Get the actual data
      data <- combined_data()
      
      if (!is.null(data) && "treatment_combined" %in% names(data)) {
        # Get actual levels from the data
        all_levels <- levels(data$treatment_combined)
        
        # Filter to only control groups (levels that start with "Control")
        control_levels <- all_levels[stringr::str_detect(all_levels, "^Control")]
        
        # If no control levels found, fall back to first level
        if (length(control_levels) == 0) {
          control_levels <- all_levels[1]
        }
        
        # Update the selectInput with only control levels
        updateSelectInput(
          session,
          "combined_reference",
          choices = control_levels,
          selected = control_levels[1]  # Default to first control
        )
        
        cat("✓ “ Reference level options updated:", paste(control_levels, collapse = ", "), "\n")
      }
    }, error = function(e) {
      # Silently fail if data not ready
      NULL
    })
  }) %>%
    bindEvent(input$combined_endpoint, input$combined_weeks)
  
  # --- Fit the model -----------------------------------------------------------
  combined_model <- eventReactive(input$run_combined_model, {
    data <- combined_data()
    req(data, nrow(data) > 0)
    
    # Get reference level from input
    ref_level <- input$combined_reference
    
    # Ensure reference level exists in data
    if (!(ref_level %in% unique(data$treatment_combined))) {
      showNotification(
        paste("Reference level", ref_level, "not found in data."),
        type = "warning"
      )
      ref_level <- levels(data$treatment_combined)[1]
    }
    
    # Relevel the factor with the selected reference
    data$treatment_combined <- stats::relevel(
      data$treatment_combined,
      ref = ref_level
    )
    
    # Check that we have enough levels
    n_lvls <- nlevels(droplevels(data$treatment_combined))
    validate(
      need(n_lvls >= 2, 
           sprintf("Not enough treatment levels after filtering (have %d). Try different weeks or endpoint.", n_lvls))
    )
    
    # Build formula
    f <- if (isTRUE(input$combined_include_week)) {
      stats::as.formula("outcome ~ treatment_combined + week")
    } else {
      stats::as.formula("outcome ~ treatment_combined")
    }
    
    # Fit model
    stats::lm(f, data = data)
  })
  
  # --- Outputs: keep ids consistent with UI and use validate() -----------------
  output$combined_treatment_levels <- renderPrint({
    data <- combined_data()
    validate(need(nrow(data) > 0, "No rows available for the selected filters."))
    level_summary <- data %>%
      dplyr::group_by(treatment_combined) %>%
      dplyr::summarise(
        n = dplyr::n(),
        mean_value = mean(outcome, na.rm = TRUE),
        sd_value   = stats::sd(outcome, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(treatment_combined)
    print(level_summary, n = Inf)
  })
  
  output$combined_model_summary <- renderPrint({
    mdl <- combined_model()
    validate(need(!is.null(mdl), "Click Run to fit the model, or check data/endpoint selection."))
    print(summary(mdl))
  })
  
  output$combined_anova <- renderPrint({
    mdl <- combined_model()
    validate(need(!is.null(mdl), "Model not available."))
    print(stats::anova(mdl))
  })
  
  # ============================================================================
  # COMBINED MODEL: EMMEANS CALCULATION (Reactive)
  # ============================================================================
  
  combined_emmeans_data <- eventReactive(input$run_combined_model, {
    mdl <- combined_model()
    validate(need(!is.null(mdl), "Model not available"))
    
    tryCatch({
      # Calculate emmeans for treatment_combined factor
      emm <- emmeans::emmeans(mdl, specs = ~ treatment_combined)
      
      # Convert to data frame
      emm_df <- as.data.frame(emm) %>%
        dplyr::mutate(dplyr::across(where(is.numeric), ~signif(.x, 4)))
      
      list(
        emmeans = emm,
        emmeans_df = emm_df
      )
    }, error = function(e) {
      message("Emmeans calculation error: ", e$message)
      list(error = e$message)
    })
  })
  
  # ============================================================================
  # COMBINED MODEL: PAIRWISE (compute once) + REACTIVE FILTERS
  # ============================================================================
  
  # 1) Compute the complete Tukey-adjusted pairwise set once per model run
  combined_pairwise_all <- eventReactive(input$run_combined_model, {
    emm_res <- combined_emmeans_data()
    validate(need(!is.null(emm_res$emmeans), "Calculate emmeans first"))
    
    tryCatch({
      emm <- emm_res$emmeans
      
      # All pairwise, no filtering here
      pw <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey")
      df <- as.data.frame(pw)
      
      # Helper parsing for filtering
      parts <- stringr::str_split_fixed(df$contrast, " - ", 2)
      a <- parts[, 1]
      b <- parts[, 2]
      
      parse_fiber <- function(x) dplyr::case_when(
        stringr::str_detect(x, regex("cotton", ignore_case = TRUE)) ~ "Cotton",
        stringr::str_detect(x, regex("\\bpet\\b|pet\\s|\\spet\\b", ignore_case = TRUE)) ~ "PET",
        TRUE ~ NA_character_
      )
      parse_trt <- function(x) dplyr::case_when(
        stringr::str_detect(x, regex("^Control| Control", ignore_case = TRUE)) ~ "Control",
        stringr::str_detect(x, regex("Untreated", ignore_case = TRUE)) ~ "Untreated",
        stringr::str_detect(x, regex("Treated", ignore_case = TRUE)) ~ "Treated",
        TRUE ~ NA_character_
      )
      parse_conc <- function(x) stringr::str_extract(x, "\\b\\d{1,5}\\b") # 0,100,1000,10000, etc.
      
      df <- df %>%
        dplyr::mutate(
          significant = dplyr::case_when(
            p.value < 0.001 ~ "***",
            p.value < 0.01  ~ "**",
            p.value < 0.05  ~ "*",
            TRUE ~ ""
          ),
          a = a, b = b,
          fiber_a = parse_fiber(a), fiber_b = parse_fiber(b),
          treat_a = parse_trt(a),   treat_b = parse_trt(b),
          conc_a  = parse_conc(a),  conc_b  = parse_conc(b)
        )
      
      df
    }, error = function(e) {
      message("Pairwise-all error: ", e$message)
      validate(need(FALSE, paste("Pairwise computation failed:", e$message)))
      NULL
    })
  })
  
  # 2) Apply lightweight, instant filters based on current UI selections
  combined_pairwise_filtered <- reactive({
    # p-adjust method for subset families; default Tukey
    p_adjust <- input$combined_p_adjust %||% "tukey"  # "tukey", "holm", "BH", etc.
    
    ftype     <- input$combined_comparison_filter %||% "all"
    ref_level <- input$combined_reference
    
    mdl <- combined_model()
    validate(need(!is.null(mdl), "Model not available"))
    
    # Helper parsers to keep downstream filters working
    parse_fiber <- function(x) dplyr::case_when(
      stringr::str_detect(x, regex("cotton", ignore_case = TRUE)) ~ "Cotton",
      stringr::str_detect(x, regex("\\bpet\\b|pet\\s|\\spet\\b", ignore_case = TRUE)) ~ "PET",
      TRUE ~ NA_character_
    )
    parse_trt <- function(x) dplyr::case_when(
      stringr::str_detect(x, regex("^Control| Control", ignore_case = TRUE)) ~ "Control",
      stringr::str_detect(x, regex("Untreated", ignore_case = TRUE)) ~ "Untreated",
      stringr::str_detect(x, regex("Treated",   ignore_case = TRUE)) ~ "Treated",
      TRUE ~ NA_character_
    )
    parse_conc <- function(x) stringr::str_extract(x, "\\b\\d{1,5}\\b")
    
    # Subset-aware recomputation so the adjustment matches the displayed family
    if (ftype %in% c("fiber", "fiber_control")) {
      emm_all  <- emmeans::emmeans(mdl, specs = ~ treatment_combined)
      all_lvls <- levels(emm_all)$treatment_combined
      fiber_pick <- input$combined_filter_fiber %||% "PET"
      
      keep_lvls <- all_lvls[stringr::str_detect(all_lvls, fiber_pick, negate = FALSE)]
      emm_fiber <- emmeans::emmeans(
        mdl, specs = ~ treatment_combined,
        exclude = setdiff(all_lvls, keep_lvls)
      )
      
      pw <- if (ftype == "fiber") {
        emmeans::contrast(emm_fiber, method = "pairwise", adjust = p_adjust)
      } else {
        fiber_control <- ifelse(
          tolower(fiber_pick) %in% c("pet","polyester","pet (polyester)"),
          "Control Pet", "Control Cotton"
        )
        emmeans::contrast(emm_fiber, method = "trt.vs.ctrl", ref = fiber_control, adjust = p_adjust)
      }
      
      out <- as.data.frame(pw)
      
      # Parse contrast into a/b and attach minimal tags
      ab <- stringr::str_split_fixed(out$contrast, " - ", 2)
      out$a <- ab[,1]
      out$b <- ab[,2]
      out$fiber_a <- parse_fiber(out$a)
      out$fiber_b <- parse_fiber(out$b)
      out$treat_a <- parse_trt(out$a)
      out$treat_b <- parse_trt(out$b)
      out$conc_a  <- parse_conc(out$a)
      out$conc_b  <- parse_conc(out$b)
      
      # Always provide 'significant' for DT styling
      out <- out %>%
        dplyr::mutate(
          significant = dplyr::case_when(
            is.na(p.value) ~ "",
            p.value < 0.001 ~ "***",
            p.value < 0.01  ~ "**",
            p.value < 0.05  ~ "*",
            TRUE ~ ""
          ),
          fiber_subset = fiber_pick
        ) %>%
        dplyr::mutate(dplyr::across(
          tidyselect::any_of(c("estimate","SE","t.ratio","p.value","lower.CL","upper.CL")),
          ~ signif(.x, 4)
        ))
      
      return(out)
    }
    
    # Otherwise: keep the original fast post-hoc display filters on the full family
    df <- combined_pairwise_all()
    validate(need(!is.null(df) && nrow(df) > 0, "No contrasts available"))
    
    out <- dplyr::as_tibble(df)
    if (ftype == "ref") {
      out <- out %>% dplyr::filter(a == ref_level | b == ref_level)
    } else if (ftype == "treat") {
      out <- out %>% dplyr::filter(!is.na(treat_a), treat_a == treat_b)
      if (!is.null(input$combined_filter_treat)) {
        out <- out %>% dplyr::filter(treat_a %in% input$combined_filter_treat)
      }
    } else if (ftype == "conc") {
      out <- out %>% dplyr::filter(conc_a == conc_b | (is.na(conc_a) & is.na(conc_b)))
      if (!is.null(input$combined_filter_conc)) {
        out <- out %>% dplyr::filter(conc_a %in% input$combined_filter_conc)
      }
    }
    
    # Guarantee 'significant' exists even after heavy filtering
    if (!"significant" %in% names(out)) {
      out <- out %>%
        dplyr::mutate(
          significant = dplyr::case_when(
            "p.value" %in% names(out) & !is.na(p.value) & p.value < 0.001 ~ "***",
            "p.value" %in% names(out) & !is.na(p.value) & p.value < 0.01  ~ "**",
            "p.value" %in% names(out) & !is.na(p.value) & p.value < 0.05  ~ "*",
            TRUE ~ ""
          )
        )
    }
    
    out %>%
      dplyr::mutate(dplyr::across(
        tidyselect::any_of(c("estimate","SE","t.ratio","p.value","lower.CL","upper.CL")),
        ~ signif(.x, 4)
      ))
  })
  
  # ============================================================================
  # OUTPUT: COMBINED EMMEANS TABLE (unchanged)
  # ============================================================================
  
  output$combined_emm <- DT::renderDataTable({
    emm_res <- combined_emmeans_data()
    validate(need(!is.null(emm_res$emmeans_df), "Click Run to calculate emmeans"))
    
    df <- emm_res$emmeans_df
    
    DT::datatable(
      df,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        dom = "Blfrtip",
        columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ),
      rownames = FALSE,
      caption = "Estimated Marginal Means for Combined Treatment Levels"
    ) %>%
      DT::formatRound(columns = c("emmean","SE","lower.CL","upper.CL"), digits = 4)
  })
  
  # Downloads for Combined EMM table (CSV/XLSX)
  output$download_combined_emm_table_csv <- create_table_download_handler(
    table_reactive = reactive({
      emm_res <- combined_emmeans_data()
      validate(need(!is.null(emm_res$emmeans_df), "No EMM data"))
      emm_res$emmeans_df
    }),
    base_filename = "combined_emm_table",
    format = "csv"
  )
  
  output$download_combined_emm_table_xlsx <- create_table_download_handler(
    table_reactive = reactive({
      emm_res <- combined_emmeans_data()
      validate(need(!is.null(emm_res$emmeans_df), "No EMM data"))
      emm_res$emmeans_df
    }),
    base_filename = "combined_emm_table",
    format = "xlsx"
  )
  
  # ============================================================================
  # OUTPUT: COMBINED EMMEANS PLOT (now honors plot filters)
  # ============================================================================
  
  combined_emm_plot_reactive <- reactive({
    emm_res <- combined_emmeans_data()
    validate(need(!is.null(emm_res$emmeans_df), "Model not available"))
    
    # Keep your existing transform + filters
    df <- dplyr::as_tibble(emm_res$emmeans_df) |>
      dplyr::mutate(
        fiber_type = dplyr::case_when(
          stringr::str_detect(treatment_combined, "Cotton") ~ "Cotton",
          stringr::str_detect(treatment_combined, regex("Pet|PET", TRUE)) ~ "PET",
          TRUE ~ "Control"
        ),
        treatment_type = dplyr::case_when(
          stringr::str_detect(treatment_combined, "Treated") ~ "Treated",
          stringr::str_detect(treatment_combined, "Untreated") ~ "Untreated",
          TRUE ~ "Control"
        ),
        conc = stringr::str_extract(treatment_combined, "\\b\\d{1,5}\\b")
      )
    
    # Same UI-driven filters
    ftype <- input$combined_emm_filter_type %||% "all"
    if (ftype == "fiber" && length(input$combined_emm_fiber)) {
      df <- df %>% dplyr::filter(fiber_type %in% input$combined_emm_fiber)
    } else if (ftype == "treat" && length(input$combined_emm_treat)) {
      df <- df %>% dplyr::filter(treatment_type %in% input$combined_emm_treat)
    } else if (ftype == "conc" && length(input$combined_emm_conc)) {
      df <- df %>% dplyr::filter(conc %in% input$combined_emm_conc | is.na(conc))
    }
    
    # Dynamic palette remains as in your code
    present <- unique(df$fiber_type)
    pal <- c()
    if ("Cotton"  %in% present) pal["Cotton"]  <- "#2E7D32"
    if ("PET"     %in% present) pal["PET"]     <- "#1565C0"
    if ("Control" %in% present) pal["Control"] <- "#757575"
    
    df$treatment_combined <- forcats::fct_reorder(df$treatment_combined, df$emmean)
    
    ggplot2::ggplot(df, ggplot2::aes(x = emmean, y = treatment_combined, color = fiber_type)) +
      ggplot2::geom_point(size = 4) +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = lower.CL, xmax = upper.CL),
        height = 0.3, orientation = "y", linewidth = 1
      ) +
      ggplot2::scale_color_manual(values = pal) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::labs(
        title = "Estimated Marginal Means with 95% Confidence Intervals",
        subtitle = paste("Reference:", input$combined_reference),
        x = "Estimated Mean (Outcome)", y = "Treatment Group", color = "Fiber Type"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 16, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 11),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "right"
      )
  })
  
  # Render (uses the new reactive)
  output$combined_emm_plot <- renderPlot({
    combined_emm_plot_reactive()
  })
  
  # Downloads (PNG/PDF/TIFF/SVG via your helper)
  output$download_combined_emm_png  <- create_plot_download_handler(combined_emm_plot_reactive, "combined_emmeans", "png")
  output$download_combined_emm_pdf  <- create_plot_download_handler(combined_emm_plot_reactive, "combined_emmeans", "pdf")
  output$download_combined_emm_tiff <- create_plot_download_handler(combined_emm_plot_reactive, "combined_emmeans", "tiff")
  output$download_combined_emm_svg  <- create_plot_download_handler(combined_emm_plot_reactive, "combined_emmeans", "svg")
  
  # ============================================================================
  # OUTPUT: COMBINED PAIRWISE TABLE (uses filtered reactive)
  # ============================================================================
  
  output$combined_pairwise <- DT::renderDataTable({
    df <- combined_pairwise_filtered()
    alpha <- input$combined_alpha
    ftype <- input$combined_comparison_filter
    
    DT::datatable(
      df %>% dplyr::select(contrast, estimate, SE, df, t.ratio, p.value, significant),
      options = list(
        pageLength = 20,
        scrollX = TRUE,
        scrollY = "500px",
        dom = "Blfrtip",
        columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ),
      rownames = FALSE,
      caption = paste0("Pairwise Comparisons (α = ", alpha, ") - Filter: ", ftype)
    ) %>%
      DT::formatRound(columns = c("estimate","SE","t.ratio","p.value"), digits = 4) %>%
      DT::formatStyle(
        'p.value',
        backgroundColor = DT::styleInterval(c(0.001, 0.01, 0.05),
                                            c('#ffebee', '#fff3e0', '#fffde7', 'white'))
      ) %>%
      DT::formatStyle(
        'significant',
        fontWeight = 'bold',
        color = DT::styleEqual(c('***','**','*',''), c('red','orange','blue','black'))
      )
  })
  
  # Downloads for Combined Pairwise table (CSV/XLSX)
  output$download_combined_pairwise_csv <- create_table_download_handler(
    table_reactive = reactive(combined_pairwise_filtered()),
    base_filename = "combined_pairwise_table",
    format = "csv"
  )
  
  output$download_combined_pairwise_xlsx <- create_table_download_handler(
    table_reactive = reactive(combined_pairwise_filtered()),
    base_filename = "combined_pairwise_table",
    format = "xlsx"
  )
  
  # ============================================================================
  # OUTPUT: COMBINED PAIRWISE PLOT (same filtered set as table)
  # ============================================================================
  
  combined_pairwise_plot_reactive <- reactive({
    df <- combined_pairwise_filtered()
    validate(need(nrow(df) > 0, "No contrasts to plot"))
    
    df_plot <- df |>
      dplyr::arrange(dplyr::desc(abs(estimate))) |>
      dplyr::slice_head(n = 30) |>
      dplyr::mutate(
        contrast = forcats::fct_reorder(contrast, estimate),
        sig_level = dplyr::case_when(
          p.value < 0.001 ~ "p < 0.001",
          p.value < 0.01  ~ "p < 0.01",
          p.value < 0.05  ~ "p < 0.05",
          TRUE            ~ "Not significant"
        )
      )
    
    ggplot2::ggplot(df_plot, ggplot2::aes(x = estimate, y = contrast, color = sig_level)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = estimate - SE * 1.96, xmax = estimate + SE * 1.96),
        height = 0.2, orientation = "y", linewidth = 0.7
      ) +
      ggplot2::scale_color_manual(
        values = c(
          "p < 0.001" = "#D32F2F",
          "p < 0.01"  = "#F57C00",
          "p < 0.05"  = "#1976D2",
          "Not significant" = "#757575"
        )
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        title = "Pairwise Comparison Estimates",
        subtitle = paste("Filter:", input$combined_comparison_filter, "| Adjustment: Tukey HSD"),
        x = "Estimated Difference", y = "Contrast", color = "Significance"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 9),
        legend.position = "bottom"
      )
  })
  
  # Render (uses the new reactive)
  output$combined_pairwise_plot <- renderPlot({
    combined_pairwise_plot_reactive()
  })
  
  # Downloads — slightly taller by default for many rows
  output$download_combined_pairwise_png  <- create_plot_download_handler(combined_pairwise_plot_reactive, "combined_pairwise", "png",  height = 10)
  output$download_combined_pairwise_pdf  <- create_plot_download_handler(combined_pairwise_plot_reactive, "combined_pairwise", "pdf",  height = 10)
  output$download_combined_pairwise_tiff <- create_plot_download_handler(combined_pairwise_plot_reactive, "combined_pairwise", "tiff", height = 10)
  output$download_combined_pairwise_svg  <- create_plot_download_handler(combined_pairwise_plot_reactive, "combined_pairwise", "svg",  height = 10)
  
  # ============================================================================
  # OUTPUT: COMBINED DIAGNOSTICS
  # ============================================================================
  
  # Diagnostics plot with safe graphics state and explicit namespaces
  output$combined_diagnostics <- renderPlot({
    mdl <- combined_model()
    validate(need(!is.null(mdl), "Run model first to see diagnostics"))
    
    # Save and restore par to avoid leaking graphics settings
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    
    # 1) Residuals vs Fitted
    plot(
      fitted(mdl), residuals(mdl),
      xlab = "Fitted Values", ylab = "Residuals",
      main = "Residuals vs Fitted", pch = 16,
      col = scales::alpha("steelblue", 0.6)   # explicit namespace
    )
    abline(h = 0, col = "red", lty = 2, lwd = 2)
    
    # 2) Normal Q-Q
    qqnorm(
      residuals(mdl), main = "Normal Q-Q Plot",
      pch = 16, col = scales::alpha("steelblue", 0.6)
    )
    qqline(residuals(mdl), col = "red", lwd = 2, lty = 2)
    
    # 3) Scale-Location
    sqrt_std_resid <- sqrt(abs(scale(residuals(mdl))))
    plot(
      fitted(mdl), sqrt_std_resid,
      xlab = "Fitted Values", ylab = expression(sqrt("|Standardized Residuals|")),
      main = "Scale-Location", pch = 16,
      col = scales::alpha("steelblue", 0.6)
    )
    lines(lowess(fitted(mdl), sqrt_std_resid), col = "red", lwd = 2)
    
    # 4) Cook's Distance
    cooks_d <- cooks.distance(mdl)
    n <- length(cooks_d)
    plot(
      seq_along(cooks_d), cooks_d, type = "h",
      xlab = "Observation Index", ylab = "Cook's Distance",
      main = "Cook's Distance (Influence)",
      col = ifelse(cooks_d > 4 / n, "red", "steelblue"), lwd = 2
    )
    abline(h = 4 / n, col = "red", lty = 2, lwd = 1.5)
    text(
      x = n * 0.7, y = max(cooks_d, na.rm = TRUE) * 0.9,
      labels = paste("Threshold:", round(4 / n, 4)),
      col = "red", cex = 0.8
    )
  })
  
  # -----------------------------------------------------------------------------
  # Download handler
  # -----------------------------------------------------------------------------
  output$download_combined_results <- downloadHandler(
    filename = function() paste0("combined_treatment_analysis_", Sys.Date(), ".txt"),
    content = function(file) {
      mdl <- combined_model()
      sink(file)
      cat("COMBINED TREATMENT ANALYSIS RESULTS\n")
      cat("Date:", as.character(Sys.Date()), "\n")
      cat("Dependent variable: outcome\n")
      cat("Reference level:", input$combined_reference, "\n\n")
      cat("MODEL SUMMARY\n")
      print(summary(mdl))
      cat("\nANOVA\n")
      print(stats::anova(mdl))
      sink()
    }
  )
  
  # ============================================================================
  # DYNAMIC ENDPOINT SELECTOR FOR COMBINED TREATMENT ANALYSIS
  # ============================================================================
  
  observe({
    req(input$active_dataset)  # CORRECT: active_dataset with underscore!
    
    cat("\n=== Updating combined_endpoint choices ===\n")
    cat("Active dataset:", input$active_dataset, "\n")
    
    choices <- character(0)
    
    if (identical(input$active_dataset, "assay")) {
      # Try to get assay endpoints from final_data
      choices <- tryCatch({
        req(final_data)
        sort(unique(final_data$assay_type))
      }, error = function(e) {
        cat("Error accessing final_data:", e$message, "\n")
        c("ACP", "AChE", "ALP", "Bradford", "CAT", "SOD")
      })
      cat("Assay endpoints:", paste(choices, collapse = ", "), "\n")
      
    } else if (identical(input$active_dataset, "physical")) {
      # Try to get physical endpoints from physical_master
      choices <- tryCatch({
        req(physical_master)
        sort(unique(physical_master$endpoint))
      }, error = function(e) {
        cat("Error accessing physical_master:", e$message, "\n")
        c("bci", "clearance_rate", "mf_counts", "respiration_rate")
      })
      cat("Physical endpoints:", paste(choices, collapse = ", "), "\n")
    }
    
    # Update the dropdown
    if (length(choices) > 0) {
      updateSelectInput(
        session, 
        "combined_endpoint", 
        choices = choices, 
        selected = choices[1]
      )
      cat("✓ “ Dropdown updated successfully with", length(choices), "choices\n")
    } else {
      updateSelectInput(
        session, 
        "combined_endpoint", 
        choices = c("No endpoints available" = "none"), 
        selected = "none"
      )
      cat("⚠  WARNING: No endpoints found!\n")
    }
  })
  
  # Status indicator
  output$combined_data_info <- renderText({
    if (is.null(input$combined_endpoint) || 
        identical(input$combined_endpoint, "none") ||
        identical(input$combined_endpoint, "loading")) {
      return("⚠  No data loaded")
    }
    
    data_count <- tryCatch({
      req(combined_base_data())
      nrow(combined_base_data())
    }, error = function(e) 0)
    
    if (data_count > 0) {
      paste0("✓ “ ", data_count, " observations loaded")
    } else {
      "⚠  No data for selected endpoint"
    }
  })
  
  # EMMEANS calculation for mixed effects models
  lmer_emmeans_results <- eventReactive(input$lmer_run_emmeans, {
    tryCatch({
      # Require that the lmer model exists first
      model <- lmer_model()
      validate(need(!is.null(model) && !inherits(model, "list"), "Run the mixed effects model first"))
      
      # Calculate dose value based on selection
      dose_val <- switch(input$lmer_emmeans_dose,
                         "0" = 0,           # Control level
                         "2" = log10(100),  # Low dose
                         "3" = log10(1000), # Medium dose
                         "4" = log10(10000)) # High dose
      
      is_ctrl <- ifelse(dose_val == 0, 1, 0)
      at_list <- list(dose_log10 = dose_val, is_control = is_ctrl)
      
      # FIXED: Better emmeans specifications to avoid convergence issues
      emmeans_spec <- switch(input$lmer_emmeans_by,
                             "treatment" = ~ chem_treatment,           # Simple: just treatment
                             "fiber_type" = ~ fiber_type,             # Simple: just fiber
                             "week" = ~ week,                         # Simple: just week
                             "treatment_week" = ~ chem_treatment + week,      # Additive, not interaction
                             "treatment_fiber" = ~ chem_treatment + fiber_type, # Additive, not interaction
                             "fiber_week" = ~ fiber_type + week,       # Additive, not interaction
                             "all_interactions" = ~ chem_treatment + fiber_type + week) # Additive only
      
      # ENHANCED: Try different approaches if emmeans fails
      emm <- tryCatch({
        # Primary approach: Use full specification
        emmeans(model, specs = emmeans_spec, at = at_list)
      }, error = function(e1) {
        message("Primary emmeans failed, trying simplified approach: ", e1$message)
        tryCatch({
          # Fallback 1: Remove the at= specification (use model's reference levels)
          emmeans(model, specs = emmeans_spec)
        }, error = function(e2) {
          message("Simplified emmeans failed, trying minimal approach: ", e2$message)
          # Fallback 2: Use simplest specification possible
          emmeans(model, specs = ~ chem_treatment)
        })
      })
      
      list(
        emmeans = emm,
        emmeans_df = as.data.frame(emm),
        specification_used = deparse(emmeans_spec)
      )
    }, error = function(e) {
      message("All emmeans approaches failed: ", e$message)
      list(error = e$message)
    })
  })
  
  # Mixed effects emmeans outputs  
  output$lmer_emmeans_table <- DT::renderDataTable({
    results <- lmer_emmeans_results()
    validate(need(!is.null(results$emmeans_df), "Click 'Calculate EM Means' for mixed effects"))
    
    df_emm <- results$emmeans_df %>%
      mutate(across(where(is.numeric), ~ signif(.x, 4)))
    
    DT::datatable(
      df_emm,
      options = list(pageLength = 10, scrollX = TRUE),
      filter = "top", rownames = FALSE
    ) %>%
      formatRound(columns = c("emmean", "SE", "lower.CL", "upper.CL"), digits = 4)
  })
  
  output$lmer_pairwise_table <- DT::renderDataTable({
    m <- lmer_model()
    validate(need(!is.null(m), "Run the mixed effects model first."))
    
    adj_method <- input$lmer_adjustment_method %||% "tukey"
    
    if (input$lmer_comparison_type == "did") {
      # DID analysis for mixed effects - simplified version
      tryCatch({
        dose_val <- as.numeric(input$lmer_emmeans_dose %||% "0")
        is_ctrl <- ifelse(dose_val == 0, 1, 0)
        
        emm_tw <- emmeans(m, 
                          specs = ~ chem_treatment | fiber_type + week,
                          at = list(dose_log10 = dose_val, is_control = is_ctrl))
        
        contrasts <- contrast(emm_tw, method = "pairwise", adjust = adj_method)
        dfp <- as.data.frame(contrasts)
        
      }, error = function(e) {
        dfp <- data.frame(
          Message = paste0("DiD Error for Mixed Effects: ", e$message),
          stringsAsFactors = FALSE
        )
      })
      
    } else {
      # Regular pairwise or treatment vs control for mixed effects
      emm_res <- lmer_emmeans_results()
      validate(need(!is.null(emm_res$emmeans), "Calculate EM Means first"))
      
      emm <- emm_res$emmeans
      
      if (input$lmer_comparison_type == "control") {
        contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = "untreated", adjust = adj_method)
      } else {
        contrasts <- contrast(emm, method = "pairwise", adjust = adj_method)
      }
      
      dfp <- as.data.frame(contrasts)
    }
    
    # Add significance column and formatting
    if (!"p.value" %in% names(dfp)) {
      dfp$p.value <- NA_real_
    }
    
    dfp <- dfp %>%
      mutate(
        across(where(is.numeric), ~ signif(.x, 5)),
        Significant = case_when(
          is.na(p.value) ~ "Unable to compute",
          p.value < 0.001 ~ "Yes (***)",
          p.value < 0.01 ~ "Yes (**)", 
          p.value < 0.05 ~ "Yes (*)",
          TRUE ~ "No"
        ),
        Adjustment = adj_method
      )
    
    # Create colored DataTable with FIXED column names
    dt <- DT::datatable(
      dfp,
      options = list(pageLength = 15, scrollX = TRUE), 
      filter = "top", rownames = FALSE
    ) %>%
      DT::formatStyle(
        "Significant",
        backgroundColor = DT::styleEqual(
          c("Yes (***)", "Yes (**)", "Yes (*)", "No", "Unable to compute"),
          c("#28a745", "#28a745", "#28a745", "#ffffff", "#ffe6cc")  # Green for significant
        ),
        color = DT::styleEqual(
          c("Yes (***)", "Yes (**)", "Yes (*)", "No", "Unable to compute"),
          c("#ffffff", "#ffffff", "#ffffff", "#000000", "#000000")   # White text on green
        )
      ) %>%
      DT::formatRound(columns = "p.value", digits = 6)
    
    # FIXED: Check which columns actually exist before formatting
    col_names <- names(dfp)
    
    if ("estimate" %in% col_names) {
      dt <- dt %>% DT::formatRound(columns = "estimate", digits = 4)
    }
    if ("SE" %in% col_names) {
      dt <- dt %>% DT::formatRound(columns = "SE", digits = 4)
    }
    if ("t.ratio" %in% col_names) {
      dt <- dt %>% DT::formatRound(columns = "t.ratio", digits = 4)
    }
    if ("z.ratio" %in% col_names) {
      dt <- dt %>% DT::formatRound(columns = "z.ratio", digits = 4)
    }
    
    dt
  })
  
  # ============================================================================
  # 3. ENHANCED MIXED EFFECTS PLOT WITH BETTER ERROR HANDLING
  # ============================================================================
  
  output$lmer_emmeans_plot <- renderPlot({
    results <- lmer_emmeans_results()
    validate(need(!is.null(results$emmeans), "Calculate EM Means for mixed effects first"))
    validate(need(is.null(results$error), paste("Emmeans error:", results$error)))
    
    tryCatch({
      # SAME AS OLD APP - use built-in emmeans plot first
      plot(results$emmeans) + 
        theme_minimal() + 
        labs(title = "Mixed Effects Estimated Marginal Means with 95% CI",
             subtitle = paste("Dose level:", input$lmer_emmeans_dose, 
                              "| Adjustment:", input$lmer_adjustment_method),
             x = "Estimated Marginal Mean") +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          panel.grid.minor = element_blank()
        )
    }, error = function(e) {
      message("Built-in emmeans plot failed, using manual approach: ", e$message)
      
      # Enhanced fallback with better error handling
      df_plot <- results$emmeans_df
      
      # Remove rows with missing or invalid confidence intervals
      valid_rows <- !is.na(df_plot$emmean) & 
        !is.na(df_plot$lower.CL) & 
        !is.na(df_plot$upper.CL) &
        is.finite(df_plot$lower.CL) & 
        is.finite(df_plot$upper.CL)
      
      if (sum(valid_rows) == 0) {
        # Fallback: create simple point plot without confidence intervals
        ggplot(df_plot, aes(x = emmean, y = rownames(df_plot))) +
          geom_point(size = 4, color = "#e74c3c") +
          theme_minimal() +
          labs(title = "Mixed Effects Estimated Marginal Means (No CI due to model issues)",
               subtitle = paste("Dose level:", input$lmer_emmeans_dose, "| Some confidence intervals invalid"),
               x = "Estimated Marginal Mean", y = "Groups") +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 11, color = "red"),
            axis.text.y = element_text(size = 10)
          )
      } else {
        # Use valid rows only with FIXED geom_errorbar
        df_clean <- df_plot[valid_rows, ]
        
        ggplot(df_clean, aes(x = emmean, y = factor(rownames(df_clean)))) +
          geom_point(size = 3, color = "#2c3e50") +
          # FIXED: Use geom_errorbar with orientation instead of deprecated geom_errorbarh
          geom_errorbar(aes(xmin = pmax(lower.CL, emmean - 3*SE), 
                            xmax = pmin(upper.CL, emmean + 3*SE)), 
                        width = 0.2, orientation = "y", color = "#3498db") +
          theme_minimal() +
          labs(title = "Mixed Effects Estimated Marginal Means with 95% CI",
               subtitle = paste("Dose level:", input$lmer_emmeans_dose,
                                "| Adjustment:", input$lmer_adjustment_method),
               x = "Estimated Marginal Mean", y = "Groups") +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 12),
            axis.text.y = element_text(size = 10),
            panel.grid.minor = element_blank()
          )
      }
    })
  })
  
  # ============================================================================
  #                           REGRESSION OUTPUTS
  # ============================================================================
  
  output$regression_model_summary <- renderPrint({
    m <- regression_model()
    validate(need(!is.null(m) && !inherits(m, "list"), "Click 'Run Regression Model' to fit the model."))
    summary(m)
  })
  
  output$regression_anova <- DT::renderDataTable({
    m <- regression_model()
    validate(need(!is.null(m) && !inherits(m, "list"), ""))
    a <- broom::tidy(car::Anova(m))
    DT::datatable(a, options = list(pageLength = 10), filter = "top", rownames = FALSE)
  })
  
  # ============================================================================
  #                          MIXED EFFECTS OUTPUTS
  # ============================================================================
  
  output$lmer_model_summary <- renderPrint({
    m <- lmer_model()
    validate(need(!is.null(m) && !inherits(m, "list"), "Click 'Run Mixed Model' to fit the model."))
    summary(m)
  })
  
  output$lmer_anova <- DT::renderDataTable({
    m <- lmer_model()
    validate(need(!is.null(m) && !inherits(m, "list"), ""))
    a <- broom.mixed::tidy(car::Anova(m))
    DT::datatable(a, options = list(pageLength = 10), filter = "top", rownames = FALSE)
  })
  
  # Model comparison output 
  output$model_compare <- DT::renderDataTable({
    # Simple model comparison between lm and lmer
    lm_model <- regression_model()
    lmer_model_obj <- lmer_model()
    
    validate(need(!is.null(lm_model) && !inherits(lm_model, "list") && 
                    !is.null(lmer_model_obj) && !inherits(lmer_model_obj, "list"), 
                  "Run both models to compare"))
    
    # Extract model comparison metrics
    lm_aic <- AIC(lm_model)
    lmer_aic <- AIC(lmer_model_obj)
    
    comparison_df <- data.frame(
      Model = c("Linear Model (lm)", "Mixed Effects (lmer)"),
      AIC = c(lm_aic, lmer_aic),
      BIC = c(BIC(lm_model), BIC(lmer_model_obj)),
      LogLik = c(as.numeric(logLik(lm_model)), as.numeric(logLik(lmer_model_obj))),
      df = c(length(coef(lm_model)), length(fixef(lmer_model_obj))),
      Better_AIC = c(lm_aic < lmer_aic, lmer_aic < lm_aic)
    )
    
    DT::datatable(comparison_df, options = list(pageLength = 5, dom = 't'), rownames = FALSE) %>%
      formatRound(columns = c("AIC", "BIC", "LogLik"), digits = 2)
  })
  
  
  # ============================================================================
  #                              RECOVERY SERVER
  # ============================================================================
  
  # state bucket if not already present
  recovery_state <- recovery_state %||% reactiveVal(NULL)
  
  # Render Recovery tissues control dynamically based on endpoint
  output$recovery_tissues_ui <- renderUI({
    req(input$recovery_dataset == "physical", input$recovery_endpoint)
    
    ep <- input$recovery_endpoint
    # Pull rows for the active endpoint
    df <- tryCatch(
      physical_master[physical_master$endpoint == ep, , drop = FALSE],
      error = function(e) NULL
    )
    if (is.null(df) || nrow(df) == 0) {
      # No data: present disabled All
      return(shinyjs::disabled(
        checkboxGroupInput("recovery_tissues", "Tissues:", choices = c("All" = "all"), selected = "all")
      ))
    }
    
    # Derive label-like values strictly from label fields
    raw_lab <- coalesce_first_existing(df, c("sample_type", "tissue_type", "tissue"))
    lab <- tolower(trimws(raw_lab))
    lab <- lab[!(is.na(lab) | lab == "" | is_numeric_like(lab))]
    choices <- sort(unique(lab))
    
    # Only mf_counts should expose tissue labels; everything else shows "All"
    has_real_labels <- identical(ep, "mf_counts") && length(choices) > 0
    
    if (has_real_labels) {
      checkboxGroupInput("recovery_tissues", "Tissues:", choices = choices, selected = choices)
    } else {
      shinyjs::disabled(
        checkboxGroupInput("recovery_tissues", "Tissues:", choices = c("All" = "all"), selected = "all")
      )
    }
  })
  
  # Optional hint to clarify behavior for non-mf_counts endpoints
  output$recovery_tissues_hint <- renderUI({
    req(input$recovery_dataset == "physical", input$recovery_endpoint)
    if (!identical(input$recovery_endpoint, "mf_counts")) {
      tags$small(style = "color:#6c757d;", "This endpoint has no tissue/sample labels; using All.")
    } else {
      NULL
    }
  })
  
  # Update week selectors whenever dataset or endpoint changes
  observeEvent(list(input$recovery_dataset, input$recovery_endpoint), ignoreInit = FALSE, {
    ds <- input$recovery_dataset
    ep <- input$recovery_endpoint
    
    # Pull weeks from the appropriate data source
    wks <- tryCatch({
      if (identical(ds, "assay")) {
        req(exists("final_data"))
        sort(unique(suppressWarnings(as.integer(final_data$week[final_data$assay_type == ep]))))
      } else {
        req(exists("physical_master"))
        sort(unique(suppressWarnings(as.integer(physical_master$week[physical_master$endpoint == ep]))))
      }
    }, error = function(e) integer())
    
    # Fallback if none found
    if (length(wks) == 0L || all(is.na(wks))) wks <- c(1L, 5L, 6L)
    
    # Update both selectors with character choices so Shiny can match selected
    updateSelectInput(session, "recovery_baseline_week",
                      choices = as.character(wks),
                      selected = as.character(min(wks, na.rm = TRUE)))
    updateSelectInput(session, "recovery_recovery_week",
                      choices = as.character(wks),
                      selected = as.character(max(wks, na.rm = TRUE)))
  })
  
  # Mini utilities (local to Recovery, no external helpers)
  canon_fiber <- function(x) {
    xl <- tolower(trimws(as.character(x)))
    dplyr::case_when(
      xl %in% c("cot", "cotton") ~ "cotton",
      xl %in% c("pet", "polyester", "poly") ~ "pet",
      TRUE ~ xl
    )
  }
  
  canon_treatment <- function(x) {
    xl <- tolower(trimws(as.character(x)))
    dplyr::case_when(
      xl %in% c("untreated", "control", "ctrl") ~ "untreated",
      xl %in% c("treated") ~ "treated",
      TRUE ~ NA_character_
    )
  }
  
  canon_conc_labels <- function(x) {
    v <- suppressWarnings(as.numeric(as.character(x)))
    allowed <- c(0, 100, 1000, 10000)
    snapped <- vapply(v, function(y) {
      if (!is.finite(y)) return(NA_real_)
      idx <- which.min(abs(allowed - y))
      tol <- max(1, 0.01 * allowed[idx])
      if (abs(allowed[idx] - y) <= tol) allowed[idx] else NA_real_
    }, numeric(1))
    as.character(as.integer(snapped))
  }
  
  nlev2 <- function(x) {
    length(unique(stats::na.omit(as.character(x))))
  }
  
  # Populate endpoints on dataset change
  observe({
    if (identical(input$recovery_dataset, "assay")) {
      req(exists("final_data"))
      eps <- sort(unique(as.character(final_data$assay_type)))
    } else {
      req(exists("physical_master"))
      eps <- sort(unique(as.character(physical_master$endpoint)))
    }
    updateSelectInput(session, "recovery_endpoint", choices = eps, selected = eps[1])
  })
  
  observe({
    req(exists("physical_master"), nrow(physical_master) > 0)
    choices <- sort(unique(physical_master$endpoint))
    selected <- if ("mf_counts" %in% choices) "mf_counts" else choices[1]
    updateSelectInput(session, "translocation_endpoint", choices = choices, selected = selected)
  })
  
  # Main recovery analysis (patched for robust physical sample_type handling)
  observeEvent(input$run_recovery, {
    message("=== RECOVERY START (patched) ===")
    recovery_state(NULL)
    
    # READ FROM UI 
    wk_base <- safe_int1(input$recovery_baseline_week)
    wk_reco <- safe_int1(input$recovery_recovery_week)
    
    # Basic validation
    validate(need(!is.na(wk_base) && !is.na(wk_reco), "Select both baseline and recovery weeks"))
    validate(need(wk_base != wk_reco, "Baseline and recovery weeks must differ"))
    
    ds <- input$recovery_dataset
    ep <- input$recovery_endpoint
    # 1) Base data per dataset with robust sample_type
    if (identical(ds, "assay")) {
      req(exists("final_data"))
      df <- final_data %>%
        dplyr::filter(.data$assay_type == ep) %>%
        dplyr::mutate(
          outcome = calculated_concentration,
          week = as.integer(week),
          fiber_type = canon_fiber(fiber_type),
          chem_treatment = canon_treatment(treatment),
          sample_type = tolower(trimws(as.character(sample_type))),
          fiber_concentration_label = canon_conc_labels(fiber_concentration)
        )
    } else {
      req(exists("physical_master"))
      df <- physical_master %>%
        dplyr::filter(.data$endpoint == ep) %>%
        { d <- .;
        # Derive label-like sample_type only from sample_type/tissue_type/tissue
        raw_lab <- coalesce_first_existing(d, c("sample_type","tissue_type","tissue"))
        lab <- tolower(trimws(raw_lab))
        # Collapse empty or numeric-only labels to "all" so downstream filters are no-ops
        lab[is.na(lab) | lab == "" | is_numeric_like(lab)] <- "all"
        d$sample_type <- lab
        # Safe fallbacks for missing columns
        if (!("fiber_concentration" %in% names(d))) d$fiber_concentration <- NA_real_
        if (!("treatment" %in% names(d)) && "chem_treatment" %in% names(d)) d$treatment <- d$chem_treatment
        d
        } %>%
        dplyr::mutate(
          outcome = value,
          week = as.integer(week),
          fiber_type = canon_fiber(fiber_type),
          chem_treatment = canon_treatment(treatment),
          fiber_concentration_label = canon_conc_labels(fiber_concentration)
        )
    }
    message("[RECOVERY] Base rows: ", nrow(df))
    
    # 2) Filters (fiber, treatment)
    if (!is.null(input$recovery_fiber_types) && length(input$recovery_fiber_types))
      df <- dplyr::filter(df, fiber_type %in% input$recovery_fiber_types)
    
    if (!is.null(input$recovery_treatments) && length(input$recovery_treatments))
      df <- dplyr::filter(df, chem_treatment %in% input$recovery_treatments)
    
    # 3) Sample/tissue filters with guards
    if (identical(ds, "assay")) {
      if (!is.null(input$recovery_samples) && length(input$recovery_samples))
        df <- dplyr::filter(df, tolower(sample_type) %in% tolower(input$recovery_samples))
    } else {
      # Only filter when real labels exist and user didn't select "all"
      label_exists <- "sample_type" %in% names(df) &&
        any(!(is.na(df$sample_type) | trimws(df$sample_type) == "")) &&
        any(df$sample_type != "all", na.rm = TRUE)
      if (label_exists && !is.null(input$recovery_tissues) && length(input$recovery_tissues)) {
        sel <- tolower(trimws(input$recovery_tissues))
        if (!("all" %in% sel)) {
          df <- dplyr::filter(df, tolower(sample_type) %in% sel)
        }
      }
    }
    
    # 4) Concentration (discrete labels) with numeric copy (guarded)
    if ("fiber_concentration_label" %in% names(df)) {
      if (identical(ds, "assay")) {
        if (!is.null(input$recovery_concentration) && length(input$recovery_concentration))
          df <- dplyr::filter(df, fiber_concentration_label %in% input$recovery_concentration)
      } else {
        if (!is.null(input$recovery_concentration_physical) && length(input$recovery_concentration_physical))
          df <- dplyr::filter(df, fiber_concentration_label %in% input$recovery_concentration_physical)
      }
      df <- df %>%
        dplyr::mutate(fiber_concentration = suppressWarnings(as.numeric(fiber_concentration_label)))
    }
    
    message("[RECOVERY] After filters: ", nrow(df))
    validate(need(nrow(df) > 0, "No data after filters"))
    
    # Subset to the two selected weeks and construct is_recovery with dynamic labels
    lbl_base <- paste0("Week", wk_base)
    lbl_reco <- paste0("Week", wk_reco)
    df <- df %>%
      dplyr::filter(week %in% c(wk_base, wk_reco)) %>%
      dplyr::mutate(
        is_recovery = factor(ifelse(week == wk_reco, lbl_reco, lbl_base),
                             levels = c(lbl_base, lbl_reco)),
        fiber_type = factor(fiber_type),
        chem_treatment = factor(chem_treatment),
        sample_type = factor(sample_type)
      ) %>%
      droplevels()
    
    validate(need(dplyr::n_distinct(df$is_recovery) == 2, "Both selected weeks must be present after filters"))
    
    # 6) Hierarchical week coverage (unchanged logic)
    present <- function(v) v %in% names(df) && nlev2(df[[v]]) > 1
    strict_vars <- c(if (present("fiber_type")) "fiber_type",
                     if (present("chem_treatment")) "chem_treatment",
                     if (present("sample_type") && identical(ds, "assay")) "sample_type")
    candidates <- list(
      strict_vars,
      setdiff(strict_vars, "sample_type"),
      setdiff(strict_vars, c("sample_type","fiber_type")),
      character(0)
    )
    
    picked <- NULL
    for (gv in candidates) {
      df_try <- if (length(gv)) {
        df %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(gv))) %>%
          dplyr::filter(dplyr::n_distinct(is_recovery) == 2) %>%
          dplyr::ungroup()
      } else {
        if (dplyr::n_distinct(df$is_recovery) == 2) df else df[0,]
      }
      if (nrow(df_try) > 0) { picked <- gv; df <- df_try; break }
    }
    msg_pick <- if (is.null(picked) || length(picked) == 0) "pooled" else paste(picked, collapse = "×")
    message("[RECOVERY] Week coverage by: ", msg_pick, " | rows: ", nrow(df))
    validate(need(nrow(df) > 0, "No strata retain both weeks; widen filters or include another dose"))
    
    # 7) Estimable formula (keep existing structure)
    include_trt  <- nlev2(df$chem_treatment) >= 2
    include_fib  <- nlev2(df$fiber_type)    >= 2
    include_samp <- nlev2(df$sample_type)   >= 2
    
    rhs <- c("is_recovery")
    if (include_trt)  rhs <- c(rhs, "chem_treatment", "is_recovery:chem_treatment")
    if (include_fib)  rhs <- c(rhs, "fiber_type", "is_recovery:fiber_type")
    if (include_samp) rhs <- c(rhs, "sample_type", "is_recovery:sample_type")
    
    form <- stats::as.formula(paste("outcome ~", paste(rhs, collapse = " + ")))
    message("[RECOVERY] Formula: ", deparse(form))
    
    # 8) Model
    mdl <- tryCatch(stats::lm(form, data = df), error = function(e) e)
    if (inherits(mdl, "error")) {
      recovery_state(list(
        error = paste("Model fitting failed:", mdl$message),
        data = df, formula = form, endpoint = ep,
        dataset = ds, baseline_week = wk_base, recovery_week = wk_reco
      ))
      return()
    }
    
    # 9) EMMeans and contrasts (unchanged; pass UI adjust)
    emm_spec <- if (include_trt && include_fib && include_samp && identical(ds, "assay")) {
      ~ is_recovery | chem_treatment + fiber_type + sample_type
    } else if (include_trt && include_fib) {
      ~ is_recovery | chem_treatment + fiber_type
    } else if (include_trt) {
      ~ is_recovery | chem_treatment
    } else if (include_fib) {
      ~ is_recovery | fiber_type
    } else if (include_samp) {
      ~ is_recovery | sample_type
    } else {
      ~ is_recovery
    }
    
    emm <- tryCatch(emmeans::emmeans(mdl, emm_spec), error = function(e) { message("[RECOVERY] emmeans failed: ", e$message); NULL })
    
    adj_method <- if (!is.null(input$recovery_p_adjust)) input$recovery_p_adjust else "BH"
    delta <- if (!is.null(emm)) tryCatch(emmeans::contrast(emm, "revpairwise", adjust = adj_method), error = function(e) NULL) else NULL
    emm_df <- if (!is.null(emm)) tryCatch(as.data.frame(emm), error = function(e) NULL) else NULL
    
    # 10) Persist
    recovery_state(list(
      model = mdl,
      data = df,
      endpoint = ep,
      dataset = ds,
      baseline_week = wk_base,
      recovery_week = wk_reco,
      emmeans = list(emmeans = emm, emmeans_df = emm_df, delta = delta),
      formula = form
    ))
    message("=== RECOVERY DONE (patched) ===")
  })
  
  # ---- Outputs ----
  output$recovery_model_summary <- renderPrint({
    res <- recovery_state()
    validate(need(!is.null(res), "Click 'Run Recovery Analysis'"))
    validate(need(is.null(res$error), res$error))
    cat("Dataset:", res$dataset, "\n")
    cat("Endpoint:", res$endpoint, "\n")
    cat("Weeks: ", res$baseline_week, " vs ", res$recovery_week, "\n", sep = "")
    cat("Formula:", deparse(res$formula), "\n\n")
    print(summary(res$model))
  })
  
  output$recovery_comparisons <- DT::renderDataTable({
    res <- recovery_state()
    validate(need(!is.null(res), ""))
    validate(need(is.null(res$error), res$error))
    validate(need(!is.null(res$emmeans$delta), "No pairwise Week6-Week5 contrasts available"))
    DT::datatable(as.data.frame(res$emmeans$delta), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
  
  output$recovery_plot <- renderPlot({
    res <- recovery_state()
    validate(need(!is.null(res), ""))
    validate(need(is.null(res$error), res$error))
    df <- res$data
    ggplot2::ggplot(df, ggplot2::aes(x = is_recovery, y = outcome, color = fiber_type, linetype = chem_treatment)) +
      ggplot2::stat_summary(fun = mean, geom = "point", size = 3) +
      ggplot2::stat_summary(fun.data = ~ c(y = mean(.), ymin = mean(.) - sd(.)/sqrt(length(.)),
                                           ymax = mean(.) + sd(.)/sqrt(length(.))), geom = "errorbar", width = 0.2) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Recovery (Week6 vs Week5)", x = "Time point", y = "Outcome")
  })
  
  output$recovery_emmeans_table <- DT::renderDataTable({
    res <- recovery_state()
    validate(need(!is.null(res), ""))
    validate(need(is.null(res$error), res$error))
    validate(need(!is.null(res$emmeans$emmeans_df), "No EM Means grid available"))
    DT::datatable(res$emmeans$emmeans_df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
  
  output$recovery_significance_table <- DT::renderDataTable({
    # Pull current recovery results
    res <- recovery_state()
    validate(need(!is.null(res), ""))                         # ensure run
    validate(need(is.null(res$error), res$error))             # show model errors if any
    validate(need(!is.null(res$emmeans$delta),
                  "No Week6 vs Week5 contrasts available for current filters."))
    
    # Coerce emmeans contrast to a data.frame
    df <- as.data.frame(res$emmeans$delta)
    
    # Numeric columns that may exist depending on model family
    numeric_cols <- intersect(
      c("estimate", "SE", "df", "t.ratio", "z.ratio", "p.value", "lower.CL", "upper.CL"),
      names(df)
    )
    
    # Chosen adjustment (default BH) and adjusted p-values
    adj_method <- if (!is.null(input$recovery_p_adjust)) input$recovery_p_adjust else "BH"
    
    df <- df %>%
      dplyr::mutate(
        dplyr::across(dplyr::all_of(numeric_cols), ~ suppressWarnings(signif(.x, 5))),
        p_value_raw = p.value,                                        # keep raw p from emmeans table
        p_value_adj = suppressWarnings(p.adjust(p_value_raw, method = adj_method)),
        stars = dplyr::case_when(                                     # star coding on adjusted p
          is.na(p_value_adj)  ~ "",
          p_value_adj < 0.001 ~ "***",
          p_value_adj < 0.01  ~ "**",
          p_value_adj < 0.05  ~ "*",
          TRUE                ~ ""
        ),
        Significant = dplyr::if_else(!is.na(p_value_adj) & p_value_adj < 0.05,
                                     paste0("Yes ", stars),
                                     paste0("No ",  stars))
      )
    
    # Friendlier label if the emmeans 'contrast' column exists
    if ("contrast" %in% names(df)) {
      df <- df %>% dplyr::rename(`Week6 vs Week5` = contrast)
    }
    
    # Columns to display (only those that exist)
    show_cols <- intersect(
      c("Week6 vs Week5", "is_recovery", "chem_treatment", "fiber_type",
        "estimate", "SE", "df", "t.ratio", "z.ratio",
        "p_value_raw", "p_value_adj", "lower.CL", "upper.CL", "Significant"),
      names(df)
    )
    out <- df %>% dplyr::select(dplyr::all_of(show_cols))
    
    # Value-specific styling for the Significant column (green for Yes, bold)
    unique_vals <- unique(out$Significant)
    color_map   <- ifelse(grepl("^Yes", unique_vals), "#1a7f37", "#444444")
    weight_map  <- ifelse(grepl("^Yes", unique_vals), "bold",     "normal")
    
    DT::datatable(
      out,
      rownames = FALSE,
      filter   = "top",
      options  = list(pageLength = 12, scrollX = TRUE, dom = "Blfrtip")
    ) %>%
      { if ("estimate"     %in% names(out)) DT::formatRound(., "estimate",     4) else . } %>%
      { if ("SE"           %in% names(out)) DT::formatRound(., "SE",           4) else . } %>%
      { if ("t.ratio"      %in% names(out)) DT::formatRound(., "t.ratio",      4) else . } %>%
      { if ("z.ratio"      %in% names(out)) DT::formatRound(., "z.ratio",      4) else . } %>%
      { if ("p_value_raw"  %in% names(out)) DT::formatRound(., "p_value_raw",  6) else . } %>%
      { if ("p_value_adj"  %in% names(out)) DT::formatRound(., "p_value_adj",  6) else . } %>%
      DT::formatStyle(
        "Significant",
        color      = DT::styleEqual(unique_vals, color_map),
        fontWeight = DT::styleEqual(unique_vals, weight_map)
      )
  })
  
  output$recovery_emmeans_plot <- renderPlot({
    res <- recovery_state()
    validate(need(!is.null(res), ""))
    validate(need(is.null(res$error), res$error))
    if (!is.null(res$emmeans$emmeans)) {
      return(plot(res$emmeans$emmeans) + ggplot2::theme_minimal() +
               ggplot2::labs(title = "Estimated Marginal Means", x = "EMMean"))
    }
    validate(need(FALSE, "No emmeans object to plot"))
  })
  
  # ----------------- Translocation: data filter & summary (mf_counts only) -----------------
  
  # Reactive to build filtered mf_counts for the selected week / concentrations / tissues
  transloc_filtered <- eventReactive(input$run_translocation, {
    req(physical_master, nrow(physical_master) > 0)
    # Require exactly 2 tissues to compare and a single week
    validate(need(!is.null(input$transloc_week), "Select a week for translocation."))
    validate(need(length(input$transloc_tissues) == 2, "Select exactly two tissues to compare."))
    validate(need(length(input$transloc_conc) >= 1, "Select at least one fiber concentration."))
    
    df <- physical_master %>%
      dplyr::filter(
        endpoint == "mf_counts",
        week == as.integer(input$transloc_week),
        fiber_concentration %in% input$transloc_conc,
        tissue_type %in% input$transloc_tissues
      )
    
    validate(need(nrow(df) > 0, "No mf_counts data for the selected week/tissues/concentrations."))
    
    # Aggregate to mean per stratum to stabilize the plots
    df %>%
      dplyr::group_by(
        fiber_type, treatment, fiber_concentration, week, tissue_type
      ) %>%
      dplyr::summarise(
        mean_count = mean(value, na.rm = TRUE),
        sd_count   = sd(value, na.rm = TRUE),
        n          = dplyr::n(),
        .groups = "drop"
      )
  })
  
  # Summarize and compute Tissue A x Tissue B difference
  transloc_summary <- reactive({
    df <- transloc_filtered()
    
    t_a <- input$transloc_tissues[[1]]
    t_b <- input$transloc_tissues[[2]]
    
    # Wide for difference
    wide <- df %>%
      dplyr::select(fiber_type, treatment, fiber_concentration, week, tissue_type, mean_count) %>%
      tidyr::pivot_wider(names_from = tissue_type, values_from = mean_count)
    
    # Guard if a tissue is missing after filtering (prevents min/max warnings)
    validate(need(all(c(t_a, t_b) %in% names(wide)),
                  "One of the tissues has no observations after filtering."))
    
    wide %>%
      dplyr::mutate(
        diff_ab = .data[[t_a]] - .data[[t_b]],
        contrast = paste0(t_a, " - ", t_b)
      )
  })
  
  # ---- Multi-week, single-tissue, multi-dose (100/1000/10000) reactive ----
  transloc_multi_filtered <- eventReactive(input$run_translocation, {
    req(physical_master, nrow(physical_master) > 0)
    # Guards: need 1+ weeks, exactly one tissue, and the target doses
    validate(need(length(input$transloc_weeks_multi) >= 1, "Select one or more weeks."))
    validate(need(!is.null(input$transloc_tissue_single) && nzchar(input$transloc_tissue_single),
                  "Select a single tissue."))
    validate(need(length(input$transloc_conc_multi) >= 1, "Select at least one concentration."))
    
    df <- physical_master %>%
      dplyr::filter(
        endpoint == "mf_counts",
        week %in% as.integer(input$transloc_weeks_multi),
        fiber_concentration %in% input$transloc_conc_multi,
        tissue_type == input$transloc_tissue_single
      ) %>%
      dplyr::group_by(fiber_type, treatment, fiber_concentration, week) %>%
      dplyr::summarise(
        mean_count = mean(value, na.rm = TRUE),
        sd_count   = sd(value,   na.rm = TRUE),
        n          = dplyr::n(),
        .groups    = "drop"
      ) %>%
      dplyr::mutate(
        # lock concentration order to 100 < 1000 < 10000
        fiber_concentration = factor(fiber_concentration, levels = c("100","1000","10000")),
        week = factor(week, levels = sort(unique(week)))
      )
    
    validate(need(nrow(df) > 0, "No mf_counts data for selected filters."))
    df
  })
  
  
  # ----------------- Translocation plots -----------------
  
  output$transloc_counts_plot <- renderPlot({
    df <- transloc_filtered()
    
    # Ensure nonempty before plotting
    validate(need(nrow(df) > 0, "No data to plot."))
    
    ggplot2::ggplot(
      df,
      ggplot2::aes(x = fiber_concentration, y = mean_count, fill = tissue_type)
    ) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = pmax(mean_count - sd_count, 0), ymax = mean_count + sd_count),
        position = ggplot2::position_dodge(width = 0.75), width = 0.2
      ) +
      ggplot2::facet_grid(fiber_type ~ treatment, labeller = ggplot2::label_both) +
      ggplot2::labs(
        title = paste0("Microfiber counts by tissue (Week ", input$transloc_week, ")"),
        x = "Fiber concentration (mF/L)", y = "Mean count", fill = "Tissue"
      ) +
      ggplot2::theme_minimal(base_size = 12)
  })
  
  output$transloc_diff_plot <- renderPlot({
    smry <- transloc_summary()
    
    validate(need(nrow(smry) > 0, "No differences to plot."))
    
    ggplot2::ggplot(
      smry,
      ggplot2::aes(x = fiber_concentration, y = diff_ab, fill = treatment)
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
      ggplot2::facet_grid(fiber_type ~ contrast, labeller = ggplot2::label_both) +
      ggplot2::labs(
        title = paste0("Tissue difference (", unique(smry$contrast), ") at Week ", input$transloc_week),
        x = "Fiber concentration (mF/L)", y = "Mean count difference (A âˆ’ B)", fill = "Treatment"
      ) +
      ggplot2::theme_minimal(base_size = 12)
  })
  
  # ---- Multi-week grouped bar plot: x = concentration, fill = week, facet by fiber/treatment ----
  output$transloc_multiweek_plot <- renderPlot({
    df <- transloc_multi_filtered()
    ggplot2::ggplot(
      df,
      ggplot2::aes(x = fiber_concentration, y = mean_count, fill = week)
    ) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = pmax(mean_count - sd_count, 0), ymax = mean_count + sd_count),
        position = ggplot2::position_dodge(width = 0.75), width = 0.2
      ) +
      ggplot2::facet_grid(fiber_type ~ treatment, labeller = ggplot2::label_both) +
      ggplot2::labs(
        title = paste0("Translocation by concentration across weeks (", input$transloc_tissue_single, ")"),
        x = "Fiber concentration (mf/L)", y = "Mean count", fill = "Week"
      ) +
      ggplot2::theme_minimal(base_size = 12)
  })
  
  output$download_transloc_multi_png  <- create_plot_download_handler(reactive(output$transloc_multiweek_plot()), "transloc_multiweek", "png")
  output$download_transloc_multi_pdf  <- create_plot_download_handler(reactive(output$transloc_multiweek_plot()), "transloc_multiweek", "pdf")
  output$download_transloc_multi_tiff <- create_plot_download_handler(reactive(output$transloc_multiweek_plot()), "transloc_multiweek", "tiff")
  output$download_transloc_multi_svg  <- create_plot_download_handler(reactive(output$transloc_multiweek_plot()), "transloc_multiweek", "svg")
  
  # ----------------- Translocation table -----------------
  
  output$transloc_summary_table <- DT::renderDataTable({
    smry <- transloc_summary()
    validate(need(nrow(smry) > 0, "No summary to show."))
    
    # Tidy table for inspection/download
    out <- smry %>%
      dplyr::arrange(fiber_type, treatment, fiber_concentration) %>%
      dplyr::select(fiber_type, treatment, week, fiber_concentration,
                    contrast, dplyr::all_of(input$transloc_tissues), diff_ab)
    
    DT::datatable(out, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # -----------------------------------------------------------------------------
  # DID BASE AND FILTERED DATA 
  # -----------------------------------------------------------------------------
  
  did_state <- reactiveVal(NULL)
  
  observeEvent(
    list(input$did_dataset, input$did_endpoint),
    ignoreInit = FALSE,
    handlerExpr = {
      df <- tryCatch(did_base_data_norm(), error = function(e) NULL)
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
        updateCheckboxGroupInput(session, "did_samples", choices = character(0), selected = NULL)
        updateCheckboxGroupInput(session, "did_tissues", choices = character(0), selected = NULL)
        shinyjs::enable("did_samples"); shinyjs::enable("did_tissues")
        output$did_samples_hint <- renderUI(NULL); output$did_tissues_hint <- renderUI(NULL)
        return(invisible(NULL))
      }
      
      if ("sample_type" %in% names(df)) {
        df$sample_type <- tolower(trimws(as.character(df$sample_type)))
      } else {
        df$sample_type <- "all"
      }
      
      if (identical(input$did_dataset, "assay")) {
        has_labels <- "sample_type" %in% names(df) && has_any_nonempty(df$sample_type)
        if (has_labels) {
          ch <- sort(unique(df$sample_type[!(is.na(df$sample_type) | df$sample_type == "")]))
          updateCheckboxGroupInput(session, "did_samples", choices = ch, selected = ch)
          shinyjs::enable("did_samples")
          output$did_samples_hint <- renderUI(NULL)
        } else {
          updateCheckboxGroupInput(session, "did_samples", choices = c("All" = "all"), selected = "all")
          shinyjs::disable("did_samples")
          output$did_samples_hint <- renderUI(
            tags$small(style = "color:#6c757d;", "This endpoint has no sample labels; using “All”.")
          )
        }
        output$did_tissues_hint <- renderUI(NULL)
      } else {
        # Physical: collapse numeric-only or excessive levels to "All"
        distinct_vals <- sort(unique(df$sample_type[!(is.na(df$sample_type) | df$sample_type == "")]))
        all_numeric <- length(distinct_vals) > 0 && all(grepl("^[0-9]+$", distinct_vals))
        too_many_levels <- length(distinct_vals) > 10
        no_labels <- length(distinct_vals) == 0 || all_numeric || (length(distinct_vals) == 1 && distinct_vals[1] == "all")
        
        if (no_labels || too_many_levels) {
          updateCheckboxGroupInput(session, "did_tissues", choices = c("All" = "all"), selected = "all")
          shinyjs::disable("did_tissues")
          output$did_tissues_hint <- renderUI(
            tags$small(style = "color:#6c757d;", "No tissue labels detected; using “All”.")
          )
        } else {
          updateCheckboxGroupInput(session, "did_tissues", choices = distinct_vals, selected = distinct_vals)
          shinyjs::enable("did_tissues")
          output$did_tissues_hint <- renderUI(NULL)
        }
        output$did_samples_hint <- renderUI(NULL)
      }
    }
  )
  
  did_base_data <- reactive({
    req(input$did_dataset, input$did_endpoint)
    
    if (identical(input$did_dataset, "assay")) {
      validate(need(exists("final_data") && nrow(final_data) > 0, "No assay data loaded"))
      final_data %>%
        create_enhanced_treatment_categories() %>%
        dplyr::filter(.data$assay_type == input$did_endpoint) %>%
        dplyr::mutate(
          outcome    = calculated_concentration,
          week       = suppressWarnings(as.integer(week)),
          fiber_type = did_canon_fiber(fiber_type),
          chem_treatment = did_canon_treatment(treatment),
          sample_type = tolower(trimws(as.character(sample_type))),
          fiber_concentration_label = did_canon_conc_labels(fiber_concentration)
        )
    } else {
      validate(need(exists("physical_master") && nrow(physical_master) > 0, "No physical data loaded"))
      physical_master %>%
        dplyr::filter(.data$endpoint == input$did_endpoint) %>%
        { d <- .;
        # Use only label-like columns; do not pull from numeric 'sample' ids
        raw_lab <- coalesce_first_existing(d, c("sample_type","tissue_type","tissue"))
        lab <- tolower(trimws(raw_lab))
        lab[is.na(lab) | lab == "" | is_numeric_like(lab)] <- "all"
        d$sample_type <- lab
        # Fallbacks so downstream logic doesn't break if columns are absent
        if (!("fiber_concentration" %in% names(d))) d$fiber_concentration <- NA_real_
        if (!("treatment" %in% names(d)) && "chem_treatment" %in% names(d)) d$treatment <- d$chem_treatment
        d
        } %>%
        create_enhanced_treatment_categories() %>%
        dplyr::mutate(
          outcome    = value,
          week       = suppressWarnings(as.integer(week)),
          fiber_type = did_canon_fiber(fiber_type),
          chem_treatment = did_canon_treatment(treatment),
          fiber_concentration_label = did_canon_conc_labels(fiber_concentration)
        )
    }
  })
  
  # Optional dose/control harmonization if present in your helpers
  did_base_data_norm <- reactive({
    df <- did_base_data()
    if ("normalize_controls_and_dose" %in% ls()) {
      df <- normalize_controls_and_dose(df)                               # keep if defined
    }
    df
  })
  
  # Apply only sample/tissue and dose filters here; mode-specific subsetting is handled later
  did_data_filtered <- reactive({
    df <- did_base_data_norm()
    
    label_exists <- "sample_type" %in% names(df) &&
      any(!(is.na(df$sample_type) | trimws(df$sample_type) == ""))
    
    if (label_exists) {
      if (identical(input$did_dataset, "assay")) {
        if (!is.null(input$did_samples) && length(input$did_samples) &&
            !("all" %in% tolower(input$did_samples))) {
          df <- df %>% dplyr::filter(sample_type %in% tolower(trimws(input$did_samples)))
        }
      } else {
        if (!is.null(input$did_tissues) && length(input$did_tissues) &&
            !("all" %in% tolower(input$did_tissues))) {
          df <- df %>% dplyr::filter(sample_type %in% tolower(trimws(input$did_tissues)))
        }
      }
    }
    
    if ("fiber_concentration_label" %in% names(df) &&
        !is.null(input$did_concentration) && length(input$did_concentration)) {
      df <- df %>%
        dplyr::filter(fiber_concentration_label %in% input$did_concentration) %>%
        dplyr::mutate(fiber_concentration = suppressWarnings(as.numeric(fiber_concentration_label)))
    }
    
    validate(need(nrow(df) > 0, "No rows after DiD pre-filters; relax your selections."))
    df %>% dplyr::filter(!is.na(outcome))
  })
  
  # -----------------------------------------------------------------------------
  # DID SERVER (consolidated, no redundant filtering)
  # -----------------------------------------------------------------------------
  
  # Populate endpoint and week choices (keeps IDs consistent with your UI)
  observe({
    if (identical(input$did_dataset, "assay")) {
      req(exists("final_data"))
      eps <- sort(unique(as.character(final_data$assay_type)))
      wks <- sort(unique(suppressWarnings(as.integer(final_data$week))))
    } else {
      req(exists("physical_master"))
      eps <- sort(unique(as.character(physical_master$endpoint)))
      wks <- sort(unique(suppressWarnings(as.integer(physical_master$week))))
    }
    if (length(eps) == 0) eps <- ""
    if (length(wks) == 0) wks <- c(1L, 3L, 5L)
    
    updateSelectInput(session, "did_endpoint", choices = eps, selected = eps[1])
    updateSelectInput(session, "did_week_select", choices = as.character(wks), selected = as.character(max(wks)))
    updateSelectInput(session, "did_week_a", choices = as.character(wks), selected = as.character(min(wks)))
    updateSelectInput(session, "did_week_b", choices = as.character(wks), selected = as.character(max(wks)))
  })
  
  # -----------------------------------------------------------------------------
  # DiD runner: handles all modes and baselines, incl. NEW untreated_ctrl baseline
  # -----------------------------------------------------------------------------
  observeEvent(input$run_did, {
    message("=== DID START (consolidated) ===")
    did_state(NULL)
    
    ds   <- input$did_dataset
    ep   <- input$did_endpoint
    mode <- input$did_mode
    if (is.null(ds) || is.null(ep) || is.null(mode)) {
      return(safe_abort("Select dataset, endpoint, and mode before running DiD."))
    }
    
    # Base (un-subset) and pre-filtered (sample/tissue + dose) frames
    df0 <- did_base_data_norm()
    df  <- did_data_filtered()
    message("[DiD] Base rows: ", nrow(df0), " | Pre-filtered rows: ", nrow(df))
    
    # Optional fiber/treatment subset filters not required by the chosen contrast
    if (!identical(mode, "across_fiber") &&
        !is.null(input$did_fibers_filter) && length(input$did_fibers_filter)) {
      df <- df %>% dplyr::filter(fiber_type %in% input$did_fibers_filter)
      message("[DiD] After optional fiber subset: ", nrow(df))
    }
    if (!(identical(mode, "across_fiber") && identical(input$did_across_fiber_contrast, "did")) &&
        !identical(mode, "across_week") &&
        !identical(mode, "across_treatment") &&
        !is.null(input$did_treatments_filter) && length(input$did_treatments_filter)) {
      df <- df %>% dplyr::filter(chem_treatment %in% input$did_treatments_filter)
      message("[DiD] After optional treatment subset: ", nrow(df))
    }
    if (nrow(df) == 0) return(safe_abort("No data after optional subsets."))
    
    # -----------------------------
    # Mode-specific subsetting
    # -----------------------------
    if (identical(mode, "across_fiber")) {
      wk <- safe_int1(input$did_week_select)
      if (is.na(wk)) return(safe_abort("Select a week for Across Fiber DiD."))
      df <- df %>%
        dplyr::filter(
          week == wk,
          fiber_type %in% c(input$did_fiber_a, input$did_fiber_b)
        )
      message("[DiD] Across_fiber rows at week ", wk, ": ", nrow(df))
      
      # Branch A: Dose-based DiD (treated(dose) - untreated(0)) difference across fibers
      if (identical(input$did_baseline_type, "dose")) {
        dose_chosen <- as.character(input$did_treated_dose)
        df <- df %>%
          dplyr::filter(
            chem_treatment %in% c("untreated","treated"),
            fiber_type %in% c(input$did_fiber_a, input$did_fiber_b),
            fiber_concentration_label %in% c("0", dose_chosen)
          ) %>%
          dplyr::mutate(
            fiber_type = factor(fiber_type, levels = c("cotton","pet")),
            chem_treatment = factor(chem_treatment, levels = c("untreated","treated")),
            fiber_concentration_label = factor(fiber_concentration_label, levels = c("0", dose_chosen)),
            sample_type = factor(sample_type)
          ) %>% droplevels()
        
        has_two <- function(v) v %in% names(df) && dplyr::n_distinct(stats::na.omit(df[[v]])) >= 2
        if (!all(c(has_two("fiber_type"), has_two("chem_treatment"), has_two("fiber_concentration_label")))) {
          return(safe_abort("Need both fibers, both treatments, and doses {0 and selected} for dose-based DiD."))
        }
        
        rhs <- "chem_treatment * fiber_type * fiber_concentration_label"
        if (has_two("sample_type")) rhs <- paste(rhs, "+ sample_type")
        form <- stats::as.formula(paste("outcome ~", rhs))
        mdl  <- tryCatch(stats::lm(form, data = df), error = function(e) e)
        if (inherits(mdl, "error")) return(safe_abort(paste("Model failed:", mdl$message)))
        
        emm3    <- safe_emmeans(mdl, ~ chem_treatment * fiber_type * fiber_concentration_label)
        emm_tbl <- if (!is.null(emm3)) tryCatch(as.data.frame(emm3), error = function(e) NULL) else NULL
        if (is.null(emm_tbl) || nrow(emm_tbl) == 0) return(safe_abort("EMMeans grid is empty for dose-based DiD."))
        
        a <- input$did_fiber_a; b <- input$did_fiber_b
        need <- rbind(
          c("treated",   a, dose_chosen),
          c("untreated", a, "0"),
          c("treated",   b, dose_chosen),
          c("untreated", b, "0")
        )
        have_key <- paste(emm_tbl$chem_treatment, emm_tbl$fiber_type, emm_tbl$fiber_concentration_label, sep = ".")
        need_key <- apply(need, 1, paste, collapse = ".")
        if (!all(need_key %in% have_key)) {
          return(safe_abort("Missing cells for dose-based DiD; widen filters or pick another dose."))
        }
        
        # Weights: (Trt,A,dose) - (Ctrl,A,0) - (Trt,B,dose) + (Ctrl,B,0)
        w <- rep(0, nrow(emm_tbl))
        sel <- function(tr, fi, doz) emm_tbl$chem_treatment == tr & emm_tbl$fiber_type == fi & emm_tbl$fiber_concentration_label == doz
        w[sel("treated",   a, dose_chosen)] <-  1
        w[sel("untreated", a, "0")]          <- -1
        w[sel("treated",   b, dose_chosen)]  <- -1
        w[sel("untreated", b, "0")]          <-  1
        
        did <- if (!is.null(emm3)) tryCatch(emmeans::contrast(emm3, method = list(did_dose = w), by = NULL), error = function(e) NULL) else NULL
        did_tbl <- if (!is.null(did)) tryCatch(as.data.frame(did), error = function(e) NULL) else NULL
        if (!is.null(did_tbl) && "p.value" %in% names(did_tbl)) {
          method <- input$did_p_adjust %||% "none"
          did_tbl$p_adjust_method <- method
          did_tbl$p_adj <- stats::p.adjust(did_tbl$p.value, method = method)
        }
        
        did_state(list(model = mdl, data = df, emmeans_tbl = emm_tbl, did_tbl = did_tbl,
                       mode = "across_fiber", baseline = "dose"))
        message("=== DID DONE (dose-based) ===")
        return(invisible(NULL))
      }
      
      # Branch B: NEW untreated_ctrl control-offset within untreated
      if (identical(input$did_baseline_type, "untreated_ctrl")) {
        dose_chosen <- as.character(input$did_treated_dose %||% "100")
        
        df_uc <- df %>%
          dplyr::filter(
            chem_treatment == "untreated",
            fiber_type %in% c(input$did_fiber_a, input$did_fiber_b),
            fiber_concentration_label %in% c("0", dose_chosen)
          ) %>%
          dplyr::mutate(
            fiber_type = factor(fiber_type, levels = c("cotton","pet")),
            fiber_concentration_label = factor(fiber_concentration_label, levels = c("0", dose_chosen)),
            sample_type = factor(sample_type)
          ) %>%
          droplevels()
        
        has_two <- function(v) v %in% names(df_uc) && dplyr::n_distinct(stats::na.omit(df_uc[[v]])) >= 2
        if (!all(c(has_two("fiber_type"), has_two("fiber_concentration_label")))) {
          return(safe_abort("Need both fibers and doses {0 and selected} for untreated control-offset DiD."))
        }
        
        rhs <- "fiber_type * fiber_concentration_label"
        if (has_two("sample_type")) rhs <- paste(rhs, "+ sample_type")
        form <- stats::as.formula(paste("outcome ~", rhs))
        mdl  <- tryCatch(stats::lm(form, data = df_uc), error = function(e) e)
        if (inherits(mdl, "error")) return(safe_abort(paste("Model failed:", mdl$message)))
        
        emm2 <- safe_emmeans(mdl, ~ fiber_type * fiber_concentration_label)
        emm_tbl <- if (!is.null(emm2)) tryCatch(as.data.frame(emm2), error = function(e) NULL) else NULL
        if (is.null(emm_tbl) || nrow(emm_tbl) == 0) {
          return(safe_abort("EMMeans grid is empty for untreated control-offset DiD."))
        }
        
        a <- input$did_fiber_a; b <- input$did_fiber_b
        have_key <- paste(emm_tbl$fiber_type, emm_tbl$fiber_concentration_label, sep = ".")
        need_key <- c(paste(a, dose_chosen, sep="."), paste(a, "0", sep="."),
                      paste(b, dose_chosen, sep="."), paste(b, "0", sep="."))
        if (!all(need_key %in% have_key)) {
          return(safe_abort("Missing untreated cells (fiber × {0, dose})."))
        }
        
        # Weights: [untreated(dose) - untreated(0)]_A - [untreated(dose) - untreated(0)]_B
        w <- rep(0, length(have_key))
        w[have_key == paste(a, dose_chosen, sep=".")] <-  1
        w[have_key == paste(a, "0", sep=".")]        <- -1
        w[have_key == paste(b, dose_chosen, sep=".")] <- -1
        w[have_key == paste(b, "0", sep=".")]        <-  1
        
        did_uc  <- if (!is.null(emm2)) tryCatch(
          emmeans::contrast(emm2, method = list(untreated_control_offset = w), by = NULL),
          error = function(e) NULL
        ) else NULL
        did_tbl <- if (!is.null(did_uc)) tryCatch(as.data.frame(did_uc), error = function(e) NULL) else NULL
        
        if (!is.null(did_tbl) && "p.value" %in% names(did_tbl)) {
          method <- input$did_p_adjust %||% "none"
          did_tbl$p_adjust_method <- method
          did_tbl$p_adj <- stats::p.adjust(did_tbl$p.value, method = method)
        }
        
        did_state(list(
          model        = mdl,
          data         = df_uc,
          emmeans_tbl  = emm_tbl,
          did_tbl      = did_tbl,
          mode         = "across_fiber",
          baseline     = "untreated_ctrl",
          dose_chosen  = dose_chosen
        ))
        message("=== DID DONE (untreated control-offset) ===")
        return(invisible(NULL))
      }
      
      # Optional single-treatment fiber contrast path (no DiD)
      if (identical(input$did_across_fiber_contrast, "single_trt") && !is.null(input$did_single_treatment)) {
        df <- df %>% dplyr::filter(chem_treatment == input$did_single_treatment)
        message("[DiD] Across_fiber single-trt subset: ", nrow(df))
      }
      
    } else if (identical(mode, "across_week")) {
      wk_a <- safe_int1(input$did_week_a); wk_b <- safe_int1(input$did_week_b)
      if (anyNA(c(wk_a, wk_b))) return(safe_abort("Select both Week A and Week B."))
      df <- df %>%
        dplyr::filter(
          week %in% c(wk_a, wk_b),
          fiber_type == input$did_fiber_select
        )
      message("[DiD] Across_week rows at weeks ", paste(c(wk_a, wk_b), collapse = ", "), ": ", nrow(df))
      
    } else { # across_treatment
      wk <- safe_int1(input$did_week_select)
      if (is.na(wk)) return(safe_abort("Select a week for Across Treatment DiD."))
      df <- df %>%
        dplyr::filter(
          week == wk,
          fiber_type == input$did_fiber_select,
          chem_treatment %in% c(input$did_trt_a, input$did_trt_b)
        )
      message("[DiD] Across_treatment rows at week ", wk, " within fiber ", input$did_fiber_select, ": ", nrow(df))
    }
    
    if (nrow(df) == 0) return(safe_abort("No rows after mode-specific subsetting."))
    
    # Smart pooling: retry without sample/tissue filter if coverage is insufficient
    has_both <- function(v) v %in% names(df) && nlev2(df[[v]]) >= 2
    needs_both_trt <- (identical(mode, "across_fiber") && identical(input$did_across_fiber_contrast, "did")) ||
      identical(mode, "across_week") || identical(mode, "across_treatment")
    needs_both_fib <- identical(mode, "across_fiber")
    
    coverage_ok <- (!needs_both_fib || has_both("fiber_type")) && (!needs_both_trt || has_both("chem_treatment"))
    if (!coverage_ok) {
      message("[DiD] Coverage insufficient; retrying without sample/tissue filter")
      df <- did_base_data_norm()
      
      # Re-apply dose filter only (if present)
      if ("fiber_concentration_label" %in% names(df) &&
          !is.null(input$did_concentration) && length(input$did_concentration)) {
        df <- df %>% dplyr::filter(fiber_concentration_label %in% input$did_concentration) %>%
          dplyr::mutate(fiber_concentration = suppressWarnings(as.numeric(fiber_concentration_label)))
      }
      
      # Re-apply mode-specific subset
      if (identical(mode, "across_fiber")) {
        wk <- safe_int1(input$did_week_select)
        df <- df %>% dplyr::filter(week == wk, fiber_type %in% c(input$did_fiber_a, input$did_fiber_b))
        if (identical(input$did_across_fiber_contrast, "single_trt") && !is.null(input$did_single_treatment)) {
          df <- df %>% dplyr::filter(chem_treatment == input$did_single_treatment)
        }
      } else if (identical(mode, "across_week")) {
        wk_a <- safe_int1(input$did_week_a); wk_b <- safe_int1(input$did_week_b)
        df <- df %>% dplyr::filter(week %in% c(wk_a, wk_b), fiber_type == input$did_fiber_select)
      } else { # across_treatment
        wk <- safe_int1(input$did_week_select)
        df <- df %>%
          dplyr::filter(week == wk, fiber_type == input$did_fiber_select,
                        chem_treatment %in% c(input$did_trt_a, input$did_trt_b))
      }
      
      coverage_ok <- (!needs_both_fib || has_both("fiber_type")) && (!needs_both_trt || has_both("chem_treatment"))
      if (!coverage_ok) return(safe_abort("Not enough data to estimate the requested contrast; widen filters."))
      message("[DiD] Coverage restored after pooling; rows: ", nrow(df))
    }
    
    # Standardize factors (drop singletons)
    df <- df %>%
      dplyr::mutate(
        week = factor(week),
        fiber_type = factor(fiber_type, levels = c("cotton","pet")),
        chem_treatment = factor(chem_treatment, levels = c("untreated","treated")),
        sample_type = factor(sample_type)
      ) %>% droplevels()
    
    # ---------------------------------------------------------------------------
    # Build formula based on varied covariates (completes the truncated section)
    # ---------------------------------------------------------------------------
    has_varied <- function(v) v %in% names(df) && length(levels(df[[v]])) >= 2
    
    if (identical(mode, "across_fiber")) {
      if (!has_varied("chem_treatment") || !has_varied("fiber_type")) {
        return(safe_abort("Need both treatments and both fibers present for Across Fiber DiD."))
      }
      rhs <- "chem_treatment * fiber_type"
      if (has_varied("sample_type")) rhs <- paste(rhs, "+ sample_type")
    } else if (identical(mode, "across_week")) {
      if (!has_varied("chem_treatment") || !has_varied("week")) {
        return(safe_abort("Need both treatments and both selected weeks present for Across Week DiD."))
      }
      rhs <- "chem_treatment * week"
      if (has_varied("sample_type")) rhs <- paste(rhs, "+ sample_type")
    } else {  # across_treatment
      if (!has_varied("chem_treatment")) {
        return(safe_abort("Need both selected treatments present for Across Treatment within this fiber."))
      }
      rhs <- "chem_treatment"
      if (has_varied("sample_type")) rhs <- paste(rhs, "+ sample_type")
    }
    form <- stats::as.formula(paste("outcome ~", rhs))
    
    # Fit model safely
    mdl <- tryCatch(stats::lm(form, data = df), error = function(e) e)
    if (inherits(mdl, "error")) return(safe_abort(paste("Model failed:", mdl$message)))
    message("[DiD] Model fit OK; n=", nrow(df))
    
    # Defensive EMMeans and contrasts for each mode
    emm_tbl <- did_tbl <- NULL
    
    if (identical(mode, "across_fiber")) {
      emm2 <- safe_emmeans(mdl, ~ chem_treatment * fiber_type)
      emm_tbl <- if (!is.null(emm2)) tryCatch(as.data.frame(emm2), error = function(e) NULL) else NULL
      
      fibers_needed <- unique(c(input$did_fiber_a, input$did_fiber_b))
      trts_needed   <- c("untreated","treated")
      have_cells <- if (!is.null(emm_tbl)) paste(emm_tbl$chem_treatment, emm_tbl$fiber_type, sep = ".") else character(0)
      need_cells <- as.vector(outer(trts_needed, fibers_needed, paste, sep = "."))
      if (!all(need_cells %in% have_cells)) {
        return(safe_abort("Selected week/filters are missing treatment-by-fiber cells; widen filters."))
      }
      
      a <- input$did_fiber_a; b <- input$did_fiber_b
      w <- rep(0, length(have_cells))
      w[have_cells == paste("treated",   a, sep = ".")] <-  1
      w[have_cells == paste("untreated", a, sep = ".")] <- -1
      w[have_cells == paste("treated",   b, sep = ".")] <- -1
      w[have_cells == paste("untreated", b, sep = ".")] <-  1
      
      did_tbl <- if (!is.null(emm2)) tryCatch(
        as.data.frame(emmeans::contrast(emm2, method = list(DiD = w), by = NULL)),
        error = function(e) NULL
      ) else NULL
      
    } else if (identical(mode, "across_week")) {
      emm2 <- safe_emmeans(mdl, ~ chem_treatment * week)
      emm_tbl <- if (!is.null(emm2)) tryCatch(as.data.frame(emm2), error = function(e) NULL) else NULL
      
      wa <- as.character(input$did_week_a); wb <- as.character(input$did_week_b)
      trts_needed <- c("untreated","treated")
      have_cells <- if (!is.null(emm_tbl)) paste(emm_tbl$chem_treatment, emm_tbl$week, sep = ".") else character(0)
      need_cells <- as.vector(outer(trts_needed, c(wa, wb), paste, sep = "."))
      if (!all(need_cells %in% have_cells)) {
        return(safe_abort("Selected weeks are missing treatment cells; widen filters."))
      }
      
      w <- rep(0, length(have_cells))
      w[have_cells == paste("treated",   wb, sep = ".")] <-  1
      w[have_cells == paste("untreated", wb, sep = ".")] <- -1
      w[have_cells == paste("treated",   wa, sep = ".")] <- -1
      w[have_cells == paste("untreated", wa, sep = ".")] <-  1
      
      did_tbl <- if (!is.null(emm2)) tryCatch(
        as.data.frame(emmeans::contrast(emm2, method = list(DiD = w), by = NULL)),
        error = function(e) NULL
      ) else NULL
      
    } else {  # across_treatment within selected fiber
      emm2 <- safe_emmeans(mdl, ~ chem_treatment)
      emm_tbl <- if (!is.null(emm2)) tryCatch(as.data.frame(emm2), error = function(e) NULL) else NULL
      
      ta <- input$did_trt_a; tb <- input$did_trt_b
      have <- if (!is.null(emm_tbl)) emm_tbl$chem_treatment else character(0)
      if (!all(c(ta, tb) %in% have)) {
        return(safe_abort("Selected treatments not both present for the chosen fiber; widen filters."))
      }
      
      w_named <- stats::setNames(c(1, -1), c(ta, tb))  # TrtA − TrtB
      did_tbl <- if (!is.null(emm2)) tryCatch(
        as.data.frame(emmeans::contrast(emm2, method = list(treatment_effect = w_named), by = NULL)),
        error = function(e) NULL
      ) else NULL
    }
    
    # p-adjust for the DiD result table (if present)
    if (!is.null(did_tbl) && "p.value" %in% names(did_tbl)) {
      method <- input$did_p_adjust %||% "none"
      did_tbl$p_adjust_method <- method
      did_tbl$p_adj <- stats::p.adjust(did_tbl$p.value, method = method)
    }
    
    # Save state for outputs
    did_state(list(model = mdl, data = df, emmeans_tbl = emm_tbl, did_tbl = did_tbl, mode = mode))
    message("=== DID DONE (consolidated) ===")
  })
  
  # Outputs (unchanged IDs)
  output$did_summary <- renderPrint({
    res <- did_state()
    if (is.null(res)) { cat("Click 'Run DiD Analysis'\n"); return() }
    if (!is.null(res$error)) { cat("Error:", res$error, "\n"); return() }
    cat("DiD Mode:", res$mode, "\n")
    print(summary(res$model))
  })
  
  output$did_table <- DT::renderDataTable({
    res <- did_state()
    if (is.null(res)) return(NULL)
    if (!is.null(res$error)) return(DT::datatable(data.frame(Message = res$error), rownames = FALSE))
    if (is.null(res$did_tbl) || nrow(res$did_tbl) == 0)
      return(DT::datatable(data.frame(Message = "No DiD rows for current settings."), rownames = FALSE))
    DT::datatable(res$did_tbl, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
  
  output$did_emm_table <- DT::renderDataTable({
    res <- did_state()
    if (is.null(res)) return(NULL)
    if (!is.null(res$error)) return(DT::datatable(data.frame(Message = res$error), rownames = FALSE))
    if (is.null(res$emmeans_tbl) || nrow(res$emmeans_tbl) == 0)
      return(DT::datatable(data.frame(Message = "No EMMeans grid for current settings."), rownames = FALSE))
    DT::datatable(res$emmeans_tbl, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
  
  output$did_plot <- renderPlot({
    res <- did_state()
    validate(need(!is.null(res), ""))
    dfp <- res$data
    validate(need(!is.null(dfp) && nrow(dfp) > 0, "No data to plot."))
    ggplot2::ggplot(dfp, ggplot2::aes(x = chem_treatment, y = outcome, color = fiber_type)) +
      ggplot2::stat_summary(fun = mean, geom = "point", size = 3,
                            position = ggplot2::position_dodge(width = 0.35)) +
      ggplot2::stat_summary(
        fun.data = ~ c(y = mean(.), ymin = mean(.) - sd(.)/sqrt(length(.)),
                       ymax = mean(.) + sd(.)/sqrt(length(.))),
        geom = "errorbar", width = 0.2,
        position = ggplot2::position_dodge(width = 0.35)
      ) +
      ggplot2::facet_wrap(~ week, nrow = 1, scales = "free_y") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "DiD visualization", x = "Treatment", y = "Outcome")
  })
    
    # Memory cleanup
    onStop(function() {
      rm(list = ls())
      gc()
    })
  
}
# ============================================================================
# RUN APPLICATION  
# ============================================================================

shinyApp(ui = ui, server = server)