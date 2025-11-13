# ============================================================================
#                                  SHINY APP
# ============================================================================

# Source all modules (loaded once at startup)
source("global.R", local = TRUE) # Libraries, data, config
source("helpers.R", local = TRUE) # Data processing functions
source("modeling.R", local = TRUE) # Statistical modeling functions

# ---------------------------------------------------------------------------
# INLINE EXPORT & STYLE HELPERS (place near top of app.R after libraries)
# ---------------------------------------------------------------------------

# Central export defaults (publication quality)
plot_export_settings <- list(
  dpi = 300, # default DPI
  width = 8, # inches
  height = 6, # inches
  units = "in" # ggsave units
)

# Simple, consistent theme for all ggplots
theme_publication <- function(
  base_size = 14,
  base_family = ""
) {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 12, color = "#666666"),
      axis.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text = ggplot2::element_text(size = 11),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 11),
      panel.grid.major = ggplot2::element_line(color = "#E0E0E0"),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#5bc0de"),
      strip.text = ggplot2::element_text(size = 13, face = "bold", color = "white")
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
  format = c("png", "pdf", "tiff", "svg"),
  width = plot_export_settings$width,
  height = plot_export_settings$height,
  dpi = plot_export_settings$dpi
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
    content = function(file) {
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
  format = c("csv", "xlsx")
) {
  format <- match.arg(format)
  shiny::downloadHandler(
    filename = function() paste0(base_filename, "_", Sys.Date(), ".", format),
    content = function(file) {
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
  if (is.null(emm)) {
    return(NULL)
  }
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
  if (length(present) == 0L) {
    return(rep(default, nrow(df)))
  }
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
  if (is.null(x)) {
    return(FALSE)
  }
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
    if (!is.finite(y)) {
      return(NA_real_)
    }
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
  if (length(v) == 0L) {
    return(NA_integer_)
  }
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
  titlePanel("Data analysis"),
  actionButton("toggleSidebar", "Toggle Sidebar"),
  fluidRow(
    column(
      width = 3, id = "sidebar",
      h4("Data Mode:"),
      radioButtons("active_dataset", NULL,
        choices = c("Assay Data" = "assay", "Physical Endpoint Data" = "physical"),
        selected = "assay"
      ),

      # Assay controls
      conditionalPanel(
        condition = "input.active_dataset == 'assay'",
        h4("Processing mode:"),
        radioButtons("mode", NULL,
          choices = c(
            "Baseline (filter only)" = "baseline",
            "Well-level CV (wrangle_v3_avg)" = "wrangle",
            "Inter-assay CV (inter_assay_CV)" = "inter"
          ),
          selected = "baseline"
        ),
        br(),
        selectizeInput("assay", "Assay Type:", choices = endpoint_choices_assay, multiple = TRUE),
        selectizeInput("fiber", "Fiber Type:", choices = c("cotton", "pet"), multiple = TRUE),
        selectizeInput("sample_type", "Sample Type:",
          choices = sample_types_assay, # From global.R: gills, gland, hemolymph
          multiple = TRUE
        ),
        selectizeInput("week", "Week:", choices = week_choices_assay, multiple = TRUE),
        selectizeInput("treatment", "Treatment:", choices = c("Treated", "Untreated", "Control"), multiple = TRUE),
        selectizeInput("fiber_concentration", "Fiber Concentration:", choices = c("0", "100", "1000", "10000"), multiple = TRUE),
        conditionalPanel(
          condition = "input.mode != 'inter'",
          selectizeInput("plate_replicate", "Plate Replicate:", choices = NULL, multiple = TRUE),
          sliderInput("cv_threshold", "CV Threshold (%)", min = 0, max = 100, value = 15)
        ),
        selectizeInput("x_axis_assay", "X-axis:",
          choices = c(
            "Fiber Group" = "fiber_group", "Week" = "week", "Sample Type" = "sample_type",
            "Treatment" = "treatment", "Fiber Concentration" = "fiber_concentration"
          ),
          selected = "fiber_group"
        ),
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
        selectizeInput("fiber_concentration_phys", "Fiber Concentration", choices = c("0", "100", "1000", "10000"), multiple = TRUE),
        selectizeInput("treatment_phys", "Treatment", choices = c("Control", "Untreated", "Treated"), multiple = TRUE),
        selectizeInput("x_axis_phys", "X Axis",
          choices = c(
            "Fiber Group" = "fiber_group", "Week" = "week", "Fiber Type" = "fiber_type",
            "Endpoint" = "endpoint", "Tissue Type" = "tissue_type"
          ), selected = "week"
        ),
        actionButton("run_phys_analysis", "Run Physical Analysis"),
        downloadButton("download_phys_data", "Download Physical Data")
      )
    ),
    column(
      9,
      tabsetPanel(
        id = "analysis_tabs",

        # Visualization Tab
        tabPanel(
          "Visualization",
          h4("Results"), DT::dataTableOutput("table"), br(),
          tags$div(
            style = "width: 100%; max-width: 1200px; aspect-ratio: 4/3;",
            plotOutput("dist_plot", width = "100%", height = "100%")
          ), br(),
          tags$div(
            style = "width: 100%; max-width: 1200px; aspect-ratio: 4/3;",
            plotOutput("dot_plot", width = "100%", height = "100%")
          ), br(),
          # --- Visualization: Box + Points (UI) ---
          wellPanel(
            fluidRow(
              column(
                3,
                selectInput(
                  inputId = "viz_y_var",
                  label   = "Outcome (Y axis)",
                  choices = c("mean_activity"),
                  selected = "mean_activity"
                )
              ),
              column(
                3,
                selectInput(
                  inputId = "viz_facet_by",
                  label   = "Split X-axis by",
                  choices = c("None" = "none", "Week" = "week", "Fiber group" = "fiber_type"),
                  selected = "none"
                )
              ),
              column(
                3,
                checkboxInput(
                  inputId = "viz_aggregate_by_tank",
                  label   = "Aggregate to tank means",
                  value   = FALSE
                )
              ),
              column(
                3,
                helpText("Boxes: fiber concentration; shapes: tank replicate (1/3, 2/3, 3/3).")
              )
            )
          ),
          plotOutput("viz_boxpoints", height = "520px"),
          # --- Export controls under the plot ---
          div(
            style = "margin-top: 10px;",
            fluidRow(
              column(
                3,
                numericInput(
                  inputId = "viz_export_dpi",
                  label   = "PNG DPI",
                  value   = 600, min = 150, max = 1200, step = 50
                )
              ),
              column(
                3,
                numericInput(
                  inputId = "viz_export_width",
                  label   = "Width (in)",
                  value   = 10, min = 4, step = 0.5
                )
              ),
              column(
                3,
                numericInput(
                  inputId = "viz_export_height",
                  label   = "Height (in)",
                  value   = 6, min = 3, step = 0.5
                )
              ),
              column(
                3,
                div(style = "margin-top: 25px;",
                    downloadButton(
                      outputId = "download_viz_boxpoints_png",
                      label    = "Download PNG"
                    )
                )
              )
            )
          ),
          br(),
          verbatimTextOutput("info")
        ),

        # Regression Analysis Tab
        tabPanel(
          "Regression Analysis",
          sidebarLayout(
            sidebarPanel(
              width = 4,
              h4("Multiple Linear Regression"),
              wellPanel(
                h5("Model Settings"),
                selectInput("regression_endpoint", "Select Endpoint:", choices = NULL),
                selectInput("regression_dataset", "Dataset:",
                  choices = c("Assay Data" = "assay", "Physical Data" = "physical"), selected = "assay"
                ),
                checkboxGroupInput("weeks_include", "Include Weeks:",
                  choices = c("1" = "1", "3" = "3", "5" = "5"),
                  selected = c("1", "3", "5"), inline = TRUE
                ),
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
                    label = "Sample types:",
                    choices = sample_types_assay, # or sample_types_assay from global.R
                    selected = c("gills"),
                    inline = TRUE
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
                  choices = list(
                    "Treatment Groups" = "treatment", "Fiber Type" = "fiber_type",
                    "Week" = "week", "Treatment × Week" = "treatment_week",
                    "Treatment × Fiber Type" = "treatment_fiber",
                    "Fiber Type × Week" = "fiber_week",
                    "Treatment × Fiber × Week" = "all_interactions"
                  ),
                  selected = "treatment"
                ),
                selectInput("emmeans_dose", "Evaluate at dose:",
                  choices = c(
                    "Controls (0 mf/L)" = "0", "Low dose (100 mf/L)" = "2",
                    "Medium dose (1000 mf/L)" = "3", "High dose (10000 mf/L)" = "4"
                  ),
                  selected = "0"
                ),
                selectInput("comparison_type", "Comparison Type:",
                  choices = c(
                    "All Pairwise" = "pairwise",
                    "Treatment vs Control" = "control"
                  ),
                  selected = "pairwise"
                ),
                selectInput("adjustment_method", "P-value Adjustment:",
                  choices = c(
                    "Tukey HSD" = "tukey", "Holm" = "holm",
                    "BH (FDR)" = "fdr", "Bonferroni" = "bonferroni"
                  ), selected = "tukey"
                ),
                actionButton("run_emmeans", "Calculate EM Means", icon = icon("chart-line"))
              )
            ),
            mainPanel(
              width = 8,
              h5("Model Summary"),
              verbatimTextOutput("regression_model_summary"),
              h5("ANOVA Table"),
              DT::dataTableOutput("regression_anova"),
              h5("Estimated Marginal Means"),
              plotOutput("emmeans_plot", height = "420px"),
              DT::dataTableOutput("emmeans_table"),
              h5("Pairwise Comparisons"),
              DT::dataTableOutput("pairwise_table"),
              # --- Diagnostics (lm) ---
              h5("Diagnostics"),
              tabsetPanel(
                id = "regression_diag_tabs",
                tabPanel("Residuals vs Fitted", plotOutput("regression_resid_fitted", height = "360px")),
                tabPanel("Normal Q-Q",          plotOutput("regression_qq",           height = "360px")),
                tabPanel("Scale-Location",      plotOutput("regression_scale_location", height = "360px")),
                tabPanel("Residuals vs Leverage", plotOutput("regression_leverage",   height = "360px")),
                tabPanel("Fit Stats",           DT::dataTableOutput("regression_fit_stats"))
              ),
              br(),
              # --- Optional adjustments driven by diagnostics ---
              wellPanel(
                h5("Optional: Refit with adjustments"),
                fluidRow(
                  column(6,
                         checkboxInput("reg_use_spline_dose", "Use spline for dose (ns(dose_log10, 2))", value = FALSE),
                         checkboxInput("reg_use_dose_factor", "Treat dose as factor (dose_factor)", value = FALSE),
                         checkboxInput("reg_log_outcome",     "Log-transform outcome (if > 0)", value = FALSE)
                  ),
                  column(6,
                         checkboxInput("reg_use_wls",   "Weighted least squares (by dose_factor variance)", value = FALSE),
                         checkboxInput("reg_use_hc3",   "Robust HC3 standard errors", value = TRUE),
                         helpText("Select only what the diagnostics suggest; then click Refit.")
                  )
                ),
                actionButton("run_regression_adjusted", "Refit Adjusted Model", class = "btn-primary")
              ),
              h5("Adjusted Model Summary"),
              verbatimTextOutput("regression_model_summary_adjusted"),
              # --- Adjusted Diagnostics (lm) ---
              h5("Adjusted Diagnostics"),
              tabsetPanel(
                id = "regression_adjusted_diag_tabs",
                tabPanel("Residuals vs Fitted",   plotOutput("regression_adjusted_resid_fitted",   height = "340px")),
                tabPanel("Normal Q-Q",            plotOutput("regression_adjusted_qq",             height = "340px")),
                tabPanel("Scale-Location",        plotOutput("regression_adjusted_scale_location", height = "340px")),
                tabPanel("Residuals vs Leverage", plotOutput("regression_adjusted_leverage",       height = "340px")),
                tabPanel("Fit Comparison",        DT::dataTableOutput("regression_adjusted_fit_stats")),
                tabPanel("Adjusted Coefficients", DT::dataTableOutput("regression_adjusted_coef_table"))
              ),

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
                  downloadButton("download_regression_trend_png", "PNG", class = "btn-sm btn-primary", icon = icon("image")),
                  downloadButton("download_regression_trend_pdf", "PDF", class = "btn-sm btn-primary", icon = icon("file-pdf")),
                  downloadButton("download_regression_trend_tiff", "TIFF", class = "btn-sm btn-primary", icon = icon("image")),
                  downloadButton("download_regression_trend_svg", "SVG", class = "btn-sm btn-primary", icon = icon("file-image"))
                )
              )
            )
          )
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
                  choices = c("Loading..." = "loading") # updated in server
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
                    label = "Sample types:",
                    choices = sample_types_assay, # or sample_types_assay
                    selected = c("gills", "gland", "hemolymph"),
                    inline = TRUE
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
                    "Treated Cotton 100 mfL" = "Treated Cotton 100 mfL"
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
                  open = NA, # set to TRUE to start expanded

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
                downloadButton("download_combined_emm_table_csv", "CSV", class = "btn-sm", icon = icon("file-csv")),
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
                          "All Levels" = "all",
                          "Fiber Type" = "fiber",
                          "Treatment" = "treat",
                          "Concentration" = "conc"
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
                          selected = c("Cotton", "PET"), inline = TRUE
                        )
                      ),
                      conditionalPanel(
                        condition = "input.combined_emm_filter_type == 'treat'",
                        checkboxGroupInput(
                          "combined_emm_treat", "Treatments:",
                          choices = c("Control", "Untreated", "Treated"),
                          selected = c("Control", "Untreated", "Treated"), inline = TRUE
                        )
                      ),
                      conditionalPanel(
                        condition = "input.combined_emm_filter_type == 'conc'",
                        checkboxGroupInput(
                          "combined_emm_conc", "Concentrations (mf/L):",
                          choices = c("0", "100", "1000", "10000"),
                          selected = c("0", "100", "1000", "10000"), inline = TRUE
                        )
                      )
                    )
                  )
                ),
                plotOutput("combined_emm_plot", height = "600px"),
                downloadButton("download_combined_emm_png", "PNG", class = "btn-sm btn-success", icon = icon("image")),
                downloadButton("download_combined_emm_pdf", "PDF", class = "btn-sm btn-success", icon = icon("file-pdf")),
                downloadButton("download_combined_emm_tiff", "TIFF", class = "btn-sm btn-success", icon = icon("image")),
                downloadButton("download_combined_emm_svg", "SVG", class = "btn-sm btn-success", icon = icon("file-image"))
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
                        "All Comparisons" = "all",
                        "vs. Reference Only" = "ref",
                        "Within Fiber vs Control" = "fiber_control",
                        "Within Fiber Type" = "fiber",
                        "Within Treatment" = "treat",
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
                        choices = c("Cotton", "PET"), selected = "Cotton"
                      )
                    ),
                    conditionalPanel(
                      condition = "input.combined_comparison_filter == 'treat'",
                      selectInput(
                        "combined_filter_treat", "Treatment:",
                        choices = c("Control", "Untreated", "Treated"), selected = "Untreated"
                      )
                    ),
                    conditionalPanel(
                      condition = "input.combined_comparison_filter == 'conc'",
                      selectInput(
                        "combined_filter_conc", "Concentration (mf/L):",
                        choices = c("0", "100", "1000", "10000"), selected = "0"
                      )
                    )
                  )
                ),
                DT::dataTableOutput("combined_pairwise"),
                downloadButton("download_combined_pairwise_csv", "CSV", class = "btn-sm", icon = icon("file-csv")),
                downloadButton("download_combined_pairwise_xlsx", "XLSX", class = "btn-sm", icon = icon("file-excel")),
                br(),
                plotOutput("combined_pairwise_plot", height = "700px"),
                downloadButton("download_combined_pairwise_png", "PNG", class = "btn-sm btn-warning", icon = icon("image")),
                downloadButton("download_combined_pairwise_pdf", "PDF", class = "btn-sm btn-warning", icon = icon("file-pdf")),
                downloadButton("download_combined_pairwise_tiff", "TIFF", class = "btn-sm btn-warning", icon = icon("image")),
                downloadButton("download_combined_pairwise_svg", "SVG", class = "btn-sm btn-warning", icon = icon("file-image"))
              ),
              # --- Diagnostics (Combined Treatment model) ---
              tabPanel(
                "Diagnostics",
                br(),
                tabsetPanel(
                  id = "combined_diag_tabs",
                  tabPanel("Residuals vs Fitted",   plotOutput("combined_resid_fitted",   height = "380px")),
                  tabPanel("Normal Q-Q",            plotOutput("combined_qq",             height = "380px")),
                  tabPanel("Scale-Location",        plotOutput("combined_scale_location", height = "380px")),
                  tabPanel("Residuals vs Leverage", plotOutput("combined_leverage",       height = "380px")),
                  tabPanel("Fit Stats",             DT::dataTableOutput("combined_fit_stats"))
                )
              )
            )
          )
        ),

        # Mixed Effects Analysis Tab
        tabPanel(
          "Mixed Effects Analysis",
          sidebarLayout(
            sidebarPanel(
              width = 4,
              h4("Linear Mixed Effects Regression"),
              # --- Mixed Effects: Data & filters (UI) ----------------
              wellPanel(
                h5("Data & filters"),
                fluidRow(
                  column(
                    6,
                    # Machine-friendly values; stable snake_case id
                    selectInput(
                      inputId = "lmer_dataset",
                      label = "Dataset:",
                      choices = c("Assay Data" = "assay", "Physical Data" = "physical"),
                      selected = "assay"
                    )
                  ),
                  column(
                    6,
                    # Start empty; server populates based on dataset
                    selectInput(
                      inputId = "lmer_endpoint",
                      label = "Endpoint:",
                      choices = character(0),
                      selected = NULL
                    )
                  )
                ),

                # Assay sample types OR Physical tissues (dynamic)
                conditionalPanel(
                  condition = "input.lmer_dataset == 'assay'",
                  checkboxGroupInput("lmer_sample_types", "Sample types:",
                    choices = sample_types_assay, selected = sample_types_assay, inline = TRUE
                  )
                ),
                conditionalPanel(
                  condition = "input.lmer_dataset == 'physical'",
                  uiOutput("lmer_tissues_ui") # shows choices for mf_counts; disabled “All” otherwise
                ),

                # Two short rows for core filters
                fluidRow(
                  column(
                    6,
                    checkboxGroupInput("lmer_fiber_types", "Fiber types:",
                      choices = c("cotton", "pet"), selected = c("cotton", "pet"), inline = TRUE
                    )
                  ),
                  column(
                    6,
                    checkboxGroupInput("lmer_treatments", "Treatments:",
                      choices = c("untreated", "treated"), selected = c("untreated", "treated"), inline = TRUE
                    )
                  )
                ),
                fluidRow(
                  column(
                    12,
                    checkboxGroupInput("lmer_concentrations", "Fiber concentration (mf/L):",
                      choices = c("0", "100", "1000", "10000"),
                      selected = c("0", "100", "1000", "10000"), inline = TRUE
                    ),
                    helpText("Weeks are controlled below in Include Weeks.")
                  )
                )
              ),

              # --- Mixed Effects: Random-effects structure (single section) ---
              wellPanel(
                h5("Random-effects structure"),
                uiOutput("model_random_effects_ui")
              ),

              # --- Mixed Effects: Model settings (remove the duplicate random-intercepts here) ---
              wellPanel(
                h5("Model Settings"),
                checkboxInput("lmer_include_three_way", "Include 3‑way: fiber × treatment × week", value = TRUE),
                checkboxInput("lmer_dose_by_fiber", "Include dose × fiber_type", value = FALSE),
                checkboxInput("lmer_dose_by_treat", "Include dose × chem_treatment", value = FALSE),
                checkboxInput("lmer_dose_as_factor", "Treat dose as a factor", value = FALSE),
                checkboxGroupInput("lmer_weeks_include", "Include Weeks:",
                  choices = c("1", "3", "5"),
                  selected = c("1", "3", "5"), inline = TRUE
                ),
                actionButton("run_lmer", "Run Mixed Model", icon = icon("play"))
              ),
              wellPanel(
                h5("Post-hoc Comparisons"),
                selectInput("lmer_emmeans_by", "Calculate EM Means by:",
                  choices = list(
                    "Treatment Groups" = "treatment", "Fiber Type" = "fiber_type",
                    "Week" = "week", "Treatment × Week" = "treatment_week",
                    "Treatment × Fiber Type" = "treatment_fiber",
                    "Fiber Type × Week" = "fiber_week",
                    "Treatment × Fiber × Week" = "all_interactions"
                  ), selected = "treatment"
                ),
                selectInput("lmer_emmeans_dose", "Evaluate at dose:",
                  choices = c(
                    "Controls (0 mf/L)" = "0", "Low dose (100 mf/L)" = "2",
                    "Medium dose (1000 mf/L)" = "3", "High dose (10000 mf/L)" = "4"
                  ), selected = "0"
                ),
                selectInput("lmer_comparison_type", "Comparison Type:",
                  choices = c(
                    "All Pairwise" = "pairwise",
                    "Treatment vs Control" = "control"
                  ),
                  selected = "pairwise"
                ),
                selectInput("lmer_adjustment_method", "P-value Adjustment:",
                  choices = c(
                    "Tukey HSD" = "tukey", "Holm" = "holm",
                    "BH (FDR)" = "fdr", "Bonferroni" = "bonferroni"
                  ), selected = "tukey"
                ),
                actionButton("lmer_run_emmeans", "Calculate EM Means", icon = icon("chart-line"))
              )
            ),
            mainPanel(
              width = 8,
              h5("Mixed Model Summary"), verbatimTextOutput("lmer_model_summary"),
              h5("ANOVA Table"), DT::dataTableOutput("lmer_anova"),
              h5("Estimated Marginal Means"), plotOutput("lmer_emmeans_plot", height = "420px"),
              DT::dataTableOutput("lmer_emmeans_table"),
              h5("Pairwise / Custom Contrasts"), DT::dataTableOutput("lmer_pairwise_table"),
              h5("Model Comparison (lm vs lmer)"), DT::dataTableOutput("model_compare"),
              # --- Diagnostics (Mixed Effects) ---
              h5("Diagnostics"),
              tabsetPanel(
                id = "lmer_diag_tabs",
                tabPanel(
                  "DHARMa residuals",
                  plotOutput("lmer_dharma", height = "420px")
                ),
                tabPanel(
                  "Residuals vs Fitted",
                  plotOutput("lmer_resid_fitted", height = "380px")
                ),
                tabPanel(
                  "Random effects (BLUPs)",
                  plotOutput("lmer_re_caterpillar", height = "420px")
                ),
                tabPanel(
                  "Fit Stats",
                  DT::dataTableOutput("lmer_fit_stats")
                )
              
              )
            )
          )
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
                  label = "Baseline Week:",
                  choices = c(), # server will update
                  selected = NULL
                ),
                selectInput(
                  inputId = "recovery_recovery_week",
                  label = "Recovery Week:",
                  choices = c(), # server will update
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
                    choices = c("Hemolymph", "Gills", "Gland"),
                    selected = c("Hemolymph", "Gills", "Gland")
                  ),
                  checkboxGroupInput(
                    "recovery_concentration", "Fiber concentration (mf/L):",
                    choices = c("0" = "0", "100" = "100", "1000" = "1000", "10000" = "10000"),
                    selected = c("0", "100", "1000", "10000"),
                    inline = TRUE
                  )
                ),

                # Physical: tissues + dose for mf_counts
                conditionalPanel(
                  condition = "input.recovery_dataset == 'physical'",
                  uiOutput("recovery_tissues_ui"),
                  uiOutput("recovery_tissues_hint"),
                  checkboxGroupInput(
                    "recovery_concentration_physical", "Fiber concentration (mf/L):",
                    choices = c("0" = "0", "100" = "100", "1000" = "1000", "10000" = "10000"),
                    selected = c("0", "100", "1000", "10000"),
                    inline = TRUE
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
                  choices = sort(unique(physical_master$week)), # assumes weeks present in physical
                  selected = max(sort(unique(physical_master$week))),
                  multiple = FALSE
                ),
                # Two-tissue selector (exactly 2)
                selectizeInput(
                  "transloc_tissues", "Tissues to compare:",
                  choices = tissue_types_physical,
                  selected = c("gland", "tissue"), # c("gills","gland","tissue") from global.R
                  multiple = TRUE, options = list(maxItems = 2)
                ),
                # Fiber concentration (multi)
                selectizeInput(
                  "transloc_conc", "Fiber concentration (mF/L):",
                  choices = c("0", "100", "1000", "10000"),
                  selected = c("100", "1000", "10000"),
                  multiple = TRUE
                ),
                selectInput(
                  inputId = "transloc_facet_rows",
                  label = "Facet rows:",
                  choices = c(
                    "None" = "none",
                    "Fiber type" = "fiber_type",
                    "Treatment" = "treatment"
                  ),
                  selected = "fiber_type"
                ),
                selectInput(
                  inputId = "transloc_facet_cols",
                  label = "Facet columns:",
                  choices = c(
                    "None" = "none",
                    "Fiber type" = "fiber_type",
                    "Treatment" = "treatment"
                  ),
                  selected = "treatment"
                ),
                helpText("Choose how to layout panels; pick none, one, or both dimensions."),

                # --- Translocation: facet filters (UI) ---
                checkboxGroupInput(
                  inputId  = "transloc_filter_fibers",
                  label    = "Facet filter: fiber_type",
                  choices  = c("cotton", "pet"),
                  selected = c("cotton", "pet"),
                  inline   = TRUE
                ),
                checkboxGroupInput(
                  inputId  = "transloc_filter_treatments",
                  label    = "Facet filter: treatment",
                  choices  = c("treated", "untreated"),
                  selected = c("treated", "untreated"),
                  inline   = TRUE
                ),
                helpText("Filters above restrict which fiber × treatment panels are shown; facets update automatically."),

                # --- Extra controls for the new multi-week translocation plot ---
                h5("Multi‑week plot controls"),

                # Weeks: allow multiple
                selectizeInput(
                  inputId = "transloc_weeks_multi",
                  label = "Weeks (multi):",
                  choices = sort(unique(physical_master$week)), # uses your loaded physical data
                  selected = sort(unique(physical_master$week)), # preselect all available
                  multiple = TRUE
                ),

                # Tissue: enforce a single tissue
                selectizeInput(
                  inputId = "transloc_tissue_single",
                  label = "Tissue (single):",
                  choices = tissue_types_physical, # from global.R
                  selected = head(tissue_types_physical, 1),
                  multiple = FALSE
                ),

                # Concentrations: default to 100/1000/10000
                checkboxGroupInput(
                  inputId  = "transloc_conc_multi",
                  label    = "Concentrations (mf/L):",
                  choices  = c("100", "1000", "10000"),
                  selected = c("100", "1000", "10000"),
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
                          label = "P-value adjust:",
                          choices = c(
                            "BH (FDR)" = "BH",
                            "Holm" = "holm",
                            "Bonferroni" = "bonferroni",
                            "None" = "none"
                          ),
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
                tabPanel(
                  "Translocation Plots",
                  plotOutput("transloc_counts_plot", height = 420),
                  plotOutput("transloc_diff_plot", height = 360),
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
                  column(
                    6,
                    selectInput(
                      "did_dataset", "Dataset:",
                      choices = c("Assay Data" = "assay", "Physical Data" = "physical"),
                      selected = "assay"
                    )
                  ),
                  column(
                    6,
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
                        label = "Sample types:",
                        choices = character(0), # filled by server
                        selected = NULL,
                        inline = TRUE
                      ),
                      uiOutput("did_samples_hint") # hint text managed by server
                    )
                  ),

                  # Physical selectors
                  conditionalPanel(
                    condition = "input.did_dataset == 'physical'",
                    div(
                      id = "did_tissues_box",
                      checkboxGroupInput(
                        inputId = "did_tissues",
                        label = "Tissues:",
                        choices = character(0), # filled by server
                        selected = NULL,
                        inline = TRUE
                      ),
                      uiOutput("did_tissues_hint") # hint text managed by server
                    )
                  )
                ),

                # Week selection - for across_fiber and across_treatment
                conditionalPanel(
                  condition = "input.did_mode == 'across_fiber' || input.did_mode == 'across_treatment'",
                  selectInput("did_week_select", "Select Week:",
                    choices = c("1", "3", "5"), selected = "5"
                  )
                ),

                # Fiber selection - for across_week and across_treatment
                conditionalPanel(
                  condition = "input.did_mode == 'across_week' || input.did_mode == 'across_treatment'",
                  selectInput("did_fiber_select", "Select Fiber:",
                    choices = c("Cotton" = "cotton", "PET" = "pet"),
                    selected = "cotton"
                  )
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
                    selected = "pet"
                  ),
                  selectInput("did_fiber_b", "Fiber B:",
                    choices = c("Cotton" = "cotton", "PET" = "pet"),
                    selected = "cotton"
                  )
                ),

                # Week A and B - only for across_week
                conditionalPanel(
                  condition = "input.did_mode == 'across_week'",
                  selectInput("did_week_a", "Week A:",
                    choices = c("1", "3", "5"), selected = "3"
                  ),
                  selectInput("did_week_b", "Week B:",
                    choices = c("1", "3", "5"), selected = "5"
                  )
                ),

                # Baseline selector (now includes untreated_ctrl)
                conditionalPanel(
                  condition = "input.did_mode == 'across_fiber'",
                  radioButtons(
                    inputId = "did_baseline_type",
                    label = "Baseline for DiD:",
                    choices = c(
                      "Untreated vs Treated (current)"        = "treatment",
                      "Dose vs 0 control (new)"               = "dose",
                      "Untreated dose vs 0 (control-offset)"  = "untreated_ctrl" # NEW
                    ),
                    selected = "treatment",
                    inline = FALSE
                  ),
                  # Dose picker shown for both 'dose' and 'untreated_ctrl'
                  conditionalPanel(
                    condition = "input.did_baseline_type == 'dose' || input.did_baseline_type == 'untreated_ctrl'",
                    selectInput(
                      inputId = "did_treated_dose",
                      label = "Dose for baseline (mf/L):",
                      choices = c("100", "1000", "10000"),
                      selected = "100"
                    )
                  )
                ),

                # Treatment A and B - only for across_treatment
                conditionalPanel(
                  condition = "input.did_mode == 'across_treatment'",
                  selectInput("did_trt_a", "Treatment A:",
                    choices = c("Treated" = "treated", "Untreated" = "untreated"),
                    selected = "treated"
                  ),
                  selectInput("did_trt_b", "Treatment B:",
                    choices = c("Treated" = "treated", "Untreated" = "untreated"),
                    selected = "untreated"
                  )
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
  observeEvent(input$toggleSidebar, {
    toggle("sidebar")
  })

  # Update plate replicate choices from pre-loaded data
  observe({
    updateSelectizeInput(session, "plate_replicate",
      choices = sort(unique(final_data$plate_replicate))
    )
  })

  # Dynamic tissue filter UI for mf_counts endpoint
  output$tissue_type_ui <- renderUI({
    # Only render in the Physical Endpoint Data mode
    req(input$active_dataset == "physical") # FIX: correct id

    # Only show when mf_counts is among the selected endpoints
    if (!is.null(input$endpoint) && "mf_counts" %in% input$endpoint) {
      selectizeInput(
        inputId = "tissue_type_phys",
        label = "Tissue Type",
        choices = tissue_types_physical, # from global.R
        multiple = TRUE,
        selected = NULL
      )
    } else {
      NULL # Hide the input for non-mf_counts endpoints
    }
  })
  
  # =============================================================================
  # Visualization: Box + Points (server logic)
  # =============================================================================
  
  # --- Data source: follow sidebar mode and filters ---
  get_visualization_data <- function() {
    if (identical(input$active_dataset, "assay")) {
      if (identical(input$mode, "baseline")) {
        df <- assay_filtered()
      } else {
        df <- summarized()
      }
    } else {
      if (identical(input$mode_phys, "baseline")) {
        df <- phys_filtered()
      } else {
        df <- phys_summarized()
      }
    }
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
    df
  }
  
  # Ensure we have consistent concentration fields and a 1/3, 2/3, 3/3 tank replicate label
  ensure_viz_concentration <- function(df) {
    # Pull numeric concentration from an existing concentration_factor if present
    conc_from_cf <- NULL
    if ("concentration_factor" %in% names(df)) {
      cf <- as.character(df[["concentration_factor"]])
      conc_from_cf <- suppressWarnings(as.numeric(cf))
      if (all(is.na(conc_from_cf))) conc_from_cf <- NULL
    }
    
    # Fallbacks for numeric concentration
    get_conc_num <- function(d) {
      cands <- list()
      if ("fiber_concentration" %in% names(d)) cands[[length(cands)+1]] <- suppressWarnings(as.numeric(d[["fiber_concentration"]]))
      if ("fiber_concentration_numeric" %in% names(d)) cands[[length(cands)+1]] <- suppressWarnings(as.numeric(d[["fiber_concentration_numeric"]]))
      if ("concentration" %in% names(d)) cands[[length(cands)+1]] <- suppressWarnings(as.numeric(d[["concentration"]]))
      if ("dose_factor" %in% names(d)) cands[[length(cands)+1]] <- suppressWarnings(as.numeric(as.character(d[["dose_factor"]])))
      if (length(cands) == 0) return(rep(0, nrow(d)))
      v <- Reduce(dplyr::coalesce, cands)
      v[is.na(v)] <- 0
      v
    }
    
    conc_num <- conc_from_cf %||% get_conc_num(df)
    
    # Build canonical factor and logs, but DO NOT overwrite any existing tank column
    df$dose_log10 <- ifelse(conc_num > 0, log10(conc_num), 0)
    df$dose_factor <- factor(as.character(as.integer(conc_num)))
    df$concentration_factor <- forcats::fct_inorder(factor(as.character(as.integer(conc_num))))
    
    # Preserve original tank ids if present; otherwise synthesize stable ids
    if (!("tank" %in% names(df))) {
      df$tank <- factor(paste0("tank_", seq_len(nrow(df))))
    } else {
      df$tank <- as.factor(df$tank)
    }
    
    # Assign 1/2/3 replicate labels within each concentration based on sorted tank ids
    df <- df %>%
      dplyr::group_by(concentration_factor) %>%
      dplyr::mutate(
        tank_rank_in_group = dplyr::dense_rank(as.character(tank)),
        tank_shape_id_num  = ((tank_rank_in_group - 1L) %% 3L) + 1L,
        tank_replicate_of_3 = factor(tank_shape_id_num, levels = c(1, 2, 3))
      ) %>%
      dplyr::ungroup()
    
    df
  }
  
  # --- Auto-update Y-axis choices when dataset/mode changes ---
  observeEvent(
    list(input$active_dataset, input$mode, input$mode_phys, input$run_analysis, input$run_phys_analysis),
    ignoreInit = FALSE, {
      df <- get_visualization_data()
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
        updateSelectInput(session, "viz_y_var", choices = "mean_activity", selected = "mean_activity")
        return(invisible(NULL))
      }
      
      preferred <- c("mean_activity", "mean_value", "calculated_concentration", "value_corr", "value", "outcome")
      preferred_present <- preferred[preferred %in% names(df)]
      other_numeric <- setdiff(names(df)[vapply(df, is.numeric, logical(1))], preferred_present)
      y_choices <- unique(c(preferred_present, other_numeric))
      if (length(y_choices) == 0) y_choices <- "mean_activity"
      
      sel <- if (!is.null(input$viz_y_var) && input$viz_y_var %in% y_choices) {
        input$viz_y_var
      } else {
        y_choices[[1]]
      }
      
      updateSelectInput(session, "viz_y_var", choices = y_choices, selected = sel)
    }
  )
  
  # --- Builder that returns the ggplot object used by the box+points viz ---
  build_viz_boxpoints_plot <- function() {
    df   <- viz_box_df()                       # existing reactive supplying the data
    yvar <- input$viz_y_var %||% "mean_activity"
    
    validate(need(nrow(df) > 0, "No rows to plot."))
    validate(need(any(is.finite(df[[yvar]])), "All values are NA/Inf for the selected outcome."))
    
    # Decide x mapping (week / fiber_group / concentration), same as in renderPlot
    facet_var <- input$viz_facet_by %||% "none"
    if (facet_var == "week" && "week" %in% names(df)) {
      df$week <- forcats::fct_inorder(as.factor(df$week))
      x_var <- "week"; x_lab <- "Week"
    } else if (facet_var == "fiber_type") {
      if ("fiber_group" %in% names(df)) {
        df$fiber_group <- forcats::fct_inorder(as.factor(df$fiber_group))
        x_var <- "fiber_group"; x_lab <- "Fiber group"
      } else {
        df$fiber_type <- forcats::fct_inorder(as.factor(df$fiber_type))
        x_var <- "fiber_type"; x_lab <- "Fiber type"
      }
    } else {
      x_var <- "concentration_factor"; x_lab <- "Fiber concentration"
    }
    
    # Colors for concentrations (same as on-screen)
    conc_levels <- levels(df$concentration_factor)
    bright_pal  <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                     "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
    conc_colors <- stats::setNames(rep(bright_pal, length.out = length(conc_levels)), conc_levels)
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = .data[[yvar]])) +
      ggplot2::geom_boxplot(
        ggplot2::aes(fill = concentration_factor),
        width = 0.65, outlier.shape = NA, alpha = 0.7,
        position = ggplot2::position_dodge(width = 0.75)
      ) +
      # Make shapes smaller and outlines thinner
      ggplot2::geom_point(
        ggplot2::aes(
          color = concentration_factor,
          shape = tank_replicate_of_3,
          group = concentration_factor
        ),
        size = 2.0,      # was 2.8
        stroke = 0.9,    # was 1.4
        fill = NA,
        position = ggplot2::position_jitterdodge(jitter.width = 0, dodge.width = 0.75)
      ) +
      ggplot2::scale_shape_manual(
        name = "Tank replicate",
        values = c("1" = 21, "2" = 22, "3" = 24),
        breaks = c("1","2","3"),
        labels = c("1/3","2/3","3/3"),
        drop = FALSE
      ) +
      ggplot2::scale_fill_manual(name = "Concentration", values = conc_colors, drop = FALSE) +
      ggplot2::scale_color_manual(name = "Concentration", values = conc_colors, drop = FALSE) +
      ggplot2::labs(title = paste0(yvar, " by ", x_lab), x = x_lab, y = yvar) +
      theme_publication()
    
    p
  }
  
  # Build a safe file-name token from a character vector
  make_safe_tag <- function(values, max_show = 3) {
    vals <- unique(na.omit(trimws(as.character(values))))
    if (length(vals) == 0) vals <- "all"
    if (length(vals) > max_show) vals <- c("mixed")
    tag  <- paste(vals, collapse = "_")
    # slugify: lowercase, replace non-alnum with '-', trim dashes
    tag  <- tolower(tag)
    tag  <- gsub("[^a-z0-9]+", "-", tag)
    tag  <- gsub("^-|-$", "", tag)
    tag
  }
  
  # Build the three tags from the currently filtered dataset
  build_filter_tags <- function(df) {
    fiber_tag    <- if ("fiber_type"  %in% names(df)) make_safe_tag(df$fiber_type)  else "all"
    treat_tag    <- if ("treatment"   %in% names(df)) make_safe_tag(df$treatment)   else "all"
    sample_tag   <- if ("sample_type" %in% names(df)) make_safe_tag(df$sample_type) else "all"
    list(fiber = fiber_tag, treat = treat_tag, sample = sample_tag)
  }

  # ---------------------------
  # Regression: tissues UI (physical)
  # ---------------------------
  output$reg_tissues_ui <- renderUI({
    req(input$regression_dataset == "physical")
    ep <- input$regression_endpoint %||% ""
    if (identical(ep, "mf_counts")) {
      checkboxGroupInput(
        inputId = "reg_tissues",
        label = "Tissues:",
        choices = tissue_types_physical, # from global.R, e.g., c("gills","gland","tissue")
        selected = tissue_types_physical,
        inline = TRUE
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
        label = "Tissues:",
        choices = tissue_types_physical,
        selected = tissue_types_physical,
        inline = TRUE
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

  # --- Helper: dynamic facetting for ggplot ---
  apply_dynamic_facets <- function(p, rows = "fiber_type", cols = "treatment") {
    # Prevent duplicated mapping (e.g., rows == cols) by collapsing to a single axis
    if (!is.null(rows) && !is.null(cols) && rows == cols && rows != "none") {
      cols <- "none"
    }
    # Build facet formula via base::reformulate to avoid NSE
    if (identical(rows, "none") && identical(cols, "none")) {
      return(p)
    } else if (!identical(rows, "none") && !identical(cols, "none")) {
      ff <- reformulate(termlabels = cols, response = rows) # rows ~ cols
      return(p + ggplot2::facet_grid(ff, labeller = ggplot2::label_both, drop = TRUE))
    } else if (!identical(rows, "none")) {
      ff <- reformulate(termlabels = ".", response = rows) # rows ~ .
      return(p + ggplot2::facet_grid(ff, labeller = ggplot2::label_both, drop = TRUE))
    } else {
      ff <- reformulate(termlabels = cols, response = ".") # . ~ cols
      return(p + ggplot2::facet_grid(ff, labeller = ggplot2::label_both, drop = TRUE))
    }
  }

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
        available_samples <- sample_types_assay # Fallback to all allowed types
      }

      # Update the choices dynamically
      updateSelectizeInput(session, "sample_type",
        choices = available_samples,
        selected = NULL
      ) # Select all by default
    }
  })

  # ---- Formula utilities (patched) ----
  
  # Collapse multi-line formula objects to a single-line formula
  collapse_formula <- function(fml) {
    # Use base::deparse; stats::deparse is not exported
    f_chr <- paste(base::deparse(fml), collapse = " ")
    stats::as.formula(f_chr)
  }
  
  # Safely replace a single symbol in a formula by editing the collapsed string
  replace_term_in_formula <- function(fml, from_sym, to_expr) {
    f_chr <- paste(base::deparse(fml), collapse = " ")
    pattern <- paste0("(?<![A-Za-z0-9_\\.])", from_sym, "(?![A-Za-z0-9_\\.])")
    if (!grepl(pattern, f_chr, perl = TRUE)) return(fml)  # no-op if term absent
    f_chr2 <- gsub(pattern, to_expr, f_chr, perl = TRUE)
    stats::as.formula(f_chr2)
  }
  
  # ---- Helper to get refit data with grouping present ----
  # Tries to use your full regression data if available; otherwise falls back to the model frame
  get_regression_refit_data <- function(base_fit) {
    # Prefer a cached data function if it exists in your app
    df_full <- tryCatch(regression_data_cached(), error = function(e) NULL)
    if (!is.null(df_full) && is.data.frame(df_full)) return(df_full)
    # Fallback: base model frame (limited to variables in the original formula)
    return(base_fit$model)
  }
  
  # --- Helpers: ensure dose columns exist on the refit data ---------------------
  ensure_dose_columns <- function(df) {
    # Prefer existing numeric concentration fields
    get_conc_num <- function(d) {
      cand <- NULL
      if ("fiber_concentration" %in% names(d)) cand <- suppressWarnings(as.numeric(d[["fiber_concentration"]]))
      if (is.null(cand) && "fiber_concentration_numeric" %in% names(d)) cand <- suppressWarnings(as.numeric(d[["fiber_concentration_numeric"]]))
      if (is.null(cand) && "concentration" %in% names(d)) cand <- suppressWarnings(as.numeric(d[["concentration"]]))
      if (is.null(cand) && "dose_factor" %in% names(d)) cand <- suppressWarnings(as.numeric(as.character(d[["dose_factor"]])))
      if (is.null(cand)) cand <- rep(NA_real_, nrow(d))
      cand[is.na(cand)] <- 0
      cand
    }
    
    # Build dose_log10 if absent: control (0) -> 0 on log10 scale
    if (!("dose_log10" %in% names(df))) {
      conc <- get_conc_num(df)
      df$is_control <- as.integer(conc == 0)
      df$dose_log10 <- ifelse(conc > 0, log10(conc), 0)
    }
    
    # Build dose_factor if absent from same numeric concentration
    if (!("dose_factor" %in% names(df))) {
      conc <- get_conc_num(df)
      df$dose_factor <- factor(as.character(as.integer(conc)))
    }
    
    df
  }
  
  # Safe package check
  has_pkg <- function(pkg) {
    isTRUE(requireNamespace(pkg, quietly = TRUE))
  }
  
  # Build a compact diagnostics list for an lmerMod
  build_lmer_diagnostics <- function(fit, df_model) {
    stopifnot(inherits(fit, "lmerMod"))
    fitted_vals <- stats::fitted(fit)
    resid_raw   <- stats::residuals(fit)            # raw residuals
    sigma_hat   <- sigma(fit)
    resid_std   <- as.numeric(resid_raw) / sigma_hat
    
    # Put key covariates (if present) for stratified plots
    diag_df <- data.frame(
      fitted   = as.numeric(fitted_vals),
      resid    = as.numeric(resid_raw),
      resid_std = resid_std,
      fiber_type      = if ("fiber_type" %in% names(df_model)) df_model$fiber_type else NA_character_,
      chem_treatment  = if ("chem_treatment" %in% names(df_model)) df_model$chem_treatment else NA_character_,
      week            = if ("week" %in% names(df_model)) as.character(df_model$week) else NA_character_,
      dose_log10      = if ("dose_log10" %in% names(df_model)) df_model$dose_log10 else NA_real_,
      stringsAsFactors = FALSE
    )
    
    # Random effects (per-group) for QQ plots
    re_list <- lme4::ranef(fit, condVar = TRUE)
    # Take first grouping (if available) for a QQ plot; code below allows specific plotting per group
    ranef_bydim <- lapply(re_list, function(re) {
      data.frame(effect = as.numeric(re[[1L]]), stringsAsFactors = FALSE)
    })
    
    # Model quality: R2 (if performance available), collinearity, convergence, singularity
    r2_vals <- NULL
    if (has_pkg("performance")) {
      r2_vals <- tryCatch(
        performance::r2_nakagawa(fit), error = function(e) NULL
      )
    }
    
    collin_tbl <- NULL
    if (has_pkg("performance")) {
      collin_tbl <- tryCatch(as.data.frame(performance::check_collinearity(fit)),
                             error = function(e) NULL)
    }
    
    conv_msg <- tryCatch({
      c(
        singular = lme4::isSingular(fit, tol = 1e-6),
        conv    = !is.null(fit@optinfo) && is.null(fit@optinfo$conv$lme4$messages)
      )
    }, error = function(e) c(singular = NA, conv = NA))
    
    list(
      diag_df    = diag_df,
      ranef_list = ranef_bydim,
      aic        = AIC(fit),
      bic        = BIC(fit),
      sigma_hat  = sigma_hat,
      r2         = r2_vals,
      collinearity = collin_tbl,
      singular   = conv_msg[["singular"]],
      converged  = conv_msg[["conv"]]
    )
  }
  
  # ============================================================================
  # HELPER: Normalize categories for Mixed Effects models
  # ============================================================================

  # -----------------------------------------------------------------------------
  # Mixed Effects: input ID compatibility shim (maps underscored and non-underscored)
  # -----------------------------------------------------------------------------
  get_input_value <- function(input, ids, default = NULL) {
    for (id in ids) {
      val <- input[[id]]
      if (!is.null(val)) {
        return(val)
      }
    }
    default
  }
  
  # ---- helpers (place once near other small utilities) ----
  col_or_na <- function(d, nm, type = c("chr","fct")) {
    type <- match.arg(type)
    if (nm %in% names(d)) {
      v <- d[[nm]]
    } else {
      v <- rep(NA_character_, nrow(d))
    }
    if (type == "chr") as.character(v) else factor(as.character(v))
  }

  # Inspect grouping structure in the data used for lmer
  inspect_grouping <- function(df) {
    n_tank  <- if ("tank" %in% names(df)) dplyr::n_distinct(df$tank) else 0L
    n_batch <- if ("experiment_batch" %in% names(df)) dplyr::n_distinct(df$experiment_batch) else 0L
    n_fiber <- if ("fiber_type" %in% names(df)) dplyr::n_distinct(df$fiber_type) else 0L
    
    # Detect perfect confounding: each batch corresponds to exactly one fiber
    confounded_batch_fiber <- ("experiment_batch" %in% names(df)) &&
      ("fiber_type" %in% names(df)) &&
      dplyr::n_distinct(df$experiment_batch) == dplyr::n_distinct(df$fiber_type) &&
      dplyr::n_distinct(interaction(df$experiment_batch, df$fiber_type, drop = TRUE)) ==
      dplyr::n_distinct(df$experiment_batch)
    
    list(
      n_tank  = n_tank,
      n_batch = n_batch,
      n_fiber = n_fiber,
      confounded_batch_fiber = confounded_batch_fiber,
      cross_tab = df %>%
        dplyr::count(experiment_batch, fiber_type, name = "n_rows") %>%
        dplyr::arrange(experiment_batch, fiber_type)
    )
  }
  
  # --- Mixed Effects: endpoint selection that preserves user choice -------------

  # Remember the user’s last explicit endpoint
  last_endpoint_sel <- reactiveVal(NULL)

  # Track user-initiated changes to the endpoint control
  observeEvent(input$lmer_endpoint, ignoreInit = TRUE, {
    last_endpoint_sel(input$lmer_endpoint)
  })

  # Helper to pick the correct data frame and endpoint column
  get_endpoint_source <- function(dataset_tag) {
    ds <- tolower(as.character(dataset_tag %||% "assay"))
    is_assay <- grepl("assay", ds)
    is_physical <- grepl("phys", ds)

    if (is_assay && exists("final_data") && is.data.frame(final_data) && nrow(final_data) > 0) {
      return(list(src = final_data, ep_col = "assay_type"))
    }
    if (is_physical && exists("physical_master") && is.data.frame(physical_master) && nrow(physical_master) > 0) {
      return(list(src = physical_master, ep_col = "endpoint"))
    }
    # Sensible fallbacks if load order varies
    if (exists("final_data") && "assay_type" %in% names(final_data) && nrow(final_data) > 0) {
      return(list(src = final_data, ep_col = "assay_type"))
    }
    if (exists("physical_master") && "endpoint" %in% names(physical_master) && nrow(physical_master) > 0) {
      return(list(src = physical_master, ep_col = "endpoint"))
    }
    list(src = NULL, ep_col = NULL)
  }

  # Update choices ONLY when dataset changes; preserve prior valid selection
  observeEvent(input$lmer_dataset, ignoreInit = TRUE, priority = 100, {
    src_info <- get_endpoint_source(input$lmer_dataset)

    ep_choices <- if (!is.null(src_info$src) && !is.null(src_info$ep_col)) {
      sort(unique(stats::na.omit(as.character(src_info$src[[src_info$ep_col]]))))
    } else {
      character(0)
    }

    prev <- isolate(last_endpoint_sel()) %||% isolate(input$lmer_endpoint)
    selected <- if (length(ep_choices) && !is.null(prev) && prev %in% ep_choices) {
      prev
    } else if (length(ep_choices)) {
      ep_choices[1]
    } else {
      NULL
    }

    # IMPORTANT: freeze input, not session
    shiny::freezeReactiveValue(input, "lmer_endpoint")
    updateSelectInput(session, "lmer_endpoint", choices = ep_choices, selected = selected)
  })

  # Initialize once after UI is ready if data are already loaded
  session$onFlushed(function() {
    isolate({
      if (is.null(input$lmer_endpoint) || identical(input$lmer_endpoint, "")) {
        src_info <- get_endpoint_source(input$lmer_dataset %||% "assay")
        ep_choices <- if (!is.null(src_info$src) && !is.null(src_info$ep_col)) {
          sort(unique(stats::na.omit(as.character(src_info$src[[src_info$ep_col]]))))
        } else {
          character(0)
        }
        if (length(ep_choices)) {
          shiny::freezeReactiveValue(input, "lmer_endpoint")
          updateSelectInput(session, "lmer_endpoint", choices = ep_choices, selected = ep_choices[1])
        }
      }
    })
  }, once = TRUE)

  build_random_terms <- function(df, include_time_slope = TRUE) {
    validate(need("tank" %in% names(df), "Random-effects require a tank column."))
    n_tanks <- dplyr::n_distinct(df$tank)
    validate(need(n_tanks >= 2L, "Random-effects require at least 2 tanks in the filtered data."))
    
    slope_ok <- isTRUE(include_time_slope) &&
      ("time_wk_z" %in% names(df)) &&
      !all(is.na(df$time_wk_z)) &&
      isTRUE(stats::sd(df$time_wk_z, na.rm = TRUE) > 0)
    
    if (!slope_ok && isTRUE(include_time_slope)) {
      message("[lmer] Tank-level week slope requested but time_wk_z is missing/constant; using intercept only.")
    }
    
    if (slope_ok) {
      "(1 + time_wk_z | tank)"
    } else {
      "(1 | tank)"
    }
  }
  
  output$model_random_effects_ui <- renderUI({
    div(
      checkboxInput(
        inputId = "include_tank_time_slope",
        label   = "Include tank time slope (1 + time_wk | tank)",
        value   = TRUE
      ),
      helpText("Adds a random intercept and a random slope for centered time within tank.", style = "margin-top:-6px;")
    )
  })

  # -----------------------------------------------------------------------------
  # Helper: robust normalization (never errors if a column is missing)
  # -----------------------------------------------------------------------------
  normalize_categories_for_lmer <- function(df) {
    # returns a column vector; if no candidate exists, returns rep(NA, nrow(df))
    safe_get <- function(d, ..., .default_type = "character") {
      cols <- c(...)
      present <- cols[cols %in% names(d)]
      if (length(present)) {
        out <- d[[present[1]]]
        return(out)
      }
      if (.default_type == "integer") {
        return(rep(NA_integer_, nrow(d)))
      }
      if (.default_type == "numeric") {
        return(rep(NA_real_, nrow(d)))
      }
      rep(NA_character_, nrow(d))
    }

    df %>%
      dplyr::mutate(
        fiber_type          = stringr::str_trim(stringr::str_to_lower(as.character(safe_get(., "fiber_type", "fiberType", "fibertype")))),
        chem_treatment      = stringr::str_trim(stringr::str_to_lower(as.character(safe_get(., "chem_treatment", "chemtreatment", "treatment")))),
        tissue_type         = stringr::str_trim(stringr::str_to_lower(as.character(safe_get(., "tissue_type", "tissue")))),
        sample_type         = stringr::str_trim(stringr::str_to_lower(as.character(safe_get(., "sample_type", "sampletype")))),
        fiber_concentration = as.character(safe_get(., "fiber_concentration", "fiberconcentration", "conc", "dose")),
        week                = suppressWarnings(as.integer(safe_get(., "week", .default_type = "integer"))),
        tank                = safe_get(., "tank", "tank_id"),
        experiment_batch    = safe_get(., "experiment_batch", "batch", "experiment_id")
      )
  }

  # --- Helpers for EMMeans ---

  # Before computing EMMeans, warn if the model is stale for current filters
  ensure_model_current <- function() {
    cur <- current_fit_state()
    last <- last_fit_state()
    is_current <- !is.null(last) && isTRUE(identical(cur, last))
    if (!is_current) {
      showNotification("Filters/random-effects changed since the last fit. Click 'Run Mixed Model' to update before EMMeans.",
        type = "warning", duration = 6
      )
    }
    is_current
  }

  # Standardize CI column names across emmeans versions/types
  standardize_emm_columns <- function(df) {
    if ("asymp.LCL" %in% names(df) && !"lower.CL" %in% names(df)) {
      df <- dplyr::rename(df, lower.CL = asymp.LCL)
    }
    if ("asymp.UCL" %in% names(df) && !"upper.CL" %in% names(df)) {
      df <- dplyr::rename(df, upper.CL = asymp.UCL)
    }
    df
  }

  # Build 'at=' list consistent with the model’s dose encoding
  build_emm_at_list <- function(selected_dose, use_dose_as_factor = FALSE) {
    selected_dose <- as.character(selected_dose %||% "0")
    if (isTRUE(use_dose_as_factor)) {
      list(dose_factor = selected_dose)
    } else {
      if (identical(selected_dose, "0")) {
        list(dose_log10 = 0, is_control = 1)
      } else {
        dose_num <- suppressWarnings(as.numeric(selected_dose))
        dose_val <- ifelse(isTRUE(is.finite(dose_num)) && dose_num > 0, log10(dose_num), 0)
        list(dose_log10 = dose_val, is_control = 0)
      }
    }
  }

  # Build a safe emmeans 'specs' based on what is actually in the fitted model
  make_safe_emm_specs <- function(request, model) {
    # What variables are in the model frame (i.e., truly in the formula)?
    mf_names <- names(stats::model.frame(model))
    # Interaction labels in the fixed-effects RHS
    term_labels <- attr(stats::terms(model), "term.labels")
    has_trt_week <- any(grepl("^chem_treatment:week$|^week:chem_treatment$", term_labels))
    has_trt_fiber <- any(grepl("^chem_treatment:fiber_type$|^fiber_type:chem_treatment$", term_labels))
    has_fiber_week <- any(grepl("^fiber_type:week$|^week:fiber_type$", term_labels))

    have_trt <- "chem_treatment" %in% mf_names
    have_fiber <- "fiber_type" %in% mf_names
    have_week <- "week" %in% mf_names

    # Helper to condition on another factor when an interaction is in the model
    simple_if_interacts <- function(lhs, rhs, has_interaction, have_rhs) {
      if (has_interaction && have_rhs) {
        # Simple effect: lhs within rhs
        return(reformulate(termlabels = paste(lhs, "|", rhs)))
      }
      reformulate(lhs)
    }

    # Map requested control to a safe spec present in the model
    out <- switch(request,
      treatment = if (have_trt) simple_if_interacts("chem_treatment", "week", has_trt_week, have_week) else reformulate("1"),
      fiber_type = if (have_fiber) simple_if_interacts("fiber_type", "week", has_fiber_week, have_week) else reformulate("1"),
      week = if (have_week) reformulate("week") else reformulate("1"),
      treatment_week = {
        vars <- c("chem_treatment", "week")
        vars <- vars[vars %in% mf_names]
        if (length(vars) == 2) reformulate(vars) else if (length(vars) == 1) reformulate(vars) else reformulate("1")
      },
      treatment_fiber = {
        vars <- c("chem_treatment", "fiber_type")
        vars <- vars[vars %in% mf_names]
        if (length(vars) == 2) reformulate(vars) else if (length(vars) == 1) reformulate(vars) else reformulate("1")
      },
      fiber_week = {
        vars <- c("fiber_type", "week")
        vars <- vars[vars %in% mf_names]
        if (length(vars) == 2) reformulate(vars) else if (length(vars) == 1) reformulate(vars) else reformulate("1")
      },
      all_interactions = {
        vars <- c("chem_treatment", "fiber_type", "week")
        vars <- vars[vars %in% mf_names]
        if (length(vars) >= 2) reformulate(vars) else if (length(vars) == 1) reformulate(vars) else reformulate("1")
      },
      # default
      {
        if (have_trt) reformulate("chem_treatment") else reformulate("1")
      }
    )

    out
  }

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
    df <- final_data # Pre-processed data

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
  }) %>% bindCache(
    input$assay, input$fiber, input$sample_type, input$week,
    input$fiber_concentration, input$treatment, input$plate_replicate, input$mode
  )

  assay_filtered <- assay_filtered_cached %>% bindEvent(input$run_analysis)

  # Physical data filter - Use reactive with proper caching
  phys_filtered_cached <- reactive({
    df <- physical_master # Pre-processed data

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
  }) %>% bindCache(
    input$fiber_type_phys, input$week_phys, input$endpoint,
    input$tissue_type_phys, input$fiber_concentration_phys, input$treatment_phys
  )

  # Now use bindEvent
  phys_filtered <- phys_filtered_cached %>% bindEvent(input$run_phys_analysis)

  # ============================================================================
  # DATA SUMMARIZATION - Using helper functions from helpers.R
  # ============================================================================

  summarized <- reactive({
    req(input$active_dataset == "assay")
    df <- assay_filtered()
    if (nrow(df) == 0) {
      return(NULL)
    }

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
    display_cols <- c(
      "fiber_type", "week", "endpoint", "tank", "sample",
      "treatment", "fiber_concentration", "fiber_group", "value"
    )

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

    if ("endpoint" %in% names(df) && "value" %in% names(df)) {
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
        facet_wrap(~assay_type) +
        theme_minimal() +
        labs(
          title = glue("Assay Activity by {str_to_title(gsub('_', ' ', x_axis_var))}"),
          x = str_to_title(gsub("_", " ", x_axis_var)), y = "Mean Activity", fill = "Fiber Concentration"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } else {
      x_axis_var <- input$x_axis_phys %||% "week"
      x_axis_factor <- paste0("factor(", x_axis_var, ")")

      if (input$mode_phys == "baseline") {
        ggplot(df, aes_string(x = x_axis_factor, y = "value", fill = "fiber_concentration")) +
          geom_boxplot(outlier.shape = NA, alpha = 0.7) +
          geom_jitter(width = 0.18, alpha = 0.6, size = 2) +
          facet_wrap(~endpoint, scales = "free_y") +
          theme_minimal() +
          labs(
            title = glue("Physical Baseline by {str_to_title(gsub('_', ' ', x_axis_var))}"),
            x = str_to_title(gsub("_", " ", x_axis_var)), y = "Value", fill = "Fiber Concentration"
          ) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      } else {
        ggplot(df, aes_string(x = x_axis_factor, y = "mean_value", fill = "fiber_concentration")) +
          geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.7) +
          geom_point(aes_string(group = "fiber_concentration"),
            position = position_dodge(width = 0.8),
            size = 2, color = "black", show.legend = FALSE
          ) +
          facet_wrap(~endpoint, scales = "free_y") +
          theme_minimal() +
          labs(
            title = glue("Physical Inter CV by {str_to_title(gsub('_', ' ', x_axis_var))}"),
            x = str_to_title(gsub("_", " ", x_axis_var)), y = "Mean Value", fill = "Fiber Concentration"
          ) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
    }
  })
  
  # --- Reactive data for the box plot ---
  viz_box_df <- reactive({
    df <- get_visualization_data()
    validate(need(!is.null(df) && is.data.frame(df), "No data available for visualization."))
    
    # Ensure we have a usable concentration_factor, but do not disturb tank
    df <- ensure_viz_concentration(df)
    
    # Make sure tank is a factor of the original ids (1..21), not 1..3 labels
    if (!("tank" %in% names(df))) {
      df$tank <- factor(paste0("tank_", seq_len(nrow(df))))
    } else {
      df$tank <- as.factor(df$tank)
    }
    
    # Pick Y
    y_var <- input$viz_y_var %||% "mean_activity"
    if (!(y_var %in% names(df))) {
      numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
      validate(need(length(numeric_cols) > 0, "No numeric outcome columns found."))
      y_var <- numeric_cols[[1]]
    }
    
    # Build the grouping columns to PRESERVE during aggregation
    group_cols <- c("concentration_factor", "tank")
    if (identical(input$viz_facet_by, "week") && "week" %in% names(df)) {
      group_cols <- c(group_cols, "week")
    } else if (identical(input$viz_facet_by, "fiber_type")) {
      if ("fiber_group" %in% names(df)) {
        group_cols <- c(group_cols, "fiber_group")
      } else if ("fiber_type" %in% names(df)) {
        group_cols <- c(group_cols, "fiber_type")
      }
    }
    
    # Optional: aggregate to one mean per (concentration, tank, [week/fiber_group])
    if (isTRUE(input$viz_aggregate_by_tank)) {
      df <- df %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
        dplyr::summarise(!!y_var := mean(.data[[y_var]], na.rm = TRUE), .groups = "drop")
    }
    
    # Rebuild the 1/3–3/3 shapes AFTER aggregation using the preserved tank ids
    df <- ensure_viz_concentration(df)
    
    # Keep only plotting columns
    keep_cols <- c("concentration_factor", "tank_replicate_of_3", y_var)
    if ("week" %in% group_cols) keep_cols <- c(keep_cols, "week")
    if ("fiber_group" %in% group_cols) keep_cols <- c(keep_cols, "fiber_group")
    if ("fiber_type" %in% group_cols) keep_cols <- c(keep_cols, "fiber_type")
    
    df %>% dplyr::select(dplyr::all_of(keep_cols))
  })
  
  # --- Renderer: box + hollow points, week/fiber_group on X with concentration grouped inside ---
  output$viz_boxpoints <- renderPlot({
    df <- viz_box_df()
    y_var <- input$viz_y_var %||% "mean_activity"
    
    validate(need(nrow(df) > 0, "No rows to plot."))
    validate(need(any(is.finite(df[[y_var]])), "All values are NA/Inf for the selected outcome."))
    
    facet_var <- input$viz_facet_by %||% "none"
    
    # Decide X-axis structure
    if (facet_var == "week" && "week" %in% names(df)) {
      # X = week factor, fill/color = concentration
      df$week <- as.factor(df$week)
      df$week <- forcats::fct_inorder(df$week)
      x_var <- "week"
      x_lab <- "Week"
      group_var <- "concentration_factor"
    } else if (facet_var == "fiber_type") {
      # Check for fiber_group first (composite column), then fall back to fiber_type
      if ("fiber_group" %in% names(df)) {
        df$fiber_group <- as.factor(df$fiber_group)
        df$fiber_group <- forcats::fct_inorder(df$fiber_group)
        x_var <- "fiber_group"
        x_lab <- "Fiber group"
      } else if ("fiber_type" %in% names(df)) {
        df$fiber_type <- as.factor(df$fiber_type)
        df$fiber_type <- forcats::fct_inorder(df$fiber_type)
        x_var <- "fiber_type"
        x_lab <- "Fiber type"
      } else {
        # Fallback if neither column exists
        x_var <- "concentration_factor"
        x_lab <- "Fiber concentration"
      }
      group_var <- "concentration_factor"
    } else {
      # No faceting: X = concentration
      x_var <- "concentration_factor"
      x_lab <- "Fiber concentration"
      group_var <- "concentration_factor"
    }
    
    # Brighter manual colors for concentration
    conc_levels <- levels(df$concentration_factor)
    n_conc <- length(conc_levels)
    bright_pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
    if (n_conc > length(bright_pal)) bright_pal <- rep(bright_pal, length.out = n_conc)
    conc_colors <- stats::setNames(bright_pal[seq_len(n_conc)], conc_levels)
    
    # Base plot
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]]))
    
    # Boxplot layer: fill by concentration, grouped by x_var
    p <- p + ggplot2::geom_boxplot(
      ggplot2::aes(fill = concentration_factor),
      width = 0.65, outlier.shape = NA, alpha = 0.7,
      position = ggplot2::position_dodge(width = 0.75)
    )
    
    # Hollow points layer: shape 21 (circle), 22 (square), 24 (triangle)
    # fill = NA (transparent), color = concentration (border)
    p <- p + ggplot2::geom_point(
      ggplot2::aes(color = concentration_factor, shape = tank_replicate_of_3, group = concentration_factor),
      size = 2.8, stroke = 1.4, fill = NA,
      position = ggplot2::position_jitterdodge(jitter.width = 0, dodge.width = 0.75)
    )
    
    # Hollow shape scale: 21, 22, 24 with NA fill
    p <- p + ggplot2::scale_shape_manual(
      name   = "Tank replicate",
      values = c("1" = 21, "2" = 22, "3" = 24),
      breaks = c("1","2","3"),
      labels = c("1/3","2/3","3/3"),
      drop   = FALSE
    )
    
    # Manual color/fill scales
    p <- p + ggplot2::scale_fill_manual(name = "Concentration", values = conc_colors, drop = FALSE)
    p <- p + ggplot2::scale_color_manual(name = "Concentration", values = conc_colors, drop = FALSE)
    
    # Labels
    p <- p + ggplot2::labs(
      title = paste0(y_var, " by ", x_lab),
      x     = x_lab,
      y     = y_var
    )
    
    p + theme_publication()
  })
  
  observe({
    df <- try(viz_box_df(), silent = TRUE)
    if (inherits(df, "data.frame")) {
      cat("viz: shape levels =", paste(levels(df$tank_replicate_of_3), collapse = ", "), "\n")
    }
  })
  
  output$download_viz_boxpoints_png <- downloadHandler(
    filename = function() {
      # Use UI selections and switch IDs by dataset to avoid "mixed" in Physical mode
      active_ds <- input$active_dataset %||% "assay"
      
      if (identical(active_ds, "assay")) {
        fiber_sel  <- tolower(input$fiber %||% character(0))          # "cotton"/"pet"
        trt_sel    <- tolower(input$treatment %||% character(0))       # "Treated"/"Untreated"/"Control"
        sample_sel <- tolower(input$sample_type %||% character(0))     # "gills"/"gland"/"hemolymph"
      } else {
        fiber_sel  <- tolower(input$fiber_type_phys %||% character(0)) # "cotton"/"pet"
        trt_sel    <- tolower(input$treatment_phys %||% character(0))  # "control"/"untreated"/"treated"
        # tissue selector only appears for mf_counts; if absent, fall back to endpoint tag
        sample_sel <- tolower(input$tissue_type_phys %||% character(0))
        if (length(sample_sel) == 0L) {
          sample_sel <- tolower(input$endpoint %||% character(0))      # e.g., "shell_length", "mf_counts"
        }
      }
      
      get_fiber_tag <- function(vals) {
        vals <- intersect(vals, c("cotton", "pet"))
        if (length(vals) == 1L) if (vals == "cotton") "Cotton" else "PET" else "mixed"
      }
      get_trt_tag <- function(vals) {
        core <- setdiff(vals, c("control", "ctrl"))
        if (length(core) == 1L && core %in% c("treated", "untreated")) core[[1L]] else "mixed"
      }
      get_sample_tag <- function(vals) {
        if (length(vals) == 1L) vals[[1L]] else "mixed"
      }
      
      fiber_tag  <- get_fiber_tag(fiber_sel)
      trt_tag    <- get_trt_tag(trt_sel)
      sample_tag <- get_sample_tag(sample_sel)
      
      paste0("boxplot_", fiber_tag, "_", trt_tag, "_", sample_tag, ".png")
    },
    content = function(file) {
      p <- build_viz_boxpoints_plot()
      ggplot2::ggsave(
        filename = file,
        plot     = p,
        width    = input$viz_export_width %||% 10,
        height   = input$viz_export_height %||% 6,
        dpi      = input$viz_export_dpi %||% 600,
        units    = "in",
        bg       = "white",
        device   = "png"
      )
    }
  )

  # Dot plot - assay only
  output$dot_plot <- renderPlot({
    req(input$active_dataset == "assay")
    df <- assay_filtered()
    req(!is.null(df), nrow(df) > 0)

    required_cols <- c("fiber_group", "calculated_concentration", "fiber_concentration", "assay_type")
    if (!all(required_cols %in% colnames(df))) {
      return(NULL)
    }

    ggplot(df, aes(x = fiber_group, y = calculated_concentration, color = fiber_concentration)) +
      geom_jitter(width = 0.3, alpha = 0.7, size = 2) +
      facet_wrap(~assay_type) +
      theme_minimal() +
      labs(
        title = "Sample-Level Distribution by Group and Fiber Concentration",
        x = "Fiber Group", y = "Calculated Concentration", color = "Fiber Concentration"
      ) +
      theme(plot.title = element_text(size = 14, face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))
  })

  # Info text
  output$info <- renderText({
    df <- active_table()
    if (is.null(df)) {
      return("No data for current filters / mode.")
    }

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
    filename = function() {
      paste0("ECOTOX_", input$mode, "_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      df <- summarized()
      writexl::write_xlsx(if (is.null(df)) data.frame(Message = "No data") else df, file)
    }
  )

  output$download_phys_data <- downloadHandler(
    filename = function() {
      paste0("ECOTOX_physical_", Sys.Date(), ".xlsx")
    },
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
    wks <- input$weeks_include %||% c("1", "3", "5")
    df <- df %>%
      dplyr::filter(as.character(week) %in% wks) %>%
      normalize_controls_and_dose() %>%
      droplevels()

    df
  }) %>%
    bindCache(
      input$regression_dataset,
      input$regression_endpoint,
      input$reg_sample_types,
      input$reg_tissues,
      input$weeks_include
    )

  # Keep existing alias if used downstream
  regression_data <- regression_data_cached

  # Regression model fitting
  regression_model <- eventReactive(input$run_regression, {
    tryCatch(
      {
        df <- regression_data()
        req(nrow(df) > 0)

        has_week <- dplyr::n_distinct(df$week) >= 2

        # Base fixed-effects formula
        form <- outcome ~ fiber_type * chem_treatment +
          is_control * fiber_type +
          dose_log10

        # Add week main effect and interactions only when estimable
        if (has_week) {
          form <- update(form, . ~ . + week + fiber_type:week + chem_treatment:week)
          if (isTRUE(input$include_three_way)) {
            form <- update(form, . ~ . + fiber_type:chem_treatment:week)
          }
        } else {
          form <- update(form, . ~ . + fiber_type:chem_treatment)
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
      },
      error = function(e) {
        list(error = e$message)
      }
    )
  })

  # -----------------------------------------------------------------------------
  # Mixed Effects endpoint choices (assay uses assay_type, physical uses endpoint)
  # -----------------------------------------------------------------------------
  observeEvent(input$lmer_dataset,
    {
      ds_raw <- input$lmer_dataset %||% "assay"
      ds_tag <- tolower(as.character(ds_raw))
      is_assay <- grepl("assay", ds_tag)
      is_physical <- grepl("phys", ds_tag)

      src <- NULL
      ep_col <- NULL

      if (is_assay && exists("final_data") && is.data.frame(final_data) && nrow(final_data) > 0) {
        src <- final_data
        ep_col <- "assay_type"
      } else if (is_physical && exists("physical_master") && is.data.frame(physical_master) && nrow(physical_master) > 0) {
        src <- physical_master
        ep_col <- "endpoint"
      } else if (exists("final_data") && "assay_type" %in% names(final_data)) {
        # Fallback: prefer assay if present
        src <- final_data
        ep_col <- "assay_type"
      } else if (exists("physical_master") && "endpoint" %in% names(physical_master)) {
        src <- physical_master
        ep_col <- "endpoint"
      }

      # Update both possible UI ids to be safe
      if (is.null(src) || is.null(ep_col)) {
        updateSelectInput(session, "lmerendpoint", choices = character(0), selected = NULL)
        updateSelectInput(session, "lmer_endpoint", choices = character(0), selected = NULL)
        return()
      }

      ep_vals <- sort(unique(as.character(src[[ep_col]])))
      sel <- if (length(ep_vals)) ep_vals[1] else NULL
      updateSelectInput(session, "lmerendpoint", choices = ep_vals, selected = sel)
      updateSelectInput(session, "lmer_endpoint", choices = ep_vals, selected = sel)
    },
    ignoreInit = FALSE
  )

  # # Also trigger once after UI is ready to clear "Loading..."
  # session$onFlushed(function() {
  #   isolate({
  #     # Reuse the observer above by poking the reactive dependency
  #
  #   })
  # }, once = TRUE)

  # 2) Tissues UI for physical mf_counts (unchanged behavior)
  output$lmer_tissues_ui <- renderUI({
    req(input$lmer_dataset == "physical")
    ep <- input$lmer_endpoint %||% ""
    if (identical(ep, "mf_counts")) {
      checkboxGroupInput(
        inputId = "lmer_tissue_types",
        label = "Tissue types:",
        choices = tissue_types_physical, # from global.R
        selected = tissue_types_physical,
        inline = TRUE
      )
    } else {
      checkboxGroupInput(
        inputId = "lmer_tissue_types",
        label = "Tissue types:",
        choices = c("All"),
        selected = "All",
        inline = TRUE
      ) %>%
        shinyjs::disable()
    }
  })

  # -----------------------------------------------------------------------------
  # Mixed Effects: model-ready data using UI shim (works for assay and physical)
  # -----------------------------------------------------------------------------
  lmer_data_cached <- reactive({
    validate(
      need(!is.null(input$lmer_dataset), "Select a dataset for Mixed Effects."),
      need(!is.null(input$lmer_endpoint) && nzchar(input$lmer_endpoint), "Select an endpoint for Mixed Effects.")
    )
    
    ds <- tolower(as.character(input$lmer_dataset))
    if (identical(ds, "assay")) {
      validate(need(exists("final_data") && is.data.frame(final_data) && nrow(final_data) > 0, "No assay data loaded."))
      df <- final_data %>%
        dplyr::rename(endpoint_ui = assay_type) %>%
        dplyr::filter(.data$endpoint_ui == input$lmer_endpoint)
      outcome_col <- "calculated_concentration"
    } else if (identical(ds, "physical")) {
      validate(need(exists("physical_master") && is.data.frame(physical_master) && nrow(physical_master) > 0, "No physical data loaded."))
      df <- physical_master %>%
        dplyr::rename(endpoint_ui = endpoint) %>%
        dplyr::filter(.data$endpoint_ui == input$lmer_endpoint)
      outcome_col <- "value"
    } else {
      validate(need(FALSE, sprintf("Unknown dataset tag: %s", ds)))
    }
    
    # Normalize key fields (tank, experiment_batch, week, factors, etc.)
    df <- normalize_categories_for_lmer(df)
    
    # Apply Mixed Effects UI filters
    to_lower <- function(x) if (is.null(x)) character(0) else tolower(as.character(x))
    if (!is.null(input$lmer_sample_types) && length(input$lmer_sample_types)) {
      df <- dplyr::filter(df, sample_type %in% to_lower(input$lmer_sample_types))
    }
    if (!is.null(input$lmer_fiber_types) && length(input$lmer_fiber_types)) {
      df <- dplyr::filter(df, fiber_type %in% to_lower(input$lmer_fiber_types))
    }
    if (!is.null(input$lmer_treatments) && length(input$lmer_treatments)) {
      df <- dplyr::filter(df, chem_treatment %in% to_lower(input$lmer_treatments))
    }
    if (!is.null(input$lmer_weeks_include) && length(input$lmer_weeks_include)) {
      df <- dplyr::filter(df, as.character(week) %in% as.character(input$lmer_weeks_include))
    }
    if (!is.null(input$lmer_concentrations) && length(input$lmer_concentrations)) {
      df <- dplyr::filter(df, as.character(fiber_concentration) %in% as.character(input$lmer_concentrations))
    }
    if (identical(ds, "physical") && identical(input$lmer_endpoint, "mf_counts") &&
        !is.null(input$lmer_tissue_types) && length(input$lmer_tissue_types)) {
      df <- dplyr::filter(df, tissue_type %in% to_lower(input$lmer_tissue_types))
    }
    
    # ---- SAFE Harmonization of experiment_batch and tank ----
    safe_col <- function(d, col_name) {
      if (col_name %in% names(d)) as.character(d[[col_name]]) else rep(NA_character_, nrow(d))
    }
    cand_batch <- list(
      safe_col(df, "experiment_batch"),
      safe_col(df, "batch"),
      safe_col(df, "batch_id"),
      safe_col(df, "experiment"),
      safe_col(df, "experiment_id"),
      safe_col(df, "exp_batch"),
      safe_col(df, "exp")
    )
    df$experiment_batch <- dplyr::coalesce(!!!cand_batch)
    df$experiment_batch <- ifelse(is.na(df$experiment_batch), NA_character_, trimws(tolower(df$experiment_batch)))
    df$tank <- if ("tank" %in% names(df)) trimws(tolower(as.character(df$tank))) else NA_character_
    
    # Availability log and gentle notice
    n_tank  <- if ("tank" %in% names(df)) dplyr::n_distinct(stats::na.omit(df$tank)) else 0L
    n_batch <- if ("experiment_batch" %in% names(df)) dplyr::n_distinct(stats::na.omit(df$experiment_batch)) else 0L
    message(sprintf("[lmer-data] groups seen post-normalize: tank=%d, batch=%d", n_tank, n_batch))
    if (n_batch < 2L) {
      showNotification(
        "experiment_batch missing or has < 2 distinct levels after filters; model will omit batch.",
        type = "message", duration = 6
      )
    }
    
    # Ensure 'time_wk' (numeric) and 'time_wk_z' (centered per week) exist for random slope
    df$time_wk <- suppressWarnings(as.numeric(df$week))
    df$time_wk_z <- as.numeric(scale(suppressWarnings(as.numeric(df$week)), 
                                     center = TRUE, scale = FALSE))
    
    # Create dose encodings and numeric outcome BEFORE finite filter
    validate(need(outcome_col %in% names(df), sprintf("Column '%s' not found in selected data.", outcome_col)))
    df <- df %>%
      dplyr::mutate(
        dose_factor = factor(fiber_concentration, levels = c("0", "100", "1000", "10000"), ordered = TRUE),
        is_control  = dplyr::if_else(fiber_concentration %in% c("0","control","ctrl"), 1L, 0L, missing = 0L),
        dose_log10  = dplyr::case_when(
          fiber_concentration %in% c("0","control","ctrl") ~ 0,
          TRUE ~ log10(suppressWarnings(as.numeric(fiber_concentration)))
        ),
        outcome = suppressWarnings(as.numeric(.data[[outcome_col]]))
      ) %>%
      dplyr::filter(is.finite(outcome)) %>%
      droplevels()
    
    validate(need(nrow(df) > 0, "No rows remain after Mixed Effects filters."))
    df
  }) %>% bindCache(
    input$lmer_dataset, input$lmer_endpoint, input$lmer_fiber_types,
    input$lmer_sample_types, input$lmer_tissue_types, input$lmer_treatments,
    input$lmer_concentrations, input$lmer_weeks_include
  )
  
  lmer_data <- lmer_data_cached

  # --- Mixed Effects model event: safe formula construction ---------------------
  lmer_model <- eventReactive(input$run_lmer, {
    df <- lmer_data_cached()  # already normalized and filtered
    
    # Optional grouping inspection
    info <- inspect_grouping(df)
    message(sprintf(
      "[lmer] groups: tank=%d, batch=%d, fiber=%d, confounded(batch,fiber)=%s",
      info$n_tank, info$n_batch, info$n_fiber, info$confounded_batch_fiber
    ))
    print(info$cross_tab)
    
    # Basic availability log
    n_tank  <- if ("tank" %in% names(df)) dplyr::n_distinct(df$tank) else 0L
    n_batch <- if ("experiment_batch" %in% names(df)) dplyr::n_distinct(df$experiment_batch) else 0L
    message(sprintf("[lmer] groups: tank=%d, batch=%d", n_tank, n_batch))
    
    validate(need(
      "tank" %in% names(df),
      "Random-effects require a 'tank' column in the selected data."
    ))
    
    # ---------- Dose encoding ----------
    use_dose_as_factor <- isTRUE(input$lmer_dose_as_factor)
    dose_term <- if (use_dose_as_factor) "dose_factor" else "dose_log10"
    
    # ---------- Guarantee time_wk_z availability (if user bypassed caching somehow) ----------
    if (!("time_wk_z" %in% names(df))) {
      time_raw <- if ("time_wk" %in% names(df)) df$time_wk else df$week
      if (!is.numeric(time_raw)) time_raw <- suppressWarnings(as.numeric(as.character(time_raw)))
      df$time_wk_z <- as.numeric(scale(time_raw, center = TRUE, scale = FALSE))
    }
    
    # ---------- Estimability gates for fixed effects ----------
    include_fiber <- "fiber_type" %in% names(df) && dplyr::n_distinct(df$fiber_type) >= 2L
    include_trt   <- "chem_treatment" %in% names(df) && dplyr::n_distinct(df$chem_treatment) >= 2L
    include_week  <- ("week" %in% names(df)) && dplyr::n_distinct(df$week) >= 2L
    if (include_week) df$week <- factor(df$week)
    
    rhs_terms <- c(
      if (include_fiber) "fiber_type",
      if (include_trt)   "chem_treatment",
      dose_term,
      if (!use_dose_as_factor) "is_control",
      if (include_week)  "week"
    )
    if (include_fiber && include_trt)  rhs_terms <- c(rhs_terms, "fiber_type:chem_treatment")
    if (include_week && include_fiber) rhs_terms <- c(rhs_terms, "fiber_type:week")
    if (include_week && include_trt)   rhs_terms <- c(rhs_terms, "chem_treatment:week")
    if (isTRUE(input$lmer_include_three_way) && include_fiber && include_trt && include_week) {
      rhs_terms <- c(rhs_terms, "fiber_type:chem_treatment:week")
    }
    if (!use_dose_as_factor && isTRUE(input$lmer_dose_by_fiber) && include_fiber) {
      rhs_terms <- c(rhs_terms, "dose_log10:fiber_type")
    }
    if (!use_dose_as_factor && isTRUE(input$lmer_dose_by_treat) && include_trt) {
      rhs_terms <- c(rhs_terms, "dose_log10:chem_treatment")
    }
    if (use_dose_as_factor && isTRUE(input$lmer_dose_by_fiber) && include_fiber) {
      rhs_terms <- c(rhs_terms, "dose_factor:fiber_type")
    }
    if (use_dose_as_factor && isTRUE(input$lmer_dose_by_treat) && include_trt) {
      rhs_terms <- c(rhs_terms, "dose_factor:chem_treatment")
    }
    if (length(rhs_terms) == 0L) rhs_terms <- "1"
    
    include_time_slope <- isTRUE(input$include_tank_time_slope)
    
    # ---------- Guarantee time_wk_z exists before building random terms ----------
    if (!("time_wk_z" %in% names(df))) {
      time_raw <- if ("time_wk" %in% names(df)) df$time_wk else df$week
      if (!is.numeric(time_raw)) time_raw <- suppressWarnings(as.numeric(as.character(time_raw)))
      df$time_wk_z <- as.numeric(scale(time_raw, center = TRUE, scale = FALSE))
      message("[lmer] created time_wk_z on the fly (mean=", round(mean(df$time_wk_z, na.rm=TRUE), 2), ")")
    }
    
    # Build random terms with your existing helper
    random_terms <- build_random_terms(df, include_time_slope = include_time_slope)
    
    fixed_str  <- paste("outcome ~", paste(rhs_terms, collapse = " + "))
    full_form  <- stats::as.formula(paste(fixed_str, "+", paste(random_terms, collapse = " + ")))
    
    fit_obj <- lme4::lmer(full_form, data = df,
                          control = lme4::lmerControl(optimizer = "bobyqa",
                                                      optCtrl = list(maxfun = 20000)))
    
    # Optional: fallback if singular with slope ON
    slope_was_used <- include_time_slope && grepl("time_wk_z", random_terms, fixed = TRUE)
    if (slope_was_used && isTRUE(lme4::isSingular(fit_obj, tol = 1e-4))) {
      showNotification("Singular fit with tank time slope; refitting without slope for stability.", type = "warning", duration = 6)
      random_terms2 <- build_random_terms(df = df, include_time_slope = FALSE)
      full_form2 <- stats::as.formula(paste(fixed_str, "+", paste(random_terms2, collapse = " + ")))
      fit_obj <- lme4::lmer(
        formula = full_form2, data = df, REML = TRUE,
        control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 20000))
      )
      message("[lmer] refit without slope; random: ", paste(random_terms2, collapse = " + "))
    }
    
    fit_obj
  })

  # Optional: quick logger to surface row counts when running lmer
  observeEvent(input$run_lmer, {
    df <- try(lmer_data_cached(), silent = TRUE)
    if (!inherits(df, "try-error") && is.data.frame(df)) {
      message(sprintf(
        "[lmer] n=%d endpoint=%s fibers={%s} trts={%s} conc={%s} weeks={%s}",
        nrow(df),
        input$lmer_endpoint %||% "NULL",
        paste(tolower(input$lmer_fiber_types %||% character(0)), collapse = ","),
        paste(tolower(input$lmer_treatments %||% character(0)), collapse = ","),
        paste(as.character(input$lmer_concentrations %||% character(0)), collapse = ","),
        paste(as.character(input$lmer_weeks_include %||% character(0)), collapse = ",")
      ))
    }
  })
  
  # ---- Mixed effects diagnostics helpers (singularity-safe) ----
  theme_diagnostics <- function() {
    ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold")
      )
  }
  
  # Build tidy diagnostics for an lm
  build_lm_diagnostics <- function(model) {
    aug <- broom::augment(model)
    aug$.sqrt_std_resid <- sqrt(abs(aug$.std.resid))
    # 1) Residuals vs Fitted
    p1 <- ggplot2::ggplot(aug, ggplot2::aes(x = .fitted, y = .resid)) +
      ggplot2::geom_hline(yintercept = 0, color = "#999999") +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "#2c7fb8") +
      ggplot2::labs(title = "Residuals vs Fitted", x = "Fitted values", y = "Residuals") +
      theme_diagnostics()
    # 2) Normal Q-Q
    p2 <- ggplot2::ggplot(aug, ggplot2::aes(sample = .std.resid)) +
      ggplot2::stat_qq(alpha = 0.8) +
      ggplot2::stat_qq_line(color = "#2c7fb8") +
      ggplot2::labs(title = "Normal Q-Q (standardized residuals)", x = "Theoretical", y = "Standardized residuals") +
      theme_diagnostics()
    # 3) Scale-Location
    p3 <- ggplot2::ggplot(aug, ggplot2::aes(x = .fitted, y = .sqrt_std_resid)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "#2c7fb8") +
      ggplot2::labs(title = "Scale-Location", x = "Fitted values", y = "√|standardized residuals|") +
      theme_diagnostics()
    # 4) Residuals vs Leverage with Cook’s distance
    p4 <- ggplot2::ggplot(aug, ggplot2::aes(x = .hat, y = .std.resid, size = .cooksd)) +
      ggplot2::geom_hline(yintercept = 0, color = "#999999") +
      ggplot2::geom_point(alpha = 0.8, color = "#2c3e50") +
      ggplot2::scale_size_continuous(name = "Cook's D") +
      ggplot2::labs(title = "Residuals vs Leverage", x = "Leverage", y = "Standardized residuals") +
      theme_diagnostics()
    
    # Fit stats and quick flags
    bp <- tryCatch(lmtest::bptest(model), error = function(e) NULL)
    sh <- tryCatch(shapiro.test(stats::residuals(model)), error = function(e) NULL)
    stats_df <- data.frame(
      metric  = c("AIC","BIC","logLik","sigma","adj_R2","nobs","BP_p","Shapiro_p"),
      value   = c(AIC(model), BIC(model), as.numeric(logLik(model)),
                  sigma(model), summary(model)$adj.r.squared,
                  stats::nobs(model),
                  if (!is.null(bp)) unname(bp$p.value) else NA_real_,
                  if (!is.null(sh)) unname(sh$p.value) else NA_real_)
    )
    
    list(plots = list(p1 = p1, p2 = p2, p3 = p3, p4 = p4), stats = stats_df)
  }
  
  # Suggest minimal adjustments from diagnostics
  suggest_lm_adjustments <- function(model, data) {
    aug <- broom::augment(model)
    bp_p <- tryCatch(lmtest::bptest(model)$p.value, error = function(e) NA_real_)
    sh_p <- tryCatch(shapiro.test(stats::residuals(model))$p.value, error = function(e) NA_real_)
    
    # Simple curvature flag via LOESS RMSE vs linear
    loess_rmse <- sqrt(mean(residuals(loess(.resid ~ .fitted, data = aug))^2, na.rm = TRUE))
    lin_rmse   <- sqrt(mean(aug$.resid^2, na.rm = TRUE))
    curvature  <- (loess_rmse + 1e-8) < (lin_rmse * 0.97)  # 3% improvement threshold
    
    list(
      heteroskedastic = is.finite(bp_p) && bp_p < 0.05,
      nonnormal       = is.finite(sh_p) && sh_p < 0.05,
      curvature       = isTRUE(curvature)
    )
  }
  
  build_lmer_diagnostics <- function(model, sims = 250L, tol = 1e-4, r2_tol = 1e-9) {
    # Singularity check up front
    singular_flag <- tryCatch(lme4::isSingular(model, tol = tol), error = function(e) NA)
    
    # 1) DHARMa residuals
    sim_res <- tryCatch(
      DHARMa::simulateResiduals(model, n = sims, seed = 123),
      error = function(e) {
        warning("DHARMa::simulateResiduals failed: ", conditionMessage(e))
        NULL
      }
    )
    
    # 2) Residuals vs fitted (Pearson), silence geom_smooth message by setting formula
    df_rf <- data.frame(
      fitted = as.numeric(stats::fitted(model)),
      resid  = as.numeric(stats::residuals(model, type = "pearson"))
    )
    p_res_fit <- ggplot2::ggplot(df_rf, ggplot2::aes(x = fitted, y = resid)) +
      ggplot2::geom_hline(yintercept = 0, color = "#999999") +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "#2c7fb8") +
      ggplot2::labs(title = "Pearson residuals vs Fitted", x = "Fitted", y = "Pearson residuals") +
      theme_diagnostics()
    
    # 3) Random effects (BLUPs) caterpillar with safe fallback
    re_df <- tryCatch(
      broom.mixed::tidy(model, effects = "ran_vals"),
      error = function(e) {
        warning("broom.mixed::tidy ran_vals failed: ", conditionMessage(e), "; falling back to lme4::ranef")
        re_list <- lme4::ranef(model, condVar = FALSE)
        out <- lapply(names(re_list), function(g) {
          mat <- as.data.frame(re_list[[g]])
          term_names <- colnames(mat)
          keep_col <- if ("(Intercept)" %in% term_names) "(Intercept)" else term_names[1L]
          data.frame(
            group    = g,
            level    = rownames(re_list[[g]]),
            term     = keep_col,
            estimate = as.numeric(mat[[keep_col]]),
            stringsAsFactors = FALSE
          )
        })
        do.call(rbind, out)
      }
    )
    p_re <- if (is.null(re_df) || !nrow(re_df)) {
      ggplot2::ggplot() + ggplot2::theme_void() +
        ggplot2::labs(title = "Random effects not available for this model")
    } else {
      ggplot2::ggplot(re_df, ggplot2::aes(x = estimate, y = forcats::fct_rev(level))) +
        ggplot2::geom_point() +
        ggplot2::facet_wrap(~ group, scales = "free_y") +
        ggplot2::labs(title = "Random effects (BLUPs) by grouping factor", x = "Estimate", y = "Level") +
        theme_diagnostics()
    }
    
    # 4) Fit statistics (exported generics + robust fallbacks)
    remle_val <- suppressWarnings(tryCatch(
      as.numeric(stats::deviance(model)),
      error = function(e) {
        s <- tryCatch(summary(model), error = function(e2) NULL)
        if (!is.null(s) && !is.null(s$REMLcrit)) as.numeric(s$REMLcrit) else NA_real_
      }
    ))
    
    # Prefer explicit Nakagawa to avoid API drift; tolerate singularity
    r2_obj <- tryCatch(performance::r2_nakagawa(model, tolerance = r2_tol), error = function(e) NULL)
    r2_marg <- tryCatch(unname(r2_obj$R2_marginal), error = function(e) NA_real_)
    # If singular, conditional may be undefined; keep NA to avoid misleading numbers
    r2_cond <- if (isTRUE(singular_flag)) NA_real_ else tryCatch(unname(r2_obj$R2_conditional), error = function(e) NA_real_)
    
    # Optional: variance components snapshot (helps interpret singular fits)
    vc_df <- tryCatch(as.data.frame(lme4::VarCorr(model)), error = function(e) NULL)
    
    stats_df <- data.frame(
      metric = c("AIC","BIC","logLik","REMLcrit","nobs","R2_marginal","R2_conditional","singular"),
      value  = c(AIC(model), BIC(model), as.numeric(logLik(model)),
                 remle_val, stats::nobs(model),
                 r2_marg, r2_cond, isTRUE(singular_flag))
    )
    
    note <- if (isTRUE(singular_flag)) {
      "Singular fit: random-effect variance ~ 0; conditional R2 is omitted and BLUPs may be near zero."
    } else {
      ""
    }
    
    list(sim = sim_res, p_res_fit = p_res_fit, p_re = p_re, stats = stats_df, note = note, varcomps = vc_df)
  }
  
  # ---- Mixed effects diagnostics (renderers; UI IDs unchanged) ----
  lmer_diag <- reactive({
    fit <- lmer_model()
    validate(need(inherits(fit, "lmerMod"), "Run Mixed Model first."))
    build_lmer_diagnostics(fit, sims = 250L)
  })
  
  output$lmer_dharma <- renderPlot({
    d <- lmer_diag()
    validate(need(!is.null(d$sim), "DHARMa residuals unavailable for this fit."))
    # Use graphics generic so the S3 method plot.DHARMa is dispatched
    graphics::plot(d$sim)
  })
  
  output$lmer_resid_fitted   <- renderPlot({ lmer_diag()$p_res_fit })
  output$lmer_re_caterpillar <- renderPlot({ lmer_diag()$p_re })
  
  output$lmer_fit_stats <- DT::renderDataTable({
    d  <- lmer_diag()
    df <- d$stats
    cap <- if (nzchar(d$note)) htmltools::tags$div(style = "color:#9a6700;", d$note) else NULL
    DT::datatable(df, options = list(dom = "t", pageLength = 10), rownames = FALSE, caption = cap) %>%
      DT::formatRound("value", 6)
  })
  
  # ============================================================================
  # ENHANCED EMMEANS WITH COLORED SIGNIFICANCE & DID ANALYSIS
  # ============================================================================

  # EMMEANS calculation for regular regression
  emmeans_results <- eventReactive(input$run_emmeans, {
    tryCatch(
      {
        mdl <- regression_model()
        validate(need(!is.null(mdl) && !inherits(mdl, "list"), "Run the regression model first"))

        # Dose settings for 'at'
        dose_val <- switch(input$emmeans_dose,
          "0" = 0,
          "2" = log10(100),
          "3" = log10(1000),
          "4" = log10(10000)
        )
        is_ctrl <- ifelse(dose_val == 0, 1, 0)
        at_list <- list(dose_log10 = dose_val, is_control = is_ctrl)

        # Is 'week' in the fitted model?
        model_terms <- attr(terms(mdl), "term.labels")
        has_week <- any(grepl("\\bweek\\b", model_terms))

        # Map UI to a valid emmeans spec, avoiding 'week' if absent
        emmeans_spec <- switch(input$emmeans_by,
          "treatment"        = ~ chem_treatment | fiber_type,
          "fiber_type"       = ~fiber_type,
          "week"             = if (has_week) ~week else ~ chem_treatment | fiber_type,
          "treatment_week"   = if (has_week) ~ chem_treatment | week + fiber_type else ~ chem_treatment | fiber_type,
          "treatment_fiber"  = ~ chem_treatment | fiber_type,
          "fiber_week"       = if (has_week) ~ fiber_type | week else ~fiber_type,
          "all_interactions" = if (has_week) ~ chem_treatment | fiber_type + week else ~ chem_treatment | fiber_type
        )

        emm <- emmeans::emmeans(mdl, specs = emmeans_spec, at = at_list)
        list(emmeans = emm, emmeans_df = as.data.frame(emm))
      },
      error = function(e) {
        message("Emmeans error: ", e$message)
        list(error = e$message)
      }
    )
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
          is.na(p.value) ~ "Unable to compute",
          p.value < 0.001 ~ "Yes (***)",
          p.value < 0.01 ~ "Yes (**) ",
          p.value < 0.05 ~ "Yes (*)",
          TRUE ~ "No"
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
          c("#28a745", "#28a745", "#28a745", "#ffffff", "#ffe6cc")
        ),
        color = DT::styleEqual(
          c("Yes (***)", "Yes (**) ", "Yes (*)", "No", "Unable to compute"),
          c("#ffffff", "#ffffff", "#ffffff", "#000000", "#000000")
        )
      ) %>%
      {
        if (all(c("estimate", "SE", "t.ratio") %in% names(dfp))) {
          DT::formatRound(., columns = c("estimate", "SE", "t.ratio"), digits = 4)
        } else {
          .
        }
      } %>%
      {
        if ("p.value" %in% names(dfp)) {
          DT::formatRound(., columns = "p.value", digits = 6)
        } else {
          .
        }
      }
  })

  # ENHANCED EMMEANS TABLE WITH SIGNIFICANCE COLORS
  output$emmeans_table <- DT::renderDataTable({
    results <- emmeans_results()
    validate(need(!is.null(results$emmeans_df), "Click 'Calculate EM Means' to generate results"))

    df_emm <- results$emmeans_df %>%
      mutate(across(where(is.numeric), ~ signif(.x, 4)))

    DT::datatable(
      df_emm,
      options = list(pageLength = 15, scrollX = TRUE, dom = "Blfrtip"),
      filter = "top",
      rownames = FALSE
    ) %>%
      DT::formatRound(columns = c("emmean", "SE", "lower.CL", "upper.CL"), digits = 4)
  })

  # ENHANCED EMMEANS PLOT
  output$emmeans_plot <- renderPlot({
    results <- emmeans_results()
    validate(need(!is.null(results$emmeans), "Calculate EM Means first"))

    tryCatch(
      {
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
      },
      error = function(e) {
        # Fallback if emmeans plot fails - use FIXED geom_errorbar (not geom_errorbarh)
        df <- results$emmeans_df
        ggplot(df, aes(x = emmean, y = rownames(df))) +
          geom_point(size = 3, color = "#2c3e50") +
          # FIXED: Use geom_errorbar with orientation instead of deprecated geom_errorbarh
          geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL),
            width = 0.2, orientation = "y"
          ) +
          theme_minimal() +
          labs(
            title = "Estimated Marginal Means",
            x = "Estimated Mean",
            y = "Groups"
          ) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
    )
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
      return(list(~chem_treatment, fiber_type = input$trend_fiber_select))
    } else if (filter_type == "treatment") {
      return(list(~fiber_type, chem_treatment = input$trend_treatment_select))
    } else if (filter_type == "both") {
      return(list(
        fiber_type = input$trend_fiber_select,
        chem_treatment = input$trend_treatment_select
      ))
    }
  })

  # Helper: map trend_concentration to at-list for emmeans/emtrends
  trend_at_list <- reactive({
    conc <- input$trend_concentration
    if (is.null(conc) || conc == "all") {
      return(NULL)
    }

    if (conc == "0") {
      # control: dose_log10 = 0, is_control = 1
      list(dose_log10 = 0, is_control = 1)
    } else {
      # non-control dose: log10(conc), is_control = 0
      conc_num <- suppressWarnings(as.numeric(conc))
      if (isTRUE(is.na(conc_num)) || conc_num <= 0) {
        return(NULL)
      }
      list(dose_log10 = log10(conc_num), is_control = 0)
    }
  })

  # ===========================
  # LINEAR TRENDS TEXT OUTPUT
  # ===========================
  output$regression_trends <- renderPrint({
    req(input$show_trends == TRUE)
    req(regression_model())

    tryCatch(
      {
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
          cat(
            "Analyzing trends for: ", toupper(input$trend_fiber_select),
            " fiber (both treatments)\n\n"
          )
        } else if (input$trend_filter_type == "treatment") {
          cat(
            "Analyzing trends for: ", toupper(input$trend_treatment_select),
            " (both fiber types)\n\n"
          )
        } else {
          cat(
            "Analyzing trends for: ", toupper(input$trend_fiber_select),
            " fiber, ", toupper(input$trend_treatment_select), " treatment\n\n"
          )
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
            var = "week", at = at_list
          )
        } else if (input$trend_filter_type == "fiber") {
          trends <- emmeans::emtrends(mdl_numeric, ~chem_treatment,
            var = "week",
            at = c(list(fiber_type = input$trend_fiber_select), at_list)
          )
        } else if (input$trend_filter_type == "treatment") {
          trends <- emmeans::emtrends(mdl_numeric, ~fiber_type,
            var = "week",
            at = c(list(chem_treatment = input$trend_treatment_select), at_list)
          )
        } else {
          trends <- emmeans::emtrends(mdl_numeric, ~1,
            var = "week",
            at = c(
              list(
                fiber_type = input$trend_fiber_select,
                chem_treatment = input$trend_treatment_select
              ),
              at_list
            )
          )
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
      },
      error = function(e) {
        cat("Error computing trends:\n")
        cat(conditionMessage(e), "\n\n")
        cat("This may occur if:\n")
        cat("- Model cannot be refit with numeric week\n")
        cat("- Insufficient data for selected filter\n")
        cat("- Model lacks necessary interaction terms\n")
      }
    )
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
    if (!any(grepl("week", model_terms, ignore.case = TRUE))) {
      return(NULL)
    }

    at_list <- trend_at_list()

    if (input$trend_filter_type == "none") {
      emm_week <- emmeans::emmeans(mdl, ~ week + fiber_type + chem_treatment, at = at_list)
    } else if (input$trend_filter_type == "fiber") {
      emm_week <- emmeans::emmeans(mdl, ~ week + chem_treatment,
        at = c(list(fiber_type = input$trend_fiber_select), at_list)
      )
    } else if (input$trend_filter_type == "treatment") {
      emm_week <- emmeans::emmeans(mdl, ~ week + fiber_type,
        at = c(list(chem_treatment = input$trend_treatment_select), at_list)
      )
    } else {
      emm_week <- emmeans::emmeans(mdl, ~week,
        at = c(list(
          fiber_type = input$trend_fiber_select,
          chem_treatment = input$trend_treatment_select
        ), at_list)
      )
    }

    emm_df <- as.data.frame(emm_week)
    emm_df$week_num <- if (is.factor(emm_df$week)) as.numeric(as.character(emm_df$week)) else emm_df$week

    dose_subtitle <- if (!is.null(input$trend_concentration) && input$trend_concentration != "all") {
      paste0("Evaluated at ", input$trend_concentration, " mf/L")
    } else {
      "Averaged across doses"
    }

    library(ggplot2)
    if (input$trend_filter_type == "none") {
      p <- ggplot(emm_df, aes(
        x = week_num, y = emmean, color = chem_treatment,
        linetype = chem_treatment, shape = chem_treatment
      )) +
        geom_line(size = 1.2) +
        geom_point(size = 3.5) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, size = 0.8, alpha = 0.7) +
        facet_wrap(~fiber_type, nrow = 1, scales = "free_y") +
        scale_color_manual(values = get_treatment_palette()) +
        theme_publication() +
        labs(
          title = "Estimated Marginal Means Trends Over Time",
          subtitle = paste("Linear trends; error bars = 95% CI •", dose_subtitle),
          x = "Week", y = "Estimated Mean Outcome",
          color = "Treatment", linetype = "Treatment", shape = "Treatment"
        )
    } else if (input$trend_filter_type == "fiber") {
      p <- ggplot(emm_df, aes(
        x = week_num, y = emmean, color = chem_treatment,
        linetype = chem_treatment, shape = chem_treatment
      )) +
        geom_line(size = 1.3) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, size = 0.9, alpha = 0.7) +
        scale_color_manual(values = get_treatment_palette()) +
        theme_publication() +
        labs(
          title = paste0("Trends Over Time (", toupper(input$trend_fiber_select), " Fiber)"),
          subtitle = paste("Linear trends; error bars = 95% CI •", dose_subtitle),
          x = "Week", y = "Estimated Mean Outcome",
          color = "Treatment", linetype = "Treatment", shape = "Treatment"
        )
    } else if (input$trend_filter_type == "treatment") {
      p <- ggplot(emm_df, aes(
        x = week_num, y = emmean, color = fiber_type,
        linetype = fiber_type, shape = fiber_type
      )) +
        geom_line(size = 1.3) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, size = 0.9, alpha = 0.7) +
        scale_color_manual(values = get_fiber_palette()) +
        theme_publication() +
        labs(
          title = paste0("Trends Over Time (", toupper(input$trend_treatment_select), " Treatment)"),
          subtitle = paste("Linear trends; error bars = 95% CI •", dose_subtitle),
          x = "Week", y = "Estimated Mean Outcome",
          color = "Fiber Type", linetype = "Fiber Type", shape = "Fiber Type"
        )
    } else {
      p <- ggplot(emm_df, aes(x = week_num, y = emmean)) +
        geom_line(size = 1.5, color = "#0073C2") +
        geom_point(size = 5, color = "#0073C2") +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, size = 1, color = "#0073C2", alpha = 0.7) +
        theme_publication() +
        labs(
          title = paste0("Trend (", toupper(input$trend_fiber_select), ", ", toupper(input$trend_treatment_select), ")"),
          subtitle = paste("Linear trend; error bars = 95% CI •", dose_subtitle),
          x = "Week", y = "Estimated Mean Outcome"
        )
    }

    p + scale_x_continuous(breaks = unique(emm_df$week_num))
  })

  output$regression_trend_plot <- renderPlot({
    p <- regression_trend_plot_reactive()
    if (is.null(p)) {
      plot.new()
      text(0.5, 0.5, "Week not included in model (check week checkboxes in 'Include Weeks')",
        cex = 1.3, col = "#d9534f", font = 2
      )
    } else {
      p
    }
  })

  output$download_regression_trend_png <- create_plot_download_handler(regression_trend_plot_reactive, "regression_trends", "png")
  output$download_regression_trend_pdf <- create_plot_download_handler(regression_trend_plot_reactive, "regression_trends", "pdf")
  output$download_regression_trend_tiff <- create_plot_download_handler(regression_trend_plot_reactive, "regression_trends", "tiff")
  output$download_regression_trend_svg <- create_plot_download_handler(regression_trend_plot_reactive, "regression_trends", "svg")

  # =============================================================================
  # Combined Treatment Analysis (server)
  # - Uses snake_case IDs
  # - Harmonizes legacy/snake_case columns
  # - Avoids circular gating by rendering even when hidden and gating by button
  # =============================================================================

  # Small helper to format a quick data status line in the UI
  output$combined_data_info <- renderText({
    if (is.null(input$combined_endpoint) || identical(input$combined_endpoint, "loading")) {
      return("No data loaded") # initial state
    }
    # Try to count rows safely
    n <- tryCatch(
      {
        df <- combined_base_data()
        if (is.null(df)) 0L else nrow(df)
      },
      error = function(e) 0L
    )
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
        create_enhanced_treatment_categories() %>% # normalize fiber/treatment fields
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
      wks <- c("1", "3", "5")
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
    bindCache(
      input$active_dataset,
      input$combined_endpoint,
      input$combined_sample_types,
      input$combined_tissues,
      input$combined_weeks
    )

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
    req(nrow(df_combined) > 0) # Also check we have data left after filtering

    df_combined
  })

  # =============================================================================
  #         OBSERVER TO DYNAMICALLY UPDATE THE REFERENCE LEVEL CHOICES
  # =============================================================================

  observe({
    # Trigger when combined_base_data changes (endpoint/week selection)
    req(input$combined_endpoint, length(input$combined_weeks) > 0)

    tryCatch(
      {
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
            selected = control_levels[1] # Default to first control
          )

          cat("✓ “ Reference level options updated:", paste(control_levels, collapse = ", "), "\n")
        }
      },
      error = function(e) {
        # Silently fail if data not ready
        NULL
      }
    )
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
      need(
        n_lvls >= 2,
        sprintf("Not enough treatment levels after filtering (have %d). Try different weeks or endpoint.", n_lvls)
      )
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
        sd_value = stats::sd(outcome, na.rm = TRUE),
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

    tryCatch(
      {
        # Calculate emmeans for treatment_combined factor
        emm <- emmeans::emmeans(mdl, specs = ~treatment_combined)

        # Convert to data frame
        emm_df <- as.data.frame(emm) %>%
          dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 4)))

        list(
          emmeans = emm,
          emmeans_df = emm_df
        )
      },
      error = function(e) {
        message("Emmeans calculation error: ", e$message)
        list(error = e$message)
      }
    )
  })

  # ============================================================================
  # COMBINED MODEL: PAIRWISE (compute once) + REACTIVE FILTERS
  # ============================================================================

  # 1) Compute the complete Tukey-adjusted pairwise set once per model run
  combined_pairwise_all <- eventReactive(input$run_combined_model, {
    emm_res <- combined_emmeans_data()
    validate(need(!is.null(emm_res$emmeans), "Calculate emmeans first"))

    tryCatch(
      {
        emm <- emm_res$emmeans

        # All pairwise, no filtering here
        pw <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey")
        df <- as.data.frame(pw)

        # Helper parsing for filtering
        parts <- stringr::str_split_fixed(df$contrast, " - ", 2)
        a <- parts[, 1]
        b <- parts[, 2]

        parse_fiber <- function(x) {
          dplyr::case_when(
            stringr::str_detect(x, regex("cotton", ignore_case = TRUE)) ~ "Cotton",
            stringr::str_detect(x, regex("\\bpet\\b|pet\\s|\\spet\\b", ignore_case = TRUE)) ~ "PET",
            TRUE ~ NA_character_
          )
        }
        parse_trt <- function(x) {
          dplyr::case_when(
            stringr::str_detect(x, regex("^Control| Control", ignore_case = TRUE)) ~ "Control",
            stringr::str_detect(x, regex("Untreated", ignore_case = TRUE)) ~ "Untreated",
            stringr::str_detect(x, regex("Treated", ignore_case = TRUE)) ~ "Treated",
            TRUE ~ NA_character_
          )
        }
        parse_conc <- function(x) stringr::str_extract(x, "\\b\\d{1,5}\\b") # 0,100,1000,10000, etc.

        df <- df %>%
          dplyr::mutate(
            significant = dplyr::case_when(
              p.value < 0.001 ~ "***",
              p.value < 0.01 ~ "**",
              p.value < 0.05 ~ "*",
              TRUE ~ ""
            ),
            a = a, b = b,
            fiber_a = parse_fiber(a), fiber_b = parse_fiber(b),
            treat_a = parse_trt(a), treat_b = parse_trt(b),
            conc_a = parse_conc(a), conc_b = parse_conc(b)
          )

        df
      },
      error = function(e) {
        message("Pairwise-all error: ", e$message)
        validate(need(FALSE, paste("Pairwise computation failed:", e$message)))
        NULL
      }
    )
  })

  # 2) Apply lightweight, instant filters based on current UI selections
  combined_pairwise_filtered <- reactive({
    # p-adjust method for subset families; default Tukey
    p_adjust <- input$combined_p_adjust %||% "tukey" # "tukey", "holm", "BH", etc.

    ftype <- input$combined_comparison_filter %||% "all"
    ref_level <- input$combined_reference

    mdl <- combined_model()
    validate(need(!is.null(mdl), "Model not available"))

    # Helper parsers to keep downstream filters working
    parse_fiber <- function(x) {
      dplyr::case_when(
        stringr::str_detect(x, regex("cotton", ignore_case = TRUE)) ~ "Cotton",
        stringr::str_detect(x, regex("\\bpet\\b|pet\\s|\\spet\\b", ignore_case = TRUE)) ~ "PET",
        TRUE ~ NA_character_
      )
    }
    parse_trt <- function(x) {
      dplyr::case_when(
        stringr::str_detect(x, regex("^Control| Control", ignore_case = TRUE)) ~ "Control",
        stringr::str_detect(x, regex("Untreated", ignore_case = TRUE)) ~ "Untreated",
        stringr::str_detect(x, regex("Treated", ignore_case = TRUE)) ~ "Treated",
        TRUE ~ NA_character_
      )
    }
    parse_conc <- function(x) stringr::str_extract(x, "\\b\\d{1,5}\\b")

    # Subset-aware recomputation so the adjustment matches the displayed family
    if (ftype %in% c("fiber", "fiber_control")) {
      emm_all <- emmeans::emmeans(mdl, specs = ~treatment_combined)
      all_lvls <- levels(emm_all)$treatment_combined
      fiber_pick <- input$combined_filter_fiber %||% "PET"

      keep_lvls <- all_lvls[stringr::str_detect(all_lvls, fiber_pick, negate = FALSE)]
      emm_fiber <- emmeans::emmeans(
        mdl,
        specs = ~treatment_combined,
        exclude = setdiff(all_lvls, keep_lvls)
      )

      pw <- if (ftype == "fiber") {
        emmeans::contrast(emm_fiber, method = "pairwise", adjust = p_adjust)
      } else {
        fiber_control <- ifelse(
          tolower(fiber_pick) %in% c("pet", "polyester", "pet (polyester)"),
          "Control Pet", "Control Cotton"
        )
        emmeans::contrast(emm_fiber, method = "trt.vs.ctrl", ref = fiber_control, adjust = p_adjust)
      }

      out <- as.data.frame(pw)

      # Parse contrast into a/b and attach minimal tags
      ab <- stringr::str_split_fixed(out$contrast, " - ", 2)
      out$a <- ab[, 1]
      out$b <- ab[, 2]
      out$fiber_a <- parse_fiber(out$a)
      out$fiber_b <- parse_fiber(out$b)
      out$treat_a <- parse_trt(out$a)
      out$treat_b <- parse_trt(out$b)
      out$conc_a <- parse_conc(out$a)
      out$conc_b <- parse_conc(out$b)

      # Always provide 'significant' for DT styling
      out <- out %>%
        dplyr::mutate(
          significant = dplyr::case_when(
            is.na(p.value) ~ "",
            p.value < 0.001 ~ "***",
            p.value < 0.01 ~ "**",
            p.value < 0.05 ~ "*",
            TRUE ~ ""
          ),
          fiber_subset = fiber_pick
        ) %>%
        dplyr::mutate(dplyr::across(
          tidyselect::any_of(c("estimate", "SE", "t.ratio", "p.value", "lower.CL", "upper.CL")),
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
            "p.value" %in% names(out) & !is.na(p.value) & p.value < 0.01 ~ "**",
            "p.value" %in% names(out) & !is.na(p.value) & p.value < 0.05 ~ "*",
            TRUE ~ ""
          )
        )
    }

    out %>%
      dplyr::mutate(dplyr::across(
        tidyselect::any_of(c("estimate", "SE", "t.ratio", "p.value", "lower.CL", "upper.CL")),
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
        columnDefs = list(list(className = "dt-center", targets = "_all"))
      ),
      rownames = FALSE,
      caption = "Estimated Marginal Means for Combined Treatment Levels"
    ) %>%
      DT::formatRound(columns = c("emmean", "SE", "lower.CL", "upper.CL"), digits = 4)
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
    if ("Cotton" %in% present) pal["Cotton"] <- "#2E7D32"
    if ("PET" %in% present) pal["PET"] <- "#1565C0"
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
  output$download_combined_emm_png <- create_plot_download_handler(combined_emm_plot_reactive, "combined_emmeans", "png")
  output$download_combined_emm_pdf <- create_plot_download_handler(combined_emm_plot_reactive, "combined_emmeans", "pdf")
  output$download_combined_emm_tiff <- create_plot_download_handler(combined_emm_plot_reactive, "combined_emmeans", "tiff")
  output$download_combined_emm_svg <- create_plot_download_handler(combined_emm_plot_reactive, "combined_emmeans", "svg")

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
        columnDefs = list(list(className = "dt-center", targets = "_all"))
      ),
      rownames = FALSE,
      caption = paste0("Pairwise Comparisons (α = ", alpha, ") - Filter: ", ftype)
    ) %>%
      DT::formatRound(columns = c("estimate", "SE", "t.ratio", "p.value"), digits = 4) %>%
      DT::formatStyle(
        "p.value",
        backgroundColor = DT::styleInterval(
          c(0.001, 0.01, 0.05),
          c("#ffebee", "#fff3e0", "#fffde7", "white")
        )
      ) %>%
      DT::formatStyle(
        "significant",
        fontWeight = "bold",
        color = DT::styleEqual(c("***", "**", "*", ""), c("red", "orange", "blue", "black"))
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
          p.value < 0.01 ~ "p < 0.01",
          p.value < 0.05 ~ "p < 0.05",
          TRUE ~ "Not significant"
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
          "p < 0.01" = "#F57C00",
          "p < 0.05" = "#1976D2",
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
  output$download_combined_pairwise_png <- create_plot_download_handler(combined_pairwise_plot_reactive, "combined_pairwise", "png", height = 10)
  output$download_combined_pairwise_pdf <- create_plot_download_handler(combined_pairwise_plot_reactive, "combined_pairwise", "pdf", height = 10)
  output$download_combined_pairwise_tiff <- create_plot_download_handler(combined_pairwise_plot_reactive, "combined_pairwise", "tiff", height = 10)
  output$download_combined_pairwise_svg <- create_plot_download_handler(combined_pairwise_plot_reactive, "combined_pairwise", "svg", height = 10)

  # ============================================================================
  # OUTPUT: COMBINED DIAGNOSTICS
  # ============================================================================

  # ---- Diagnostics for Combined Treatment (lm) ----
  build_lm_diagnostics <- function(model) {
    aug <- broom::augment(model)
    aug$.sqrt_std_resid <- sqrt(abs(aug$.std.resid))
    
    p1 <- ggplot2::ggplot(aug, ggplot2::aes(x = .fitted, y = .resid)) +
      ggplot2::geom_hline(yintercept = 0, color = "#999999") +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#2c7fb8") +
      ggplot2::labs(title = "Residuals vs Fitted", x = "Fitted", y = "Residuals") +
      theme_diagnostics()
    
    p2 <- ggplot2::ggplot(aug, ggplot2::aes(sample = .std.resid)) +
      ggplot2::stat_qq(alpha = 0.8) +
      ggplot2::stat_qq_line(color = "#2c7fb8") +
      ggplot2::labs(title = "Normal Q-Q (standardized residuals)", x = "Theoretical", y = "Standardized residuals") +
      theme_diagnostics()
    
    p3 <- ggplot2::ggplot(aug, ggplot2::aes(x = .fitted, y = .sqrt_std_resid)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#2c7fb8") +
      ggplot2::labs(title = "Scale-Location", x = "Fitted", y = "√|standardized residuals|") +
      theme_diagnostics()
    
    p4 <- ggplot2::ggplot(aug, ggplot2::aes(x = .hat, y = .std.resid, size = .cooksd)) +
      ggplot2::geom_hline(yintercept = 0, color = "#999999") +
      ggplot2::geom_point(alpha = 0.8, color = "#2c3e50") +
      ggplot2::scale_size_continuous(name = "Cook's D") +
      ggplot2::labs(title = "Residuals vs Leverage", x = "Leverage", y = "Standardized residuals") +
      theme_diagnostics()
    
    bp <- tryCatch(lmtest::bptest(model), error = function(e) NULL)
    
    stats_df <- data.frame(
      metric  = c("AIC","BIC","logLik","sigma","adj_R2","nobs","BP_p"),
      value   = c(AIC(model), BIC(model), as.numeric(logLik(model)),
                  sigma(model), summary(model)$adj.r.squared,
                  stats::nobs(model),
                  if (!is.null(bp)) unname(bp$p.value) else NA_real_)
    )
    list(plots = list(p1 = p1, p2 = p2, p3 = p3, p4 = p4), stats = stats_df)
  }
  
  combined_diag <- reactive({
    m <- combined_model()
    validate(need(!is.null(m), "Run Combined Treatment Model first."))
    build_lm_diagnostics(m)
  })
  
  output$combined_resid_fitted   <- renderPlot({ combined_diag()$plots$p1 })
  output$combined_qq             <- renderPlot({ combined_diag()$plots$p2 })
  output$combined_scale_location <- renderPlot({ combined_diag()$plots$p3 })
  output$combined_leverage       <- renderPlot({ combined_diag()$plots$p4 })
  
  output$combined_fit_stats <- DT::renderDataTable({
    df <- combined_diag()$stats
    DT::datatable(df, options = list(dom = "t", pageLength = 10), rownames = FALSE) %>%
      DT::formatRound("value", 5)
  })
  
  # # Diagnostics plot with safe graphics state and explicit namespaces
  # output$combined_diagnostics <- renderPlot({
  #   mdl <- combined_model()
  #   validate(need(!is.null(mdl), "Run model first to see diagnostics"))
  # 
  #   # Save and restore par to avoid leaking graphics settings
  #   op <- par(no.readonly = TRUE)
  #   on.exit(par(op), add = TRUE)
  # 
  #   par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  # 
  #   # 1) Residuals vs Fitted
  #   plot(
  #     fitted(mdl), residuals(mdl),
  #     xlab = "Fitted Values", ylab = "Residuals",
  #     main = "Residuals vs Fitted", pch = 16,
  #     col = scales::alpha("steelblue", 0.6) # explicit namespace
  #   )
  #   abline(h = 0, col = "red", lty = 2, lwd = 2)
  # 
  #   # 2) Normal Q-Q
  #   qqnorm(
  #     residuals(mdl),
  #     main = "Normal Q-Q Plot",
  #     pch = 16, col = scales::alpha("steelblue", 0.6)
  #   )
  #   qqline(residuals(mdl), col = "red", lwd = 2, lty = 2)
  # 
  #   # 3) Scale-Location
  #   sqrt_std_resid <- sqrt(abs(scale(residuals(mdl))))
  #   plot(
  #     fitted(mdl), sqrt_std_resid,
  #     xlab = "Fitted Values", ylab = expression(sqrt("|Standardized Residuals|")),
  #     main = "Scale-Location", pch = 16,
  #     col = scales::alpha("steelblue", 0.6)
  #   )
  #   lines(lowess(fitted(mdl), sqrt_std_resid), col = "red", lwd = 2)
  # 
  #   # 4) Cook's Distance
  #   cooks_d <- cooks.distance(mdl)
  #   n <- length(cooks_d)
  #   plot(
  #     seq_along(cooks_d), cooks_d,
  #     type = "h",
  #     xlab = "Observation Index", ylab = "Cook's Distance",
  #     main = "Cook's Distance (Influence)",
  #     col = ifelse(cooks_d > 4 / n, "red", "steelblue"), lwd = 2
  #   )
  #   abline(h = 4 / n, col = "red", lty = 2, lwd = 1.5)
  #   text(
  #     x = n * 0.7, y = max(cooks_d, na.rm = TRUE) * 0.9,
  #     labels = paste("Threshold:", round(4 / n, 4)),
  #     col = "red", cex = 0.8
  #   )
  # })

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
    req(input$active_dataset) # CORRECT: active_dataset with underscore!

    cat("\n=== Updating combined_endpoint choices ===\n")
    cat("Active dataset:", input$active_dataset, "\n")

    choices <- character(0)

    if (identical(input$active_dataset, "assay")) {
      # Try to get assay endpoints from final_data
      choices <- tryCatch(
        {
          req(final_data)
          sort(unique(final_data$assay_type))
        },
        error = function(e) {
          cat("Error accessing final_data:", e$message, "\n")
          c("ACP", "AChE", "ALP", "Bradford", "CAT", "SOD")
        }
      )
      cat("Assay endpoints:", paste(choices, collapse = ", "), "\n")
    } else if (identical(input$active_dataset, "physical")) {
      # Try to get physical endpoints from physical_master
      choices <- tryCatch(
        {
          req(physical_master)
          sort(unique(physical_master$endpoint))
        },
        error = function(e) {
          cat("Error accessing physical_master:", e$message, "\n")
          c("bci", "clearance_rate", "mf_counts", "respiration_rate")
        }
      )
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

    data_count <- tryCatch(
      {
        req(combined_base_data())
        nrow(combined_base_data())
      },
      error = function(e) 0
    )

    if (data_count > 0) {
      paste0("✓ “ ", data_count, " observations loaded")
    } else {
      "⚠  No data for selected endpoint"
    }
  })

  # --- Mixed Effects: EMMeans calculation (safe spec + robust) ---
  lmer_emmeans_results <- eventReactive(input$lmer_run_emmeans, {
    tryCatch(
      {
        # Require a fitted mixed model first
        mdl <- lmer_model()
        validate(need(
          !is.null(mdl) && !inherits(mdl, "list"),
          "Run the mixed effects model first"
        ))

        # Match the model's dose encoding and build 'at='
        use_dose_as_factor <- isTRUE(input$lmer_dose_as_factor)
        dose_choice <- as.character(input$lmer_emmeans_dose %||% "0")
        at_list <- build_emm_at_list(dose_choice, use_dose_as_factor)

        # Build a spec that only uses variables present in the model
        requested <- as.character(input$lmer_emmeans_by %||% "treatment")
        emm_specs <- make_safe_emm_specs(requested, mdl)

        # If the spec reduced to ~1 (nothing estimable), fall back to treatment if available
        if (identical(format(emm_specs), "~1")) {
          mf_names <- names(stats::model.frame(mdl))
          if ("chem_treatment" %in% mf_names) emm_specs <- ~chem_treatment
        }

        # Try full emmeans with 'at' when appropriate
        emm <- tryCatch(
          {
            # Only pass 'at' if the model actually has the dose variables in use
            if (use_dose_as_factor && "dose_factor" %in% names(stats::model.frame(mdl))) {
              emmeans::emmeans(mdl, specs = emm_specs, at = at_list)
            } else if (!use_dose_as_factor && all(c("dose_log10", "is_control") %in% names(stats::model.frame(mdl)))) {
              emmeans::emmeans(mdl, specs = emm_specs, at = at_list)
            } else {
              emmeans::emmeans(mdl, specs = emm_specs)
            }
          },
          error = function(e1) {
            message("Primary emmeans failed, retrying without 'at=': ", e1$message)
            tryCatch(
              {
                emmeans::emmeans(mdl, specs = emm_specs)
              },
              error = function(e2) {
                message("Fallback emmeans failed, using simplest spec: ", e2$message)
                # Final fallback: single available factor, or stop
                mf_names <- names(stats::model.frame(mdl))
                if ("chem_treatment" %in% mf_names) {
                  emmeans::emmeans(mdl, specs = ~chem_treatment)
                } else if ("week" %in% mf_names) {
                  emmeans::emmeans(mdl, specs = ~week)
                } else {
                  stop("No suitable factor is available in the model for EMMeans.")
                }
              }
            )
          }
        )

        # Summarize via S3 method (do not namespace), then standardize CI names
        emm_tbl <- as.data.frame(summary(emm, infer = TRUE))
        emm_tbl <- standardize_emm_columns(emm_tbl)

        # Numeric cleanup for safety
        num_cols <- intersect(c("emmean", "SE", "df", "lower.CL", "upper.CL"), names(emm_tbl))
        emm_tbl[num_cols] <- lapply(emm_tbl[num_cols], function(x) suppressWarnings(as.numeric(x)))

        # Plot data: keep finite rows
        emm_plot <- dplyr::filter(
          emm_tbl,
          is.finite(emmean),
          (!"lower.CL" %in% names(emm_tbl)) | is.finite(lower.CL),
          (!"upper.CL" %in% names(emm_tbl)) | is.finite(upper.CL)
        )

        # Log the spec actually used
        req_txt <- paste(deparse(emm_specs), collapse = " ")
        message("[emmeans] specs used: ", req_txt)

        list(
          emmeans = emm,
          emmeans_df = emm_tbl,
          emmeans_plot = emm_plot,
          specification_used = req_txt
        )
      },
      error = function(e) {
        # IMPORTANT: error handler must be inside this same tryCatch(...) call
        message("EMMeans failed: ", e$message)
        list(error = e$message)
      }
    )
  })

  # ---- Mixed effects emmeans table (DT) ----
  output$lmer_emmeans_table <- DT::renderDataTable({
    results <- lmer_emmeans_results()
    validate(need(!is.null(results$emmeans_df), "Click 'Calculate EM Means' for mixed effects"))

    df_emm <- results$emmeans_df %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 4)))

    lower_col <- if ("lower.CL" %in% names(df_emm)) "lower.CL" else NULL
    upper_col <- if ("upper.CL" %in% names(df_emm)) "upper.CL" else NULL
    round_cols <- c("emmean", "SE", lower_col, upper_col)
    round_cols <- round_cols[round_cols %in% names(df_emm)]

    DT::datatable(
      df_emm,
      options = list(pageLength = 10, scrollX = TRUE),
      filter = "top",
      rownames = FALSE
    ) %>%
      (function(dtbl) {
        if (length(round_cols)) DT::formatRound(dtbl, columns = round_cols, digits = 4) else dtbl
      })
  })

  # ---- Mixed effects pairwise table (DT) ----
  output$lmer_pairwise_table <- DT::renderDataTable({
    m <- lmer_model()
    validate(need(!is.null(m), "Run the mixed effects model first."))

    adj_method <- input$lmer_adjustment_method %||% "tukey"

    if (identical(input$lmer_comparison_type, "did")) {
      # Simple DiD at selected dose
      dfp <- tryCatch(
        {
          use_dose_as_factor <- isTRUE(input$lmer_dose_as_factor)
          dose_choice <- as.character(input$lmer_emmeans_dose %||% "0")
          at_list <- build_emm_at_list(dose_choice, use_dose_as_factor)

          emm_tw <- emmeans::emmeans(
            m,
            specs = ~ chem_treatment | fiber_type + week, at = at_list
          )
          as.data.frame(emmeans::contrast(emm_tw, method = "pairwise", adjust = adj_method))
        },
        error = function(e) {
          data.frame(Message = paste0("DiD Error (Mixed Effects): ", e$message), stringsAsFactors = FALSE)
        }
      )
    } else {
      # Regular pairwise or trt vs ctrl using current EMMeans object
      emm_res <- lmer_emmeans_results()
      validate(need(!is.null(emm_res$emmeans), "Calculate EM Means first"))
      emm <- emm_res$emmeans

      if (identical(input$lmer_comparison_type, "control")) {
        dfp <- as.data.frame(emmeans::contrast(emm,
          method = "trt.vs.ctrl",
          ref = "untreated", adjust = adj_method
        ))
      } else {
        dfp <- as.data.frame(emmeans::contrast(emm, method = "pairwise", adjust = adj_method))
      }
    }

    # Standardize CI names if present, then format
    dfp <- standardize_emm_columns(dfp)

    if (!"p.value" %in% names(dfp)) dfp$p.value <- NA_real_

    dfp <- dfp %>%
      dplyr::mutate(
        dplyr::across(where(is.numeric), ~ signif(.x, 5)),
        Significant = dplyr::case_when(
          is.na(p.value) ~ "Unable to compute",
          p.value < 0.001 ~ "Yes (***)",
          p.value < 0.01 ~ "Yes (**)",
          p.value < 0.05 ~ "Yes (*)",
          TRUE ~ "No"
        ),
        Adjustment = adj_method
      )

    dt <- DT::datatable(
      dfp,
      options = list(pageLength = 15, scrollX = TRUE),
      filter = "top",
      rownames = FALSE
    ) %>%
      DT::formatStyle(
        "Significant",
        backgroundColor = DT::styleEqual(
          c("Yes (***)", "Yes (**)", "Yes (*)", "No", "Unable to compute"),
          c("#28a745", "#28a745", "#28a745", "#ffffff", "#ffe6cc")
        ),
        color = DT::styleEqual(
          c("Yes (***)", "Yes (**)", "Yes (*)", "No", "Unable to compute"),
          c("#ffffff", "#ffffff", "#ffffff", "#000000", "#000000")
        )
      ) %>%
      (function(dtbl) {
        for (nm in c("p.value", "estimate", "SE", "t.ratio", "z.ratio")) {
          if (nm %in% names(dfp)) dtbl <- DT::formatRound(dtbl, columns = nm, digits = if (nm == "p.value") 6 else 4)
        }
        dtbl
      })

    dt
  })

  # ============================================================================
  # 3. ENHANCED MIXED EFFECTS PLOT WITH BETTER ERROR HANDLING
  # ============================================================================

  # ---- Mixed effects emmeans plot (manual, robust) ----
  output$lmer_emmeans_plot <- renderPlot({
    # Require EMMeans results and no upstream error
    results <- lmer_emmeans_results()
    validate(need(!is.null(results$emmeans_df), "Calculate EM Means for mixed effects first"))
    validate(need(is.null(results$error), paste("Emmeans error:", results$error)))

    # Use pre-filtered, standardized data prepared by lmer_emmeans_results()
    df <- results$emmeans_plot
    validate(need(nrow(df) > 0, "No finite EMMeans rows to plot."))

    # Build readable y-axis labels from available grouping columns
    group_cols <- intersect(c("chem_treatment", "fiber_type", "week", "dose_factor"), names(df))
    if (length(group_cols) == 0L) group_cols <- "chem_treatment"
    df$group_label <- do.call(paste, c(df[group_cols], sep = " | "))

    # Draw point estimates and horizontal CIs (if available)
    p <- ggplot2::ggplot(df, ggplot2::aes(y = group_label, x = emmean)) +
      ggplot2::geom_point(size = 2.5, color = "#2c3e50") +
      ggplot2::labs(
        title = "Mixed Effects Estimated Marginal Means with 95% CI",
        x = "Estimated marginal mean",
        y = paste(group_cols, collapse = " | ")
      ) +
      ggplot2::theme_bw()

    # Add CIs only when both bounds are present
    if (all(c("lower.CL", "upper.CL") %in% names(df))) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(xmin = lower.CL, xmax = upper.CL),
        width = 0.15, orientation = "y", color = "#3498db"
      )
    }

    p
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
  
  # ---- Base regression diagnostics (lm) ----
  regression_diag <- reactive({
    # Pull the fitted lm from your existing eventReactive
    m <- regression_model()
    validate(need(!is.null(m) && inherits(m, "lm"),
                  "Run Regression Model first."))
    # Reuse the shared diagnostics helper already defined above
    build_lm_diagnostics(m)
  })
  
  # Plots for the base regression model
  output$regression_resid_fitted   <- renderPlot({ regression_diag()$plots$p1 })
  output$regression_qq             <- renderPlot({ regression_diag()$plots$p2 })
  output$regression_scale_location <- renderPlot({ regression_diag()$plots$p3 })
  output$regression_leverage       <- renderPlot({ regression_diag()$plots$p4 })
  
  # Fit stats table for the base regression model
  output$regression_fit_stats <- DT::renderDataTable({
    df <- regression_diag()$stats
    DT::datatable(df, options = list(dom = "t", pageLength = 10), rownames = FALSE) %>%
      DT::formatRound("value", 6)
  })
  
  # Diagnostics for the refitted model
  regression_diag_adjusted <- reactive({
    fit <- regression_model_adjusted()
    validate(need(!is.null(fit) && inherits(fit$model, "lm"), "Refit the Adjusted Model first."))
    build_lm_diagnostics(fit$model)
  })
  
  # Plots
  output$regression_adjusted_resid_fitted   <- renderPlot({ regression_diag_adjusted()$plots$p1 })
  output$regression_adjusted_qq             <- renderPlot({ regression_diag_adjusted()$plots$p2 })
  output$regression_adjusted_scale_location <- renderPlot({ regression_diag_adjusted()$plots$p3 })
  output$regression_adjusted_leverage       <- renderPlot({ regression_diag_adjusted()$plots$p4 })
  
  # Fit comparison (unchanged) plus a small banner if WLS weights were applied
  output$regression_adjusted_fit_stats <- DT::renderDataTable({
    fit <- regression_model_adjusted()
    cmp <- fit$comparison
    DT::datatable(cmp, options = list(dom = "t", pageLength = 8), rownames = FALSE) %>%
      DT::formatRound(c("base","adjusted","delta"), 6)
  })
  
  # Coefficient comparison (base vs adjusted)
  output$regression_adjusted_coef_table <- DT::renderDataTable({
    fit <- regression_model_adjusted()
    tbl <- fit$coef_cmp
    DT::datatable(tbl, options = list(dom = "t", pageLength = 20), rownames = FALSE) %>%
      DT::formatRound(c("base_estimate","base_se","base_p","adj_estimate","adj_se","adj_p"), 6)
  })
  
  # Human-readable summary of the refit
  output$regression_model_summary_adjusted <- renderPrint({
    fit <- regression_model_adjusted()
    validate(need(!is.null(fit), "Select options and click 'Refit Adjusted Model'."))
    cat("Adjusted model formula:\n", format(fit$formula_final), "\n\n")
    print(summary(fit$model))
    if (isTRUE(fit$used_weights)) {
      cat("\nNote: Weighted least squares was applied using dose_factor-level residual variances.\n")
    }
  })
  
  # --- Optional: Refit adjusted model when user clicks ---
  regression_model_adjusted <- eventReactive(input$run_regression_adjusted, {
    base_m <- regression_model()
    validate(need(!is.null(base_m) && inherits(base_m, "lm"), "Run Regression Model first."))
    
    # Use fuller data so both dose encodings are available if present
    df <- get_regression_refit_data(base_m)
    
    # Start from the base model's formula (single line)
    fml <- stats::formula(base_m) |> collapse_formula()
    
    # ---- 1) Outcome transform (log) ----
    if (isTRUE(input$reg_log_outcome)) {
      y_name <- all.vars(update(fml, . ~ 0))[1]
      if (all(df[[y_name]] > 0, na.rm = TRUE)) {
        df[[y_name]] <- log(df[[y_name]])
      } else {
        showNotification("Log transform skipped: outcome has non-positive values.", type = "warning", duration = 5)
      }
    }
    
    # Make sure dose columns exist before any formula surgery
    df  <- ensure_dose_columns(df)  # guarantees dose_log10 and dose_factor
    fml <- stats::formula(base_m) |> collapse_formula()
    
    # 1) Apply requested encoding with strict guards
    has_dose_factor <- "dose_factor" %in% names(df)
    has_dose_log10  <- "dose_log10"  %in% names(df)
    
    want_factor <- isTRUE(input$reg_use_dose_factor)
    want_spline <- isTRUE(input$reg_use_spline_dose)
    
    # Factor takes precedence if both toggles are on
    if (want_factor && want_spline) {
      want_spline <- FALSE
      showNotification("Both dose adjustments selected; using dose_factor and skipping spline.", type = "message", duration = 5)
    }
    
    # Only allow spline if dose_log10 truly exists after ensure_dose_columns()
    if (want_spline && !has_dose_log10) {
      want_spline <- FALSE
      showNotification("Spline skipped: dose_log10 not available in refit data.", type = "warning", duration = 6)
    }
    
    # String form to test for existing terms
    f_chr <- paste(base::deparse(fml), collapse = " ")
    
    if (want_factor && has_dose_factor) {
      # Replace any existing numeric dose term with factor
      if (grepl("splines::ns\\s*\\(\\s*dose_log10\\s*,\\s*2\\s*\\)", f_chr)) {
        fml <- update(fml, . ~ . - splines::ns(dose_log10, 2) + dose_factor)
      } else if (grepl("\\bdose_log10\\b", f_chr)) {
        fml <- update(fml, . ~ . - dose_log10 + dose_factor)
      } else if (!grepl("\\bdose_factor\\b", f_chr)) {
        fml <- update(fml, . ~ . + dose_factor)
      }
    } else if (want_spline && has_dose_log10) {
      # Replace plain dose terms with spline, or add if absent
      if (grepl("\\bdose_factor\\b", f_chr)) {
        fml <- update(fml, . ~ . - dose_factor + splines::ns(dose_log10, 2))
      } else if (grepl("\\bdose_log10\\b", f_chr)) {
        fml <- update(fml, . ~ . - dose_log10 + splines::ns(dose_log10, 2))
      } else if (!grepl("splines::ns\\s*\\(\\s*dose_log10\\s*,\\s*2\\s*\\)", f_chr)) {
        fml <- update(fml, . ~ . + splines::ns(dose_log10, 2))
      }
    }
    # Otherwise keep the base model’s encoding as-is.
    
    # ---- 3) Build model frame FIRST (rows used in refit) ----
    mf <- stats::model.frame(fml, data = df, na.action = stats::na.exclude)
    
    # ---- 4) Weighted least squares on the model frame ----
    use_wls <- FALSE
    if (isTRUE(input$reg_use_wls)) {
      # Grouping priority for heteroskedasticity modeling
      if ("dose_factor" %in% names(mf)) {
        grp <- mf$dose_factor
      } else if ("dose_log10" %in% names(mf)) {
        qs <- stats::quantile(mf$dose_log10, probs = seq(0, 1, 0.25), na.rm = TRUE)
        qs <- unique(qs)
        grp <- if (length(qs) >= 2L) cut(mf$dose_log10, breaks = qs, include.lowest = TRUE, right = TRUE)
        else factor(rep("all", nrow(mf)))
      } else {
        # Fallback: bins of base predictions on the same rows
        y_hat_base <- tryCatch(stats::predict(base_m, newdata = mf), error = function(e) rep(NA_real_, nrow(mf)))
        qs <- stats::quantile(y_hat_base, probs = seq(0, 1, 0.25), na.rm = TRUE)
        qs <- unique(qs)
        grp <- if (length(qs) >= 2L && any(is.finite(y_hat_base)))
          cut(y_hat_base, breaks = qs, include.lowest = TRUE, right = TRUE)
        else factor(rep("all", nrow(mf)))
      }
      
      # Proxy residuals from base model aligned to mf
      y_mf  <- stats::model.response(mf)
      y_hat <- tryCatch(stats::predict(base_m, newdata = mf), error = function(e) rep(NA_real_, length(y_mf)))
      res2  <- (y_mf - y_hat)^2
      
      # Group variances -> inverse-variance weights
      v_by_g <- tapply(res2, grp, function(z) mean(z, na.rm = TRUE))
      v_by_g <- ifelse(is.finite(v_by_g), v_by_g, median(v_by_g, na.rm = TRUE))
      wts    <- 1 / v_by_g[as.character(grp)]
      
      if (length(wts) == nrow(mf) && all(is.finite(wts))) {
        mf$._wts <- as.numeric(wts) / mean(wts, na.rm = TRUE)
        use_wls  <- TRUE
      } else {
        showNotification("WLS disabled: could not construct finite weights on the model frame.", type = "warning", duration = 6)
      }
    }
    
    # ---- 5) Refit on the model frame (._wts attached when needed) ----
    adj_fit <- if (use_wls) stats::lm(fml, data = mf, weights = ._wts) else stats::lm(fml, data = mf)
    
    # ---- 6) Tidy outputs ----
    coef_tidy_base <- broom::tidy(base_m, conf.int = TRUE)
    coef_tidy_adj  <- if (isTRUE(input$reg_use_hc3)) {
      broom::tidy(adj_fit, conf.int = TRUE, vcov = sandwich::vcovHC(adj_fit, type = "HC3"))
    } else {
      broom::tidy(adj_fit, conf.int = TRUE)
    }
    
    coef_cmp <- dplyr::full_join(
      coef_tidy_base %>% dplyr::select(term, estimate, std.error, p.value) %>%
        dplyr::rename(base_estimate = estimate, base_se = std.error, base_p = p.value),
      coef_tidy_adj %>% dplyr::select(term, estimate, std.error, p.value) %>%
        dplyr::rename(adj_estimate = estimate, adj_se = std.error, adj_p = p.value),
      by = "term"
    )
    
    # Fit metrics comparison (patched)
    ll_base <- as.numeric(stats::logLik(base_m))     # extract numeric logLik
    ll_adj  <- as.numeric(stats::logLik(adj_fit))    # extract numeric logLik
    
    cmp <- data.frame(
      metric   = c("AIC","BIC","logLik","sigma","adj_R2","nobs"),
      base     = c(AIC(base_m), BIC(base_m), ll_base,
                   sigma(base_m), summary(base_m)$adj.r.squared, stats::nobs(base_m)),
      adjusted = c(AIC(adj_fit), BIC(adj_fit), ll_adj,
                   sigma(adj_fit), summary(adj_fit)$adj.r.squared, stats::nobs(adj_fit))
    )
    cmp$delta <- cmp$adjusted - cmp$base
    
    list(
      model          = adj_fit,
      coef_cmp       = coef_cmp,
      used_weights   = isTRUE(use_wls),
      weight_group   = if (isTRUE(use_wls)) levels(grp) else NULL,
      comparison     = cmp,
      formula_final  = collapse_formula(fml)
    )
  })
  
  output$regression_model_summary_adjusted <- renderPrint({
    fit <- regression_model_adjusted()
    validate(need(!is.null(fit), "Select options and click 'Refit Adjusted Model'."))
    cat("Adjusted model formula:\n", format(fit$formula), "\n\n")
    print(summary(fit$model))
    if (!is.null(fit$hc3)) {
      cat("\nRobust HC3 coefficient table:\n")
      print(fit$hc3)
    }
    if (isTRUE(fit$used_weights)) {
      cat("\nNote: Weighted least squares was applied using dose_factor-level residual variances.\n")
    }
  })

  # ============================================================================
  #                          MIXED EFFECTS OUTPUTS
  # ============================================================================

  output$lmer_model_summary <- renderPrint({
    fit <- lmer_model()
    validate(need(!inherits(fit, "list"), "Model did not fit."))
    print(summary(fit))
    cat("\nRandom-effects singular? ", lme4::isSingular(fit, tol = 1e-4), "\n")
  })

  # ---- Mixed effects ANOVA (robust DT) ----
  output$lmer_anova <- DT::renderDataTable({
    fit <- lmer_model()
    validate(need(inherits(fit, "lmerMod"), "Fit the mixed model first."))
    
    # Type-III Wald Chi-square ANOVA for merMod
    a3 <- car::Anova(fit, type = 3)
    a3 <- as.data.frame(a3)
    a3 <- tibble::rownames_to_column(a3, "term")
    
    # Normalize potential column name differences across car versions
    names(a3) <- sub("^Pr\\(>Chisq\\)$", "PrChisq", names(a3))
    df_col  <- intersect(c("Df", "df"), names(a3))[1]
    chi_col <- intersect(c("Chisq", "LR Chisq", "Wald", "Wald Chisq"), names(a3))[1]
    
    # Drop rows with non-positive df or non-finite test statistic (prevents blank rows)
    if (!is.null(df_col) && !is.null(chi_col)) {
      a3 <- a3 %>%
        dplyr::filter((!!rlang::sym(df_col)) > 0, is.finite((!!rlang::sym(chi_col))))
    }
    
    # Render with formatting
    DT::datatable(
      a3,
      options = list(pageLength = 10, scrollX = TRUE),
      filter  = "top",
      rownames = FALSE
    ) %>%
      (function(dtbl) {
        cols_to_round <- c(chi_col, "PrChisq")
        cols_to_round <- cols_to_round[cols_to_round %in% names(a3)]
        if (length(cols_to_round)) {
          DT::formatRound(dtbl, columns = cols_to_round, digits = 4)
        } else {
          dtbl
        }
      })
  })

  # Model comparison output
  output$model_compare <- DT::renderDataTable({
    # Simple model comparison between lm and lmer
    lm_model <- regression_model()
    lmer_model_obj <- lmer_model()

    validate(need(
      !is.null(lm_model) && !inherits(lm_model, "list") &&
        !is.null(lmer_model_obj) && !inherits(lmer_model_obj, "list"),
      "Run both models to compare"
    ))

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

    DT::datatable(comparison_df, options = list(pageLength = 5, dom = "t"), rownames = FALSE) %>%
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
    wks <- tryCatch(
      {
        if (identical(ds, "assay")) {
          req(exists("final_data"))
          sort(unique(suppressWarnings(as.integer(final_data$week[final_data$assay_type == ep]))))
        } else {
          req(exists("physical_master"))
          sort(unique(suppressWarnings(as.integer(physical_master$week[physical_master$endpoint == ep]))))
        }
      },
      error = function(e) integer()
    )

    # Fallback if none found
    if (length(wks) == 0L || all(is.na(wks))) wks <- c(1L, 5L, 6L)

    # Update both selectors with character choices so Shiny can match selected
    updateSelectInput(session, "recovery_baseline_week",
      choices = as.character(wks),
      selected = as.character(min(wks, na.rm = TRUE))
    )
    updateSelectInput(session, "recovery_recovery_week",
      choices = as.character(wks),
      selected = as.character(max(wks, na.rm = TRUE))
    )
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
      if (!is.finite(y)) {
        return(NA_real_)
      }
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
        {
          d <- .
          # Derive label-like sample_type only from sample_type/tissue_type/tissue
          raw_lab <- coalesce_first_existing(d, c("sample_type", "tissue_type", "tissue"))
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
    if (!is.null(input$recovery_fiber_types) && length(input$recovery_fiber_types)) {
      df <- dplyr::filter(df, fiber_type %in% input$recovery_fiber_types)
    }

    if (!is.null(input$recovery_treatments) && length(input$recovery_treatments)) {
      df <- dplyr::filter(df, chem_treatment %in% input$recovery_treatments)
    }

    # 3) Sample/tissue filters with guards
    if (identical(ds, "assay")) {
      if (!is.null(input$recovery_samples) && length(input$recovery_samples)) {
        df <- dplyr::filter(df, tolower(sample_type) %in% tolower(input$recovery_samples))
      }
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
        if (!is.null(input$recovery_concentration) && length(input$recovery_concentration)) {
          df <- dplyr::filter(df, fiber_concentration_label %in% input$recovery_concentration)
        }
      } else {
        if (!is.null(input$recovery_concentration_physical) && length(input$recovery_concentration_physical)) {
          df <- dplyr::filter(df, fiber_concentration_label %in% input$recovery_concentration_physical)
        }
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
          levels = c(lbl_base, lbl_reco)
        ),
        fiber_type = factor(fiber_type),
        chem_treatment = factor(chem_treatment),
        sample_type = factor(sample_type)
      ) %>%
      droplevels()

    validate(need(dplyr::n_distinct(df$is_recovery) == 2, "Both selected weeks must be present after filters"))

    # 6) Hierarchical week coverage (unchanged logic)
    present <- function(v) v %in% names(df) && nlev2(df[[v]]) > 1
    strict_vars <- c(
      if (present("fiber_type")) "fiber_type",
      if (present("chem_treatment")) "chem_treatment",
      if (present("sample_type") && identical(ds, "assay")) "sample_type"
    )
    candidates <- list(
      strict_vars,
      setdiff(strict_vars, "sample_type"),
      setdiff(strict_vars, c("sample_type", "fiber_type")),
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
        if (dplyr::n_distinct(df$is_recovery) == 2) df else df[0, ]
      }
      if (nrow(df_try) > 0) {
        picked <- gv
        df <- df_try
        break
      }
    }
    msg_pick <- if (is.null(picked) || length(picked) == 0) "pooled" else paste(picked, collapse = "×")
    message("[RECOVERY] Week coverage by: ", msg_pick, " | rows: ", nrow(df))
    validate(need(nrow(df) > 0, "No strata retain both weeks; widen filters or include another dose"))

    # 7) Estimable formula (keep existing structure)
    include_trt <- nlev2(df$chem_treatment) >= 2
    include_fib <- nlev2(df$fiber_type) >= 2
    include_samp <- nlev2(df$sample_type) >= 2

    rhs <- c("is_recovery")
    if (include_trt) rhs <- c(rhs, "chem_treatment", "is_recovery:chem_treatment")
    if (include_fib) rhs <- c(rhs, "fiber_type", "is_recovery:fiber_type")
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
      ~is_recovery
    }

    emm <- tryCatch(emmeans::emmeans(mdl, emm_spec), error = function(e) {
      message("[RECOVERY] emmeans failed: ", e$message)
      NULL
    })

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
      ggplot2::stat_summary(fun.data = ~ c(
        y = mean(.), ymin = mean(.) - sd(.) / sqrt(length(.)),
        ymax = mean(.) + sd(.) / sqrt(length(.))
      ), geom = "errorbar", width = 0.2) +
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
    validate(need(!is.null(res), "")) # ensure run
    validate(need(is.null(res$error), res$error)) # show model errors if any
    validate(need(
      !is.null(res$emmeans$delta),
      "No Week6 vs Week5 contrasts available for current filters."
    ))

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
        p_value_raw = p.value, # keep raw p from emmeans table
        p_value_adj = suppressWarnings(p.adjust(p_value_raw, method = adj_method)),
        stars = dplyr::case_when( # star coding on adjusted p
          is.na(p_value_adj) ~ "",
          p_value_adj < 0.001 ~ "***",
          p_value_adj < 0.01 ~ "**",
          p_value_adj < 0.05 ~ "*",
          TRUE ~ ""
        ),
        Significant = dplyr::if_else(!is.na(p_value_adj) & p_value_adj < 0.05,
          paste0("Yes ", stars),
          paste0("No ", stars)
        )
      )

    # Friendlier label if the emmeans 'contrast' column exists
    if ("contrast" %in% names(df)) {
      df <- df %>% dplyr::rename(`Week6 vs Week5` = contrast)
    }

    # Columns to display (only those that exist)
    show_cols <- intersect(
      c(
        "Week6 vs Week5", "is_recovery", "chem_treatment", "fiber_type",
        "estimate", "SE", "df", "t.ratio", "z.ratio",
        "p_value_raw", "p_value_adj", "lower.CL", "upper.CL", "Significant"
      ),
      names(df)
    )
    out <- df %>% dplyr::select(dplyr::all_of(show_cols))

    # Value-specific styling for the Significant column (green for Yes, bold)
    unique_vals <- unique(out$Significant)
    color_map <- ifelse(grepl("^Yes", unique_vals), "#1a7f37", "#444444")
    weight_map <- ifelse(grepl("^Yes", unique_vals), "bold", "normal")

    DT::datatable(
      out,
      rownames = FALSE,
      filter   = "top",
      options  = list(pageLength = 12, scrollX = TRUE, dom = "Blfrtip")
    ) %>%
      {
        if ("estimate" %in% names(out)) DT::formatRound(., "estimate", 4) else .
      } %>%
      {
        if ("SE" %in% names(out)) DT::formatRound(., "SE", 4) else .
      } %>%
      {
        if ("t.ratio" %in% names(out)) DT::formatRound(., "t.ratio", 4) else .
      } %>%
      {
        if ("z.ratio" %in% names(out)) DT::formatRound(., "z.ratio", 4) else .
      } %>%
      {
        if ("p_value_raw" %in% names(out)) DT::formatRound(., "p_value_raw", 6) else .
      } %>%
      {
        if ("p_value_adj" %in% names(out)) DT::formatRound(., "p_value_adj", 6) else .
      } %>%
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

  # helper: normalize categorical strings before grouping/factoring
  normalize_categories <- function(df) {
    df %>%
      dplyr::mutate(
        fiber_type = stringr::str_trim(stringr::str_to_lower(fiber_type)),
        treatment = stringr::str_trim(stringr::str_to_lower(treatment)),
        tissue_type = stringr::str_trim(stringr::str_to_lower(tissue_type)),
        fiber_concentration = as.character(fiber_concentration) # keep as ordered labels
      )
  }

  # ---- Single-week, two-tissue, multi-dose reactive (fixed) ----
  transloc_filtered <- eventReactive(input$run_translocation, {
    req(physical_master, nrow(physical_master) > 0)

    validate(need(!is.null(input$transloc_week), "Select a week for translocation."))
    validate(need(length(input$transloc_tissues) == 2, "Select exactly two tissues to compare."))
    validate(need(length(input$transloc_conc) >= 1, "Select at least one fiber concentration."))

    fibers_sel <- if (!is.null(input$transloc_filter_fibers) && length(input$transloc_filter_fibers) > 0) {
      tolower(input$transloc_filter_fibers)
    } else {
      c("cotton", "pet")
    }
    trts_sel <- if (!is.null(input$transloc_filter_treatments) && length(input$transloc_filter_treatments) > 0) {
      tolower(input$transloc_filter_treatments)
    } else {
      c("treated", "untreated")
    }

    df <- physical_master %>%
      dplyr::filter(
        endpoint == "mf_counts",
        week == as.integer(input$transloc_week),
        fiber_concentration %in% input$transloc_conc,
        tissue_type %in% input$transloc_tissues,
        tolower(fiber_type) %in% fibers_sel,
        tolower(treatment) %in% trts_sel
      ) %>%
      normalize_categories() %>% # <- NEW
      dplyr::group_by(fiber_type, treatment, fiber_concentration, week, tissue_type) %>%
      dplyr::summarise(
        mean_count = mean(value, na.rm = TRUE),
        sd_count   = sd(value, na.rm = TRUE),
        n          = dplyr::n(),
        .groups    = "drop"
      ) %>%
      dplyr::mutate(
        fiber_type = factor(fiber_type, levels = c("cotton", "pet")),
        treatment  = factor(treatment, levels = c("treated", "untreated"))
      )

    validate(need(nrow(df) > 0, "No mf_counts data for the selected week/tissues/concentrations/filters."))
    df
  })

  observeEvent(list(input$transloc_facet_rows, input$transloc_facet_cols),
    {
      r <- input$transloc_facet_rows %||% "fiber_type"
      c <- input$transloc_facet_cols %||% "treatment"
      if (!identical(r, "none") && identical(r, c)) {
        # Prefer keeping rows and clear columns
        updateSelectInput(session, "transloc_facet_cols", selected = "none")
      }
    },
    ignoreInit = TRUE
  )

  # ---- Tissue difference summary (fixed) ----
  transloc_summary <- reactive({
    df <- transloc_filtered()

    # names of the two tissues, normalized
    ta <- stringr::str_to_lower(input$transloc_tissues[1])
    tb <- stringr::str_to_lower(input$transloc_tissues[2])

    # Make sure both tissues exist post-filter
    available_tissues <- unique(df$tissue_type)
    validate(need(
      all(c(ta, tb) %in% available_tissues),
      paste0(
        "Difference plot unavailable: selected facet filters exclude one or both tissues (",
        ta, ", ", tb, ")."
      )
    ))

    # Ensure one row per stratum x tissue; average if anything slipped through
    df_uni <- df %>%
      dplyr::group_by(fiber_type, treatment, fiber_concentration, week, tissue_type) %>%
      dplyr::summarise(mean_count = mean(as.numeric(mean_count), na.rm = TRUE), .groups = "drop")

    wide <- df_uni %>%
      tidyr::pivot_wider(
        id_cols = c(fiber_type, treatment, fiber_concentration, week),
        names_from = tissue_type,
        values_from = mean_count,
        values_fn = mean, # final guard; keeps numeric if rare duplicates remain
        values_fill = NA_real_
      )

    # Compute A - B
    wide %>%
      dplyr::mutate(
        diff_ab  = .data[[ta]] - .data[[tb]],
        contrast = paste0(ta, " - ", tb)
      ) %>%
      dplyr::filter(!is.na(diff_ab))
  })

  # ---- Multi-week, single-tissue, multi-dose reactive (fixed) ----
  transloc_multi_filtered <- eventReactive(input$run_translocation, {
    req(physical_master, nrow(physical_master) > 0)

    validate(need(length(input$transloc_weeks_multi) >= 1, "Select one or more weeks."))
    validate(need(
      !is.null(input$transloc_tissue_single) && nzchar(input$transloc_tissue_single),
      "Select a single tissue."
    ))
    validate(need(length(input$transloc_conc_multi) >= 1, "Select at least one concentration."))

    fibers_sel <- if (!is.null(input$transloc_filter_fibers) && length(input$transloc_filter_fibers) > 0) {
      tolower(input$transloc_filter_fibers)
    } else {
      c("cotton", "pet")
    }
    trts_sel <- if (!is.null(input$transloc_filter_treatments) && length(input$transloc_filter_treatments) > 0) {
      tolower(input$transloc_filter_treatments)
    } else {
      c("treated", "untreated")
    }

    df <- physical_master %>%
      dplyr::filter(
        endpoint == "mf_counts",
        week %in% as.integer(input$transloc_weeks_multi),
        fiber_concentration %in% input$transloc_conc_multi,
        tissue_type == input$transloc_tissue_single,
        tolower(fiber_type) %in% fibers_sel,
        tolower(treatment) %in% trts_sel
      ) %>%
      normalize_categories() %>% # <- NEW
      dplyr::group_by(fiber_type, treatment, fiber_concentration, week) %>%
      dplyr::summarise(
        mean_count = mean(value, na.rm = TRUE),
        sd_count   = sd(value, na.rm = TRUE),
        n          = dplyr::n(),
        .groups    = "drop"
      ) %>%
      dplyr::mutate(
        fiber_concentration = factor(fiber_concentration, levels = c("100", "1000", "10000")),
        week = factor(week, levels = sort(unique(week))),
        fiber_type = factor(fiber_type, levels = c("cotton", "pet")),
        treatment = factor(treatment, levels = c("treated", "untreated"))
      )

    validate(need(nrow(df) > 0, "No mf_counts data for the selected multi-week filters."))
    df
  })


  # ----------------- Translocation counts plot (dynamic facets) -----------------
  output$transloc_counts_plot <- renderPlot({
    df <- transloc_filtered() # existing eventReactive for mf_counts
    validate(need(nrow(df) > 0, "No data to plot."))

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = fiber_concentration, y = mean_count, fill = tissue_type)
    ) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = pmax(mean_count - sd_count, 0), ymax = mean_count + sd_count),
        position = ggplot2::position_dodge(width = 0.75),
        width = 0.2
      ) +
      ggplot2::labs(
        title = paste0("Microfiber counts by tissue (Week ", input$transloc_week, ")"),
        x = "Fiber concentration (mf/L)", y = "Mean count", fill = "Tissue"
      ) +
      ggplot2::theme_minimal(base_size = 12)

    # Apply chosen facets
    rows <- input$transloc_facet_rows %||% "fiber_type"
    cols <- input$transloc_facet_cols %||% "treatment"
    apply_dynamic_facets(p, rows = rows, cols = cols)
  })


  # ----------------- Translocation difference plot (dynamic facets) -----------------
  output$transloc_diff_plot <- renderPlot({
    # Try to compute summary; if it fails, skip silently
    smry <- tryCatch(transloc_summary(), error = function(e) NULL)
    if (is.null(smry) || nrow(smry) == 0) {
      validate(need(FALSE, "Difference plot unavailable: selected facet filters exclude one or both tissues."))
    }

    p <- ggplot2::ggplot(
      smry,
      ggplot2::aes(x = fiber_concentration, y = diff_ab, fill = treatment)
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
      ggplot2::labs(
        title = paste0("Tissue difference (", unique(smry$contrast), ") at Week ", input$transloc_week),
        x = "Fiber concentration (mf/L)",
        y = "Mean count difference (A − B)",
        fill = "Treatment"
      ) +
      ggplot2::theme_minimal(base_size = 12)

    rows <- input$transloc_facet_rows %||% "fiber_type"
    cols <- input$transloc_facet_cols %||% "treatment"
    apply_dynamic_facets(p, rows = rows, cols = cols)
  })

  # ---- Multi-week grouped bar plot: x = concentration, fill = week, dynamic faceting ----
  output$transloc_multiweek_plot <- renderPlot({
    # Existing filtered summary for multi-week, single-tissue, multi-dose
    df <- transloc_multi_filtered() # keep your current reactive name

    # Defensive guard (matches pattern used elsewhere)
    validate(need(!is.null(df) && nrow(df) > 0, "No mf_counts data for selected filters."))

    # Base plot stays identical to your current implementation
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = fiber_concentration, y = mean_count, fill = week)
    ) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = pmax(mean_count - sd_count, 0), ymax = mean_count + sd_count),
        position = ggplot2::position_dodge(width = 0.75),
        width = 0.2
      ) +
      ggplot2::labs(
        title = paste0("Translocation by concentration across weeks (", input$transloc_tissue_single, ")"),
        x = "Fiber concentration (mf/L)",
        y = "Mean count",
        fill = "Week"
      ) +
      ggplot2::theme_minimal(base_size = 12)

    # NEW: dynamic facets using the same pickers as the first plot
    # Defaults mirror your current layout (rows = fiber_type, cols = treatment)
    rows <- input$transloc_facet_rows %||% "fiber_type"
    cols <- input$transloc_facet_cols %||% "treatment"

    # Reuse your helper introduced for the first plot
    apply_dynamic_facets(p, rows = rows, cols = cols)
  })

  # Keep your existing download handlers unchanged
  output$download_transloc_multi_png <- create_plot_download_handler(reactive(output$transloc_multiweek_plot()), "transloc_multiweek", "png")
  output$download_transloc_multi_pdf <- create_plot_download_handler(reactive(output$transloc_multiweek_plot()), "transloc_multiweek", "pdf")
  output$download_transloc_multi_tiff <- create_plot_download_handler(reactive(output$transloc_multiweek_plot()), "transloc_multiweek", "tiff")
  output$download_transloc_multi_svg <- create_plot_download_handler(reactive(output$transloc_multiweek_plot()), "transloc_multiweek", "svg")


  # ----------------- Translocation table -----------------

  output$transloc_summary_table <- DT::renderDataTable({
    smry <- transloc_summary()
    validate(need(nrow(smry) > 0, "No summary to show."))

    # Tidy table for inspection/download
    out <- smry %>%
      dplyr::arrange(fiber_type, treatment, fiber_concentration) %>%
      dplyr::select(
        fiber_type, treatment, week, fiber_concentration,
        contrast, dplyr::all_of(input$transloc_tissues), diff_ab
      )

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
        shinyjs::enable("did_samples")
        shinyjs::enable("did_tissues")
        output$did_samples_hint <- renderUI(NULL)
        output$did_tissues_hint <- renderUI(NULL)
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
          outcome = calculated_concentration,
          week = suppressWarnings(as.integer(week)),
          fiber_type = did_canon_fiber(fiber_type),
          chem_treatment = did_canon_treatment(treatment),
          sample_type = tolower(trimws(as.character(sample_type))),
          fiber_concentration_label = did_canon_conc_labels(fiber_concentration)
        )
    } else {
      validate(need(exists("physical_master") && nrow(physical_master) > 0, "No physical data loaded"))
      physical_master %>%
        dplyr::filter(.data$endpoint == input$did_endpoint) %>%
        {
          d <- .
          # Use only label-like columns; do not pull from numeric 'sample' ids
          raw_lab <- coalesce_first_existing(d, c("sample_type", "tissue_type", "tissue"))
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
          outcome = value,
          week = suppressWarnings(as.integer(week)),
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
      df <- normalize_controls_and_dose(df) # keep if defined
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

    ds <- input$did_dataset
    ep <- input$did_endpoint
    mode <- input$did_mode
    if (is.null(ds) || is.null(ep) || is.null(mode)) {
      return(safe_abort("Select dataset, endpoint, and mode before running DiD."))
    }

    # Base (un-subset) and pre-filtered (sample/tissue + dose) frames
    df0 <- did_base_data_norm()
    df <- did_data_filtered()
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
    if (nrow(df) == 0) {
      return(safe_abort("No data after optional subsets."))
    }

    # -----------------------------
    # Mode-specific subsetting
    # -----------------------------
    if (identical(mode, "across_fiber")) {
      wk <- safe_int1(input$did_week_select)
      if (is.na(wk)) {
        return(safe_abort("Select a week for Across Fiber DiD."))
      }
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
            chem_treatment %in% c("untreated", "treated"),
            fiber_type %in% c(input$did_fiber_a, input$did_fiber_b),
            fiber_concentration_label %in% c("0", dose_chosen)
          ) %>%
          dplyr::mutate(
            fiber_type = factor(fiber_type, levels = c("cotton", "pet")),
            chem_treatment = factor(chem_treatment, levels = c("untreated", "treated")),
            fiber_concentration_label = factor(fiber_concentration_label, levels = c("0", dose_chosen)),
            sample_type = factor(sample_type)
          ) %>%
          droplevels()

        has_two <- function(v) v %in% names(df) && dplyr::n_distinct(stats::na.omit(df[[v]])) >= 2
        if (!all(c(has_two("fiber_type"), has_two("chem_treatment"), has_two("fiber_concentration_label")))) {
          return(safe_abort("Need both fibers, both treatments, and doses {0 and selected} for dose-based DiD."))
        }

        rhs <- "chem_treatment * fiber_type * fiber_concentration_label"
        if (has_two("sample_type")) rhs <- paste(rhs, "+ sample_type")
        form <- stats::as.formula(paste("outcome ~", rhs))
        mdl <- tryCatch(stats::lm(form, data = df), error = function(e) e)
        if (inherits(mdl, "error")) {
          return(safe_abort(paste("Model failed:", mdl$message)))
        }

        emm3 <- safe_emmeans(mdl, ~ chem_treatment * fiber_type * fiber_concentration_label)
        emm_tbl <- if (!is.null(emm3)) tryCatch(as.data.frame(emm3), error = function(e) NULL) else NULL
        if (is.null(emm_tbl) || nrow(emm_tbl) == 0) {
          return(safe_abort("EMMeans grid is empty for dose-based DiD."))
        }

        a <- input$did_fiber_a
        b <- input$did_fiber_b
        need <- rbind(
          c("treated", a, dose_chosen),
          c("untreated", a, "0"),
          c("treated", b, dose_chosen),
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
        w[sel("treated", a, dose_chosen)] <- 1
        w[sel("untreated", a, "0")] <- -1
        w[sel("treated", b, dose_chosen)] <- -1
        w[sel("untreated", b, "0")] <- 1

        did <- if (!is.null(emm3)) tryCatch(emmeans::contrast(emm3, method = list(did_dose = w), by = NULL), error = function(e) NULL) else NULL
        did_tbl <- if (!is.null(did)) tryCatch(as.data.frame(did), error = function(e) NULL) else NULL
        if (!is.null(did_tbl) && "p.value" %in% names(did_tbl)) {
          method <- input$did_p_adjust %||% "none"
          did_tbl$p_adjust_method <- method
          did_tbl$p_adj <- stats::p.adjust(did_tbl$p.value, method = method)
        }

        did_state(list(
          model = mdl, data = df, emmeans_tbl = emm_tbl, did_tbl = did_tbl,
          mode = "across_fiber", baseline = "dose"
        ))
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
            fiber_type = factor(fiber_type, levels = c("cotton", "pet")),
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
        mdl <- tryCatch(stats::lm(form, data = df_uc), error = function(e) e)
        if (inherits(mdl, "error")) {
          return(safe_abort(paste("Model failed:", mdl$message)))
        }

        emm2 <- safe_emmeans(mdl, ~ fiber_type * fiber_concentration_label)
        emm_tbl <- if (!is.null(emm2)) tryCatch(as.data.frame(emm2), error = function(e) NULL) else NULL
        if (is.null(emm_tbl) || nrow(emm_tbl) == 0) {
          return(safe_abort("EMMeans grid is empty for untreated control-offset DiD."))
        }

        a <- input$did_fiber_a
        b <- input$did_fiber_b
        have_key <- paste(emm_tbl$fiber_type, emm_tbl$fiber_concentration_label, sep = ".")
        need_key <- c(
          paste(a, dose_chosen, sep = "."), paste(a, "0", sep = "."),
          paste(b, dose_chosen, sep = "."), paste(b, "0", sep = ".")
        )
        if (!all(need_key %in% have_key)) {
          return(safe_abort("Missing untreated cells (fiber × {0, dose})."))
        }

        # Weights: [untreated(dose) - untreated(0)]_A - [untreated(dose) - untreated(0)]_B
        w <- rep(0, length(have_key))
        w[have_key == paste(a, dose_chosen, sep = ".")] <- 1
        w[have_key == paste(a, "0", sep = ".")] <- -1
        w[have_key == paste(b, dose_chosen, sep = ".")] <- -1
        w[have_key == paste(b, "0", sep = ".")] <- 1

        did_uc <- if (!is.null(emm2)) {
          tryCatch(
            emmeans::contrast(emm2, method = list(untreated_control_offset = w), by = NULL),
            error = function(e) NULL
          )
        } else {
          NULL
        }
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
      wk_a <- safe_int1(input$did_week_a)
      wk_b <- safe_int1(input$did_week_b)
      if (anyNA(c(wk_a, wk_b))) {
        return(safe_abort("Select both Week A and Week B."))
      }
      df <- df %>%
        dplyr::filter(
          week %in% c(wk_a, wk_b),
          fiber_type == input$did_fiber_select
        )
      message("[DiD] Across_week rows at weeks ", paste(c(wk_a, wk_b), collapse = ", "), ": ", nrow(df))
    } else { # across_treatment
      wk <- safe_int1(input$did_week_select)
      if (is.na(wk)) {
        return(safe_abort("Select a week for Across Treatment DiD."))
      }
      df <- df %>%
        dplyr::filter(
          week == wk,
          fiber_type == input$did_fiber_select,
          chem_treatment %in% c(input$did_trt_a, input$did_trt_b)
        )
      message("[DiD] Across_treatment rows at week ", wk, " within fiber ", input$did_fiber_select, ": ", nrow(df))
    }

    if (nrow(df) == 0) {
      return(safe_abort("No rows after mode-specific subsetting."))
    }

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
        df <- df %>%
          dplyr::filter(fiber_concentration_label %in% input$did_concentration) %>%
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
        wk_a <- safe_int1(input$did_week_a)
        wk_b <- safe_int1(input$did_week_b)
        df <- df %>% dplyr::filter(week %in% c(wk_a, wk_b), fiber_type == input$did_fiber_select)
      } else { # across_treatment
        wk <- safe_int1(input$did_week_select)
        df <- df %>%
          dplyr::filter(
            week == wk, fiber_type == input$did_fiber_select,
            chem_treatment %in% c(input$did_trt_a, input$did_trt_b)
          )
      }

      coverage_ok <- (!needs_both_fib || has_both("fiber_type")) && (!needs_both_trt || has_both("chem_treatment"))
      if (!coverage_ok) {
        return(safe_abort("Not enough data to estimate the requested contrast; widen filters."))
      }
      message("[DiD] Coverage restored after pooling; rows: ", nrow(df))
    }

    # Standardize factors (drop singletons)
    df <- df %>%
      dplyr::mutate(
        week = factor(week),
        fiber_type = factor(fiber_type, levels = c("cotton", "pet")),
        chem_treatment = factor(chem_treatment, levels = c("untreated", "treated")),
        sample_type = factor(sample_type)
      ) %>%
      droplevels()

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
    } else { # across_treatment
      if (!has_varied("chem_treatment")) {
        return(safe_abort("Need both selected treatments present for Across Treatment within this fiber."))
      }
      rhs <- "chem_treatment"
      if (has_varied("sample_type")) rhs <- paste(rhs, "+ sample_type")
    }
    form <- stats::as.formula(paste("outcome ~", rhs))

    # Fit model safely
    mdl <- tryCatch(stats::lm(form, data = df), error = function(e) e)
    if (inherits(mdl, "error")) {
      return(safe_abort(paste("Model failed:", mdl$message)))
    }
    message("[DiD] Model fit OK; n=", nrow(df))

    # Defensive EMMeans and contrasts for each mode
    emm_tbl <- did_tbl <- NULL

    if (identical(mode, "across_fiber")) {
      emm2 <- safe_emmeans(mdl, ~ chem_treatment * fiber_type)
      emm_tbl <- if (!is.null(emm2)) tryCatch(as.data.frame(emm2), error = function(e) NULL) else NULL

      fibers_needed <- unique(c(input$did_fiber_a, input$did_fiber_b))
      trts_needed <- c("untreated", "treated")
      have_cells <- if (!is.null(emm_tbl)) paste(emm_tbl$chem_treatment, emm_tbl$fiber_type, sep = ".") else character(0)
      need_cells <- as.vector(outer(trts_needed, fibers_needed, paste, sep = "."))
      if (!all(need_cells %in% have_cells)) {
        return(safe_abort("Selected week/filters are missing treatment-by-fiber cells; widen filters."))
      }

      a <- input$did_fiber_a
      b <- input$did_fiber_b
      w <- rep(0, length(have_cells))
      w[have_cells == paste("treated", a, sep = ".")] <- 1
      w[have_cells == paste("untreated", a, sep = ".")] <- -1
      w[have_cells == paste("treated", b, sep = ".")] <- -1
      w[have_cells == paste("untreated", b, sep = ".")] <- 1

      did_tbl <- if (!is.null(emm2)) {
        tryCatch(
          as.data.frame(emmeans::contrast(emm2, method = list(DiD = w), by = NULL)),
          error = function(e) NULL
        )
      } else {
        NULL
      }
    } else if (identical(mode, "across_week")) {
      emm2 <- safe_emmeans(mdl, ~ chem_treatment * week)
      emm_tbl <- if (!is.null(emm2)) tryCatch(as.data.frame(emm2), error = function(e) NULL) else NULL

      wa <- as.character(input$did_week_a)
      wb <- as.character(input$did_week_b)
      trts_needed <- c("untreated", "treated")
      have_cells <- if (!is.null(emm_tbl)) paste(emm_tbl$chem_treatment, emm_tbl$week, sep = ".") else character(0)
      need_cells <- as.vector(outer(trts_needed, c(wa, wb), paste, sep = "."))
      if (!all(need_cells %in% have_cells)) {
        return(safe_abort("Selected weeks are missing treatment cells; widen filters."))
      }

      w <- rep(0, length(have_cells))
      w[have_cells == paste("treated", wb, sep = ".")] <- 1
      w[have_cells == paste("untreated", wb, sep = ".")] <- -1
      w[have_cells == paste("treated", wa, sep = ".")] <- -1
      w[have_cells == paste("untreated", wa, sep = ".")] <- 1

      did_tbl <- if (!is.null(emm2)) {
        tryCatch(
          as.data.frame(emmeans::contrast(emm2, method = list(DiD = w), by = NULL)),
          error = function(e) NULL
        )
      } else {
        NULL
      }
    } else { # across_treatment within selected fiber
      emm2 <- safe_emmeans(mdl, ~chem_treatment)
      emm_tbl <- if (!is.null(emm2)) tryCatch(as.data.frame(emm2), error = function(e) NULL) else NULL

      ta <- input$did_trt_a
      tb <- input$did_trt_b
      have <- if (!is.null(emm_tbl)) emm_tbl$chem_treatment else character(0)
      if (!all(c(ta, tb) %in% have)) {
        return(safe_abort("Selected treatments not both present for the chosen fiber; widen filters."))
      }

      w_named <- stats::setNames(c(1, -1), c(ta, tb)) # TrtA − TrtB
      did_tbl <- if (!is.null(emm2)) {
        tryCatch(
          as.data.frame(emmeans::contrast(emm2, method = list(treatment_effect = w_named), by = NULL)),
          error = function(e) NULL
        )
      } else {
        NULL
      }
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
    if (is.null(res)) {
      cat("Click 'Run DiD Analysis'\n")
      return()
    }
    if (!is.null(res$error)) {
      cat("Error:", res$error, "\n")
      return()
    }
    cat("DiD Mode:", res$mode, "\n")
    print(summary(res$model))
  })

  output$did_table <- DT::renderDataTable({
    res <- did_state()
    if (is.null(res)) {
      return(NULL)
    }
    if (!is.null(res$error)) {
      return(DT::datatable(data.frame(Message = res$error), rownames = FALSE))
    }
    if (is.null(res$did_tbl) || nrow(res$did_tbl) == 0) {
      return(DT::datatable(data.frame(Message = "No DiD rows for current settings."), rownames = FALSE))
    }
    DT::datatable(res$did_tbl, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  output$did_emm_table <- DT::renderDataTable({
    res <- did_state()
    if (is.null(res)) {
      return(NULL)
    }
    if (!is.null(res$error)) {
      return(DT::datatable(data.frame(Message = res$error), rownames = FALSE))
    }
    if (is.null(res$emmeans_tbl) || nrow(res$emmeans_tbl) == 0) {
      return(DT::datatable(data.frame(Message = "No EMMeans grid for current settings."), rownames = FALSE))
    }
    DT::datatable(res$emmeans_tbl, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  output$did_plot <- renderPlot({
    res <- did_state()
    validate(need(!is.null(res), ""))
    dfp <- res$data
    validate(need(!is.null(dfp) && nrow(dfp) > 0, "No data to plot."))
    ggplot2::ggplot(dfp, ggplot2::aes(x = chem_treatment, y = outcome, color = fiber_type)) +
      ggplot2::stat_summary(
        fun = mean, geom = "point", size = 3,
        position = ggplot2::position_dodge(width = 0.35)
      ) +
      ggplot2::stat_summary(
        fun.data = ~ c(
          y = mean(.), ymin = mean(.) - sd(.) / sqrt(length(.)),
          ymax = mean(.) + sd(.) / sqrt(length(.))
        ),
        geom = "errorbar", width = 0.2,
        position = ggplot2::position_dodge(width = 0.35)
      ) +
      ggplot2::facet_wrap(~week, nrow = 1, scales = "free_y") +
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
