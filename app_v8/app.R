# ============================================================================
# FIXED SHINY APP - Corrected bindCache/bindEvent Order  
# ============================================================================

# Source all modules (loaded once at startup)
source("global.R", local = TRUE)      # Libraries, data, config
source("helpers.R", local = TRUE)     # Data processing functions  
source("modeling.R", local = TRUE)    # Statistical modeling functions

# ============================================================================
# USER INTERFACE (Keep your existing UI - no changes needed)
# ============================================================================

ui <- fluidPage(
  useShinyjs(),
  titlePanel("ECOTOX Assay Dashboard"),
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
                       choices = if(exists("config") && !is.null(config$samples)) {
                         tolower(config$samples)
                       } else if(exists("final_data") && nrow(final_data) > 0) {
                         sort(unique(tolower(final_data$sample_type)))
                       } else {
                         c("tissue", "gland", "gills", "hemolymph")
                       }, 
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
                                                  choices = c("1" = "1", "3" = "3", "5" = "5", "6 (recovery)" = "6"),
                                                  selected = c("1","3","5"), inline = TRUE),
                                checkboxInput("include_recovery_interaction", 
                                             "Include recovery × treatment interaction (for week 6)", FALSE),
                                checkboxInput("include_three_way", "Include 3-way: fiber × treatment × week", TRUE),
                                checkboxInput("dose_by_fiber", "Include dose × fiber_type", FALSE),
                                checkboxInput("dose_by_treat", "Include dose × chem_treatment", FALSE),
                                actionButton("run_regression", "Run Regression Model", icon = icon("play"))
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
                                           choices = c("All Pairwise" = "pairwise", "Treatment vs Control" = "control",
                                                      "Difference of Differences (vs control)" = "did"), selected = "pairwise"),
                                conditionalPanel(condition = "input.comparison_type == 'did'", uiOutput("did_controls")),
                                selectInput("adjustment_method", "P-value Adjustment:",
                                           choices = c("Tukey HSD" = "tukey", "Holm" = "holm", 
                                                      "BH (FDR)" = "fdr", "Bonferroni" = "bonferroni"), selected = "tukey"),
                                actionButton("run_emmeans", "Calculate EM Means", icon = icon("chart-line"))
                              )),
                  mainPanel(width = 8,
                           h5("Model Summary"), verbatimTextOutput("regression_model_summary"),
                           h5("ANOVA Table"), DT::dataTableOutput("regression_anova"),
                           h5("Estimated Marginal Means"), plotOutput("emmeans_plot", height = "420px"),
                           DT::dataTableOutput("emmeans_table"),
                           h5("Pairwise Comparisons"), DT::dataTableOutput("pairwise_table")))),
        
        # Mixed Effects Analysis Tab
        tabPanel("Mixed Effects Analysis",
                sidebarLayout(
                  sidebarPanel(width = 4,
                              h4("Linear Mixed Effects Regression"),
                              wellPanel(
                                h5("Model Settings"),
                                checkboxInput("lmer_include_three_way", "Include 3-way: fiber × treatment × week", TRUE),
                                checkboxInput("lmer_dose_by_fiber", "Include dose × fiber_type", FALSE),
                                checkboxInput("lmer_dose_by_treat", "Include dose × chem_treatment", FALSE),
                                checkboxGroupInput("lmer_weeks_include", "Include Weeks:",
                                                  choices = c("1" = "1", "3" = "3", "5" = "5", "6 (recovery)" = "6"),
                                                  selected = c("1","3","5"), inline = TRUE),
                                checkboxInput("lmer_include_recovery_interaction",
                                             "Include recovery × treatment interaction (for week 6)", FALSE),
                                checkboxGroupInput("lmer_random_terms", "Random intercepts:",
                                                  choices = c("Tank" = "tank", "Experiment batch" = "experiment_batch"),
                                                  selected = c("tank","experiment_batch")),
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
                                           choices = c("All Pairwise" = "pairwise", "Treatment vs Control" = "control",
                                                      "Difference of Differences (vs control)" = "did"), selected = "pairwise"),
                                conditionalPanel(condition = "input.lmer_comparison_type == 'did'", uiOutput("lmer_did_controls")),
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
                           h5("Model Comparison (lm vs lmer)"), DT::dataTableOutput("model_compare")))),
        
        # Recovery & Translocation Tab  
        tabPanel("Recovery & Translocation Analysis",
                sidebarLayout(
                  sidebarPanel(width = 4,
                              h4("Recovery Analysis (Week 6 vs Week 5)"),
                              wellPanel(
                                selectInput("recovery_endpoint", "Select Endpoint:", choices = NULL),
                                selectInput("recovery_dataset", "Dataset:",
                                           choices = c("Assay Data" = "assay", "Physical Data" = "physical"), selected = "physical"),
                                helpText("Recovery analysis compares week 6 (clean water) to week 5 (last exposure)."),
                                actionButton("run_recovery", "Run Recovery Analysis")
                              ), br(),
                              h4("Translocation Analysis (mf_counts only)"),
                              wellPanel(
                                helpText("Translocation analysis shows microfiber distribution across tissue types."),
                                selectInput("translocation_tissue_compare", "Compare Tissues:",
                                           choices = c("All tissues" = "all", "Tissue vs Gland" = "tissue_gland",
                                                      "Tissue vs Gills" = "tissue_gills", "Gland vs Gills" = "gland_gills"), selected = "all"),
                                actionButton("run_translocation", "Run Translocation Analysis")
                              )),
                  mainPanel(width = 8,
                           tabsetPanel(
                             tabPanel("Recovery Results",
                                     h5("Recovery Model Summary"), verbatimTextOutput("recovery_model_summary"),
                                     h5("Recovery Comparisons"), DT::dataTableOutput("recovery_comparisons"),
                                     plotOutput("recovery_plot", height = "400px")),
                             tabPanel("Translocation Results",
                                     h5("Microfiber Translocation Summary"), DT::dataTableOutput("translocation_summary"),
                                     plotOutput("translocation_plot", height = "400px"), br(),
                                     h5("Tissue Comparison Plot"), plotOutput("tissue_comparison_plot", height = "400px"))))))
      )
    )
  )
)

# ============================================================================
# SERVER FUNCTION - FIXED CACHING ISSUE
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
    req(input$endpoint)
    if ("mf_counts" %in% input$endpoint) {
      selectizeInput("tissue_type_phys", "Tissue Type", 
                    choices = unique(na.omit(physical_master$tissue_type)), multiple = TRUE)
    }
  })
  
  # ADD THIS NEW BLOCK HERE - STEP 2: Dynamic sample type choices
  observe({
    # Get available sample types from actual data
    if (exists("final_data") && nrow(final_data) > 0) {
      actual_samples <- sort(unique(tolower(final_data$sample_type)))
      
      # Get config samples (convert to lowercase for consistency)  
      config_samples <- if(exists("config") && !is.null(config$samples)) {
        tolower(config$samples)
      } else {
        c("tissue", "gland", "gills", "hemolymph")
      }
      
      # Use intersection of config and actual data, or just actual data
      available_samples <- if(any(config_samples %in% actual_samples)) {
        intersect(config_samples, actual_samples)
      } else {
        actual_samples
      }
      
      # Update the choices
      updateSelectizeInput(session, "sample_type", 
                           choices = available_samples)
    }
  })
  
  
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
  # # FIXED: Assay data filter - Use reactive with proper caching
  # assay_filtered_cached <- reactive({
  #   df <- final_data  # Pre-processed data
  #   
  #   if (length(input$assay)) df <- df %>% filter(assay_type %in% input$assay)
  #   if (length(input$fiber)) df <- df %>% filter(fiber_type %in% input$fiber)
  #   if (length(input$sample_type)) df <- df %>% filter(sample_type %in% input$sample_type)
  #   if (length(input$week)) df <- df %>% filter(week %in% as.integer(input$week))
  #   if (length(input$fiber_concentration)) df <- df %>% filter(fiber_concentration %in% input$fiber_concentration)
  #   if (length(input$treatment)) {
  #     if ("Control" %in% input$treatment) {
  #       df <- df %>% filter(fiber_concentration == "0" | treatment %in% setdiff(input$treatment, "Control"))
  #     } else {
  #       df <- df %>% filter(treatment %in% input$treatment)
  #     }
  #   }
  #   if (input$mode != "inter" && length(input$plate_replicate)) {
  #     df <- df %>% filter(plate_replicate %in% input$plate_replicate)
  #   }
  #   df
  # }) %>% bindCache(input$assay, input$fiber, input$sample_type, input$week, 
  #                  input$fiber_concentration, input$treatment, input$plate_replicate, input$mode)
  # 
  # # Now use bindEvent
  # assay_filtered <- assay_filtered_cached %>% bindEvent(input$run_analysis)
  
  # FIXED: Physical data filter - Use reactive with proper caching  
  phys_filtered_cached <- reactive({
    df <- physical_master  # Pre-processed data
    
    if (length(input$fiber_type_phys)) df <- df %>% filter(fiber_type %in% input$fiber_type_phys)
    if (length(input$week_phys)) df <- df %>% filter(week %in% input$week_phys)
    if (length(input$endpoint)) df <- df %>% filter(endpoint %in% input$endpoint)
    if ("mf_counts" %in% input$endpoint && length(input$tissue_type_phys)) {
      df <- df %>% filter(tissue_type %in% input$tissue_type_phys)
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
    
    # Apply blank correction for mf_counts
    df <- df %>%
      group_by(fiber_type, week, tissue_type, endpoint) %>%
      mutate(
        avg_control = if_else(endpoint == "mf_counts", mean(value[treatment == "Control"], na.rm = TRUE), 0),
        value_corr = if_else(endpoint == "mf_counts" & treatment != "Control", value - avg_control, value),
        value_corr = if_else(endpoint == "mf_counts", round(value_corr), value_corr),
        value = if_else(endpoint == "mf_counts", round(value), value)
      ) %>% ungroup()
    
    if (input$mode_phys == "baseline") {
      df
    } else {
      df %>%
        mutate(use_val = if_else(endpoint == "mf_counts", value_corr, value)) %>%
        group_by(fiber_group, week, endpoint, tissue_type, fiber_concentration) %>%
        summarise(
          mean_value = mean(use_val, na.rm = TRUE), sd_value = sd(use_val, na.rm = TRUE),
          n = n(), cv = ifelse(mean_value != 0, 100 * sd_value / mean_value, NA_real_),
          .groups = "drop"
        ) %>%
        mutate(
          mean_value = if_else(abs(mean_value) < 1, round(mean_value, 4), round(mean_value, 2)),
          sd_value = signif(sd_value, 5), cv = round(cv, 2)
        )
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
  
  # FIXED: Model-ready data for regression with proper caching
  regression_data_cached <- reactive({
    req(input$regression_dataset, input$regression_endpoint)
    
    if (input$regression_dataset == "assay") {
      df <- final_data %>%
        create_enhanced_treatment_categories() %>%
        mutate(outcome = calculated_concentration) %>%
        filter(assay_type == input$regression_endpoint) %>%
        average_assay_replicates(outcome)
    } else {
      df <- physical_master %>%
        create_enhanced_treatment_categories() %>%
        mutate(outcome = value) %>%
        filter(endpoint == input$regression_endpoint)
    }
    
    wks <- input$weeks_include %||% c("1","3","5")
    df <- df %>% filter(as.character(week) %in% wks)
    
    if ("6" %in% wks) {
      df <- df %>% mutate(is_recovery = ifelse(week == 6, 1, 0))
    } else {
      df <- df %>% mutate(is_recovery = 0)
    }
    
    df %>% normalize_controls_and_dose() %>% droplevels()
  }) %>% bindCache(input$regression_dataset, input$regression_endpoint, input$weeks_include)
  
  regression_data <- regression_data_cached
  
  # FIXED: Mixed effects model data with proper caching
  lmer_data_cached <- reactive({
    req(input$regression_dataset, input$regression_endpoint)
    
    if (input$regression_dataset == "assay") {
      df <- final_data %>%
        create_experiment_batch() %>%
        create_enhanced_treatment_categories() %>%
        mutate(
          outcome = calculated_concentration,
          tank = as.factor(tank),
          week = as.character(week),
          sample_type = as.factor(sample_type)
        ) %>%
        filter(assay_type == input$regression_endpoint) %>%
        average_assay_replicates(outcome) %>%
        mutate(grouping_var = tank)
    } else {
      df <- physical_master %>%
        create_experiment_batch() %>%
        create_enhanced_treatment_categories() %>%
        filter(endpoint == input$regression_endpoint) %>%
        mutate(
          outcome = value,
          week = as.character(week),
          tissue_type = as.factor(tissue_type)
        )

      if (identical(input$regression_endpoint, "mf_counts")) {
        df <- df %>% mutate(
          individual_id = stringr::str_extract(sample, "^[0-9]+"),
          tissue_replicate = sample
        )
      } else {
        df <- df %>% mutate(
          individual_id = as.character(sample)
        )
      }
      df <- df %>% mutate(grouping_var = as.factor(individual_id))
    }

    wks <- input$lmer_weeks_include %||% c("1","3","5")
    df <- df %>%
      filter(as.character(week) %in% wks) %>%
      mutate(
        is_recovery = ifelse(as.character(week) == "6", 1, 0),
        week = factor(as.character(week), levels = sort(unique(as.character(week))))
      )

    normalize_controls_and_dose(df) %>% droplevels()
  }) %>% bindCache(input$regression_dataset, input$regression_endpoint, input$lmer_weeks_include)
  
  lmer_data <- lmer_data_cached
  
  # Regression model fitting
  regression_model <- eventReactive(input$run_regression, {
    tryCatch({
      df <- regression_data()
      req(nrow(df) > 0)
      
      # Build formula
      form <- outcome ~ fiber_type * chem_treatment + week + 
        fiber_type:week + chem_treatment:week + is_control * fiber_type + dose_log10
      
      if (isTRUE(input$include_three_way)) form <- update(form, . ~ . + fiber_type:chem_treatment:week)
      if (isTRUE(input$dose_by_fiber)) form <- update(form, . ~ . + dose_log10:fiber_type)
      if (isTRUE(input$dose_by_treat)) form <- update(form, . ~ . + dose_log10:chem_treatment)
      
      wks_selected <- input$weeks_include %||% c("1","3","5")
      if (isTRUE(input$include_recovery_interaction) && "6" %in% wks_selected) {
        form <- update(form, . ~ . + is_recovery + is_recovery:chem_treatment + is_recovery:dose_log10)
      }
      
      lm(form, data = df)
    }, error = function(e) { 
      list(error = e$message) 
    })
  })
  
  # Mixed effects model fitting  
  lmer_model <- eventReactive(input$run_lmer, {
    tryCatch({
      df <- lmer_data()
      req(nrow(df) > 0)
      
      message("Building lmer model with ", nrow(df), " observations")
      
      # ENHANCED: More stable model building approach
      
      # Start with basic formula
      fixed_formula <- "outcome ~ fiber_type + chem_treatment + week + dose_log10 + is_control"
      
      # Add interactions carefully to avoid over-parameterization
      if (input$lmer_include_three_way && length(unique(df$week)) > 2) {
        fixed_formula <- paste(fixed_formula, "+ fiber_type:chem_treatment + fiber_type:week + chem_treatment:week")
        
        # Only add 3-way if we have enough data
        if (nrow(df) > 100) {
          fixed_formula <- paste(fixed_formula, "+ fiber_type:chem_treatment:week")
        }
      } else {
        # Add basic 2-way interactions
        fixed_formula <- paste(fixed_formula, "+ fiber_type:chem_treatment")
      }
      
      if (input$lmer_dose_by_fiber) {
        fixed_formula <- paste(fixed_formula, "+ dose_log10:fiber_type")
      }
      
      if (input$lmer_dose_by_treat) {
        fixed_formula <- paste(fixed_formula, "+ dose_log10:chem_treatment")
      }
      
      # Handle recovery interaction
      if (input$lmer_include_recovery_interaction && "6" %in% input$lmer_weeks_include) {
        fixed_formula <- paste(fixed_formula, "+ is_recovery + is_recovery:chem_treatment")
      }
      
      # SMART: Build random effects based on data availability
      random_terms <- input$lmer_random_terms
      random_formula <- ""
      
      if ("tank" %in% random_terms && "tank" %in% names(df)) {
        # Check if we have enough tanks
        n_tanks <- length(unique(df$tank))
        if (n_tanks > 3) {
          random_formula <- paste(random_formula, "+ (1|tank)")
        } else {
          message("Too few tanks (", n_tanks, ") for random effect, skipping tank")
        }
      }
      
      if ("experiment_batch" %in% random_terms && "experiment_batch" %in% names(df)) {
        # Check if we have enough batches  
        n_batches <- length(unique(df$experiment_batch))
        if (n_batches > 1) {
          random_formula <- paste(random_formula, "+ (1|experiment_batch)")
        } else {
          message("Only one experiment batch, skipping batch random effect")
        }
      }
      
      # Fallback random effect
      if (random_formula == "" && "grouping_var" %in% names(df)) {
        n_groups <- length(unique(df$grouping_var))
        if (n_groups > 3) {
          random_formula <- "+ (1|grouping_var)"
          message("Using grouping_var as fallback random effect")
        }
      }
      
      # Construct final formula
      if (random_formula != "") {
        full_formula <- paste(fixed_formula, random_formula)
      } else {
        # If no random effects possible, fall back to simple lm
        message("No suitable random effects found, this may cause issues")
        full_formula <- fixed_formula
      }
      
      formula_obj <- as.formula(full_formula)
      message("Final formula: ", deparse(formula_obj))
      
      # Fit model with enhanced control parameters for stability
      model <- lmer(formula_obj, data = df, 
                    control = lmerControl(
                      optimizer = "bobyqa",      # More stable optimizer
                      optCtrl = list(maxfun = 20000),  # More iterations
                      check.conv.singular = "warning",  # Don't fail on singular fits
                      calc.derivs = FALSE        # Skip expensive derivative calculations
                    ))
      
      # Check for convergence warnings
      if (length(model@optinfo$conv$lme4$messages) > 0) {
        message("Model convergence warnings: ", paste(model@optinfo$conv$lme4$messages, collapse = "; "))
      }
      
      model
      
    }, error = function(e) {
      message("LMER model fitting error: ", e$message)
      list(error = e$message)
    })
  })
  
  
  # ============================================================================
  # ENHANCED EMMEANS WITH COLORED SIGNIFICANCE & DID ANALYSIS
  # ============================================================================
  
  # EMMEANS calculation for regular regression 
  emmeans_results <- eventReactive(input$run_emmeans, {
    tryCatch({
      # Require that the regression model exists first
      model <- regression_model()
      validate(need(!is.null(model) && !inherits(model, "list"), "Run the regression model first"))
      
      # Calculate dose value based on selection
      dose_val <- switch(input$emmeans_dose,
                         "0" = 0,           # Control level (log10 scale)
                         "2" = log10(100),  # Low dose
                         "3" = log10(1000), # Medium dose  
                         "4" = log10(10000)) # High dose
      
      is_ctrl <- ifelse(dose_val == 0, 1, 0)
      at_list <- list(dose_log10 = dose_val, is_control = is_ctrl)
      
      # Get the grouping specification
      emmeans_spec <- switch(input$emmeans_by,
                             "treatment" = ~ chem_treatment | fiber_type,
                             "fiber_type" = ~ fiber_type,
                             "week" = ~ week,
                             "treatment_week" = ~ chem_treatment | week + fiber_type,
                             "treatment_fiber" = ~ chem_treatment | fiber_type,
                             "fiber_week" = ~ fiber_type | week,
                             "all_interactions" = ~ chem_treatment | fiber_type + week)
      
      # Get the emmeans
      emm <- emmeans(model, specs = emmeans_spec, at = at_list)
      
      list(
        emmeans = emm,
        emmeans_df = as.data.frame(emm)
      )
    }, error = function(e) {
      message("Emmeans error: ", e$message)
      list(error = e$message)
    })
  })
  
  # ENHANCED PAIRWISE COMPARISONS WITH DID SUPPORT
  output$pairwise_table <- DT::renderDataTable({
    m <- regression_model()
    validate(need(!is.null(m), "Run the regression model first."))
    
    adj_method <- input$adjustment_method %||% "tukey"
    
    if (input$comparison_type == "did") {
      # DIFFERENCE OF DIFFERENCES ANALYSIS
      
      # Get required inputs
      baseline_week <- input$did_baseline_week %||% "1"
      target_week <- input$did_target_week %||% "3"  
      scope <- input$did_scope %||% "by_fiber"
      dose_val <- as.numeric(input$emmeans_dose %||% "0")
      is_ctrl <- ifelse(dose_val == 0, 1, 0)
      
      # Calculate DiD contrasts
      tryCatch({
        # Step 1: Get EM means by treatment within fiber and week
        emm_tw <- emmeans(m, 
                          specs = ~ chem_treatment | fiber_type + week,
                          at = list(dose_log10 = dose_val, is_control = is_ctrl))
        
        # Step 2: Calculate treatment effects within each fiber × week
        trt_effects <- contrast(emm_tw, method = "trt.vs.ctrl", ref = "untreated")
        
        # Step 3: Filter to baseline and target weeks
        df_te <- as.data.frame(trt_effects) %>%
          filter(week %in% c(baseline_week, target_week))
        
        # Step 4: Calculate DiD (target - baseline) 
        if (scope == "by_fiber") {
          did_results <- df_te %>%
            arrange(fiber_type, week) %>%
            group_by(fiber_type) %>%
            summarise(
              contrast = paste0("DiD (", fiber_type[1], "): ", target_week, " - ", baseline_week),
              estimate = diff(estimate),  # target - baseline
              SE = sqrt(sum(SE^2)),      # Conservative SE approximation
              df = mean(df),
              t.ratio = estimate / SE,
              p.value = 2 * pt(-abs(t.ratio), df = df),
              .groups = "drop"
            )
        } else {
          # Pooled across fibers
          pooled <- df_te %>%
            group_by(week) %>%
            summarise(
              emmean = mean(estimate),
              SE = sqrt(mean(SE^2)),
              df = mean(df),
              .groups = "drop"
            ) %>%
            arrange(match(week, c(baseline_week, target_week)))
          
          did_results <- tibble(
            contrast = paste0("DiD (pooled): ", target_week, " - ", baseline_week),
            estimate = pooled$emmean[2] - pooled$emmean[1],
            SE = sqrt(sum(pooled$SE^2)),
            df = mean(pooled$df),
            t.ratio = estimate / SE,
            p.value = 2 * pt(-abs(t.ratio), df = df)
          )
        }
        
        # Apply adjustment
        if (adj_method != "none" && nrow(did_results) > 1) {
          did_results$p.value <- p.adjust(did_results$p.value, method = adj_method)
        }
        
        dfp <- did_results
        
      }, error = function(e) {
        dfp <- data.frame(
          Message = paste0("DiD Error: ", e$message),
          stringsAsFactors = FALSE
        )
      })
      
    } else {
      # REGULAR PAIRWISE OR TREATMENT VS CONTROL
      emm_res <- emmeans_results()
      validate(need(!is.null(emm_res$emmeans), "Calculate EM Means first"))
      
      emm <- emm_res$emmeans
      
      if (input$comparison_type == "control") {
        # Treatment vs Control
        primary <- names(emm@levels)[1] 
        validate(need(identical(primary, "chem_treatment"), 
                      "Treatment vs Control requires EM means by Treatment Groups."))
        
        lvls <- emm@levels[[primary]]
        validate(need("untreated" %in% lvls, "No 'untreated' level present in EM grid."))
        
        contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = "untreated", adjust = adj_method)
      } else {
        # All pairwise
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
    
    # Create colored DataTable
    DT::datatable(
      dfp,
      options = list(pageLength = 20, scrollX = TRUE, dom = 'Blfrtip'),
      filter = "top",
      rownames = FALSE
    ) %>%
      DT::formatStyle(
        "Significant",
        backgroundColor = DT::styleEqual(
          c("Yes (***)", "Yes (**)", "Yes (*)", "No", "Unable to compute"),
          c("#28a745", "#28a745", "#28a745", "#ffffff", "#ffe6cc")  # Green for significant, white for non-significant
        ),
        color = DT::styleEqual(
          c("Yes (***)", "Yes (**)", "Yes (*)", "No", "Unable to compute"),
          c("#ffffff", "#ffffff", "#ffffff", "#000000", "#000000")   # White text on green, black on white
        )
      ) %>%
      DT::formatRound(columns = c("estimate", "SE", "t.ratio"), digits = 4) %>%
      DT::formatRound(columns = "p.value", digits = 6)
  })
  
  # UI CONTROLS FOR DIFFERENCE OF DIFFERENCES
  output$did_controls <- renderUI({
    df <- regression_data()
    req(nrow(df) > 0)
    
    wks <- sort(unique(as.character(df$week)))
    validate(need(length(wks) >= 2, "At least two weeks are required for DiD."))
    
    tagList(
      selectInput("did_baseline_week", "Baseline week:", 
                  choices = wks, selected = wks[1]),
      selectInput("did_target_week", "Target week:", 
                  choices = wks, selected = wks[min(2, length(wks))]),
      radioButtons("did_scope", "Scope:",
                   choices = c("By fiber type" = "by_fiber", "Pooled across fibers" = "pooled"),
                   selected = "by_fiber", inline = TRUE),
      helpText("DiD compares the treatment effect between two time points. 
             Positive values indicate treatment effect increased over time.")
    )
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
      # Enhanced plot with better styling
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
      # Fallback to basic ggplot
      df <- results$emmeans_df
      ggplot(df, aes(x = emmean, y = rownames(df))) +
        geom_point(size = 3, color = "#2c3e50") +
        geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2) +
        theme_minimal() +
        labs(
          title = "Estimated Marginal Means", 
          x = "Estimated Mean", 
          y = "Groups"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
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
      # ENHANCED: Better plot with error handling for negative variances
      
      # Check if emmeans has issues with confidence intervals
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
        # Use valid rows only
        df_clean <- df_plot[valid_rows, ]
        
        # Create the plot manually to handle negative variances better
        ggplot(df_clean, aes(x = emmean, y = factor(rownames(df_clean)))) +
          geom_point(size = 3, color = "#2c3e50") +
          geom_errorbarh(aes(xmin = pmax(lower.CL, emmean - 3*SE), 
                             xmax = pmin(upper.CL, emmean + 3*SE)), 
                         height = 0.2, color = "#3498db") +
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
    }, error = function(e) {
      # Ultimate fallback plot
      message("Plot error: ", e$message)
      ggplot(data.frame(x = 1, y = 1, message = "Plot failed due to model convergence issues")) +
        geom_text(aes(x, y, label = message), size = 5, color = "red") +
        theme_void() +
        labs(title = "Mixed Effects Plot Error", subtitle = "Try simpler model settings")
    })
  })
  
  output$lmer_did_controls <- renderUI({
    df <- lmer_data()
    req(nrow(df) > 0)
    
    wks <- sort(unique(as.character(df$week)))
    validate(need(length(wks) >= 2, "At least two weeks are required for DiD."))
    
    tagList(
      selectInput("lmer_did_baseline_week", "Baseline week:", 
                  choices = wks, selected = wks[1]),
      selectInput("lmer_did_target_week", "Target week:", 
                  choices = wks, selected = wks[min(2, length(wks))]),
      radioButtons("lmer_did_scope", "Scope:",
                   choices = c("By fiber type" = "by_fiber", "Pooled across fibers" = "pooled"),
                   selected = "by_fiber", inline = TRUE),
      helpText("DiD compares the treatment effect between two time points for mixed effects models.")
    )
  })
  
  # ============================================================================
  # REGRESSION OUTPUTS
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
  # MIXED EFFECTS OUTPUTS
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
  # RECOVERY & TRANSLOCATION OUTPUTS
  # ============================================================================
  
  output$recovery_model_summary <- renderPrint({
    analysis <- recovery_analysis()
    validate(need(!is.null(analysis$model), "Click 'Run Recovery Analysis' to fit the model."))
    summary(analysis$model)
  })
  
  output$recovery_comparisons <- DT::renderDataTable({
    analysis <- recovery_analysis()
    validate(need(!is.null(analysis$model), ""))
    
    # Create recovery comparisons
    emm <- emmeans(analysis$model, ~ treatment_category | is_recovery)
    contrasts <- contrast(emm, method = "pairwise", by = "is_recovery")
    
    DT::datatable(as.data.frame(contrasts), options = list(pageLength = 10), filter = "top", rownames = FALSE)
  })
  
  output$recovery_plot <- renderPlot({
    analysis <- recovery_analysis()
    validate(need(!is.null(analysis$data), ""))
    
    ggplot(analysis$data, aes(x = treatment_category, y = outcome, fill = is_recovery)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      theme_minimal() +
      labs(title = "Recovery Analysis: Week 6 vs Week 5",
           x = "Treatment Category", y = "Outcome Value", fill = "Period") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$translocation_summary <- DT::renderDataTable({
    analysis <- translocation_analysis()
    validate(need(!is.null(analysis), "Click 'Run Translocation Analysis' to generate summary."))
    
    DT::datatable(analysis, options = list(pageLength = 15), filter = "top", rownames = FALSE)
  })
  
  output$translocation_plot <- renderPlot({
    analysis <- translocation_analysis()
    validate(need(!is.null(analysis), ""))
    
    ggplot(analysis, aes(x = week, y = mean_count, fill = tissue_type)) +
      geom_col(position = position_dodge(width = 0.8), alpha = 0.7) +
      facet_grid(fiber_type ~ treatment_category) +
      theme_minimal() +
      labs(title = "Microfiber Translocation Across Tissues",
           x = "Week", y = "Mean Count", fill = "Tissue Type") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$tissue_comparison_plot <- renderPlot({
    analysis <- translocation_analysis()
    validate(need(!is.null(analysis), ""))
    
    ggplot(analysis, aes(x = tissue_type, y = prop_nonzero, fill = treatment_category)) +
      geom_col(position = position_dodge(width = 0.8), alpha = 0.7) +
      facet_wrap(~ week) +
      theme_minimal() +
      labs(title = "Proportion of Non-zero Counts by Tissue",
           x = "Tissue Type", y = "Proportion Non-zero", fill = "Treatment") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
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