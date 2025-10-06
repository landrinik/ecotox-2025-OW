library(shiny)
library(tidyverse)
library(readxl)
library(janitor)
library(here)
library(DT)
library(glue)
library(yaml)
library(stringr)
library(shinyjs)
library(emmeans)
library(broom)        
library(broom.mixed)
library(car)
library(lme4)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- Load config and data ---

config <- yaml::read_yaml(here("config.yml"))
filter_choices <- function(x) unlist(x[!startsWith(as.character(x), "#")])
tissue_weights <- read_excel(here("tissue_weights_clean.xlsx")) %>% clean_names()

# -- Assay Data --

final_data_raw <- read_excel(here("Assay_endpoint_DATA_cleaned.xlsx")) %>% clean_names()


create_experiment_batch <- function(df) {
  df %>%
    mutate(
      experiment_batch = case_when(
        tolower(fiber_type) == "cotton" ~ "April_Cotton",
        tolower(fiber_type) == "pet" ~ "January_PET",
        TRUE ~ NA_character_
      ),
      experiment_batch = as.factor(experiment_batch)
    )
}

# Annotate assay-specific columns with robust parsing
annotate_common <- function(df) {
  df %>%
    mutate(
      week_start = as.numeric(str_extract(sample_week, "^[0-9]+")),
      week_end = as.numeric(str_extract(sample_week, "[0-9]+$")),
      well_number = as.numeric(str_extract(well_name, "\\d+")),
      well_letter = str_extract(well_name, "[A-H]"),
      week = case_when(
        well_number == 6 & well_letter %in% c("A", "B") ~ week_start,
        well_number == 6 & well_letter %in% c("C", "D", "E", "F", "G", "H") ~ week_end,
        well_number < 6 ~ week_start,
        well_number > 6 ~ week_end,
        TRUE ~ NA_real_
      ),
      tank = case_when(
        well_number <= 5 & well_letter %in% c("A", "B") ~ (well_number - 1) * 4 + 1,
        well_number <= 5 & well_letter %in% c("C", "D") ~ (well_number - 1) * 4 + 2,
        well_number <= 5 & well_letter %in% c("E", "F") ~ (well_number - 1) * 4 + 3,
        well_number <= 5 & well_letter %in% c("G", "H") ~ (well_number - 1) * 4 + 4,
        well_number == 6 & well_letter %in% c("A", "B") ~ 21,
        well_number == 6 & well_letter %in% c("C", "D") ~ 1,
        well_number == 6 & well_letter %in% c("E", "F") ~ 2,
        well_number == 6 & well_letter %in% c("G", "H") ~ 3,
        well_number > 6 & well_letter %in% c("A", "B") ~ (well_number -7) * 4 + 4,
        well_number > 6 & well_letter %in% c("C", "D") ~ (well_number -7) * 4 + 5,
        well_number > 6 & well_letter %in% c("E", "F") ~ (well_number -7) * 4 + 6,
        well_number > 6 & well_letter %in% c("G", "H") ~ (well_number -7) * 4 + 7,
        TRUE ~ NA_real_
      ),
      fiber_concentration = case_when(
        tank %in% c(1, 2, 3) ~ "0",
        tank %in% c(4, 5, 6, 13, 14, 15) ~ "100",
        tank %in% c(7, 8, 9, 16, 17, 18) ~ "1000",
        tank %in% c(10, 11, 12, 19, 20, 21) ~ "10000",
        TRUE ~ NA_character_
      ),
      treatment = case_when(
        tank %in% c(1, 2, 3) ~ "Control",
        tank <= 12 & !(tank %in% c(1, 2, 3)) ~ "Untreated",
        tank > 12 ~ "Treated",
        TRUE ~ NA_character_
      ),
      fiber_group = if_else(treatment == "Control" & fiber_concentration == "0", "Control", paste(fiber_type, treatment))
    ) %>%
    # suppress unused columns if desired
    dplyr::select(-week_start, -week_end, -well_number)
}

final_data <- annotate_common(final_data_raw) %>%
  mutate(calculated_concentration = suppressWarnings(as.numeric(calculated_concentration)))

week_choices_assay <- sort(unique(na.omit(final_data$week)))


# -- Physical Data --

physical_master <- read_excel(here("physical_endpoint_master_sheet.xlsx")) %>% clean_names() %>%
  mutate(
    tank = as.integer(str_extract(sample, "^[0-9]+")),
    fiber_concentration = case_when(
      tank %in% c(1, 2, 3) ~ "0",
      tank %in% c(4, 5, 6, 13, 14, 15) ~ "100",
      tank %in% c(7, 8, 9, 16, 17, 18) ~ "1000",
      tank %in% c(10, 11, 12, 19, 20, 21) ~ "10000",
      TRUE ~ NA_character_
    ),
    treatment = case_when(
      tank %in% c(1, 2, 3) ~ "Control",
      tank <= 12 & !(tank %in% c(1, 2, 3)) ~ "Untreated",
      tank > 12 ~ "Treated",
      TRUE ~ NA_character_
    ),
    fiber_group = if_else(treatment == "Control" & fiber_concentration == "0","Control", paste(fiber_type, treatment))
  )

summarize_baseline <- function(df, threshold) {
  df %>%
    group_by(fiber_group, assay_type, week, sample_type, tank, well_letter,
             fiber_concentration, well_name, plate_replicate) %>%
    summarise(
      vals = list(calculated_concentration[is.finite(calculated_concentration)]),
      mean_activity = ifelse(length(vals[[1]]) == 0, NA_real_, mean(vals[[1]], na.rm = TRUE)),
      sd_activity = ifelse(length(vals[[1]]) == 0, NA_real_, sd(vals[[1]], na.rm = TRUE)),
      has_inf = any(is.infinite(calculated_concentration)),
      has_zero = any(calculated_concentration == 0),
      has_na = any(is.na(calculated_concentration)),
      cv_percent = round(if_else(mean_activity != 0 & !is.na(mean_activity),
                                 (sd_activity / mean_activity) * 100, NA_real_), 2),
      cv_flag = cv_percent > threshold,
      .groups = "drop"
    ) %>%
    mutate(
      # Create display version of mean_activity with asterisk if any flag true
      mean_activity_display = if_else(
        has_inf | has_zero | has_na,
        paste0(round(mean_activity, 2), "*"),
        as.character(round(mean_activity, 2))
      )
    ) %>%
    # Remove columns vars you don't want to display after summarise
    dplyr::select(-has_inf, -has_zero, -has_na, -vals)
}

summarize_wrangle <- function(df, threshold) {
  df %>%
    left_join(tissue_weights,
              by = c("assay_type", "fiber_type", "sample_type", "sample_week", "well_name")) %>%
    group_by(fiber_group, assay_type, week, sample_type, tank, well_letter,
             fiber_concentration, well_name) %>%
    summarise(
      vals = list(calculated_concentration[is.finite(calculated_concentration) & calculated_concentration != 0]),
      mean_activity = ifelse(length(vals[[1]]) == 0, NA_real_, mean(unlist(vals), na.rm = TRUE)),
      sd_activity_raw = ifelse(length(vals[[1]]) == 0, NA_real_, sd(unlist(vals), na.rm = TRUE)),
      has_inf = any(is.infinite(calculated_concentration)),
      has_zero = any(calculated_concentration == 0),
      has_na = any(is.na(calculated_concentration)),
      .groups = "drop"
    ) %>%
    mutate(
      sd_activity = round(sd_activity_raw, 4),
      cv_activity = if_else(mean_activity != 0 & !is.na(mean_activity),
                            (sd_activity / mean_activity) * 100,
                            NA_real_),
      cv_percent = round(cv_activity, 2),
      cv_flag = cv_percent > threshold,
      mean_activity_display = if_else(has_inf | has_na,
                                      paste0(round(mean_activity, 2), "*"),
                                      as.character(round(mean_activity, 2)))
    ) %>%
    dplyr::select(-has_inf, -has_zero, -has_na, -vals) %>%
    dplyr::select(fiber_group, assay_type, week, sample_type, tank,
                  well_letter, fiber_concentration, well_name,
                  mean_activity, mean_activity_display, sd_activity, cv_percent, cv_flag)
}

summarize_inter <- function(df, threshold) {
  df %>%
    group_by(fiber_group, assay_type, week, sample_type, fiber_concentration, plate_replicate) %>%
    summarise(
      vals = list(calculated_concentration[is.finite(calculated_concentration)]),
      mean_activity = ifelse(length(vals[[1]]) == 0, NA_real_, mean(vals[[1]], na.rm = TRUE)),
      sd_activity = ifelse(length(vals[[1]]) == 0, NA_real_, sd(vals[[1]], na.rm = TRUE)),
      n_samples = n(),
      has_inf = any(is.infinite(calculated_concentration)),
      has_zero = any(calculated_concentration == 0),
      has_na = any(is.na(calculated_concentration)),
      .groups = "drop"
    ) %>%
    mutate(
      cv = if_else(mean_activity == 0, NA_real_, 100 * sd_activity / mean_activity),
      cv = round(cv, 2),
      cv_flag = cv > threshold,
      mean_activity_display = if_else(
        has_inf | has_zero | has_na,
        paste0(format(mean_activity, digits = 4, scientific = FALSE), "*"),
        format(mean_activity, digits = 4, scientific = FALSE)
      ),
      sd_activity = signif(sd_activity, 5)
    ) %>%
    dplyr::select(-has_inf, -has_zero, -has_na, -vals)
}

# --- UI ---

ui <- fluidPage(
  useShinyjs(),
  titlePanel("ECOTOX Assay Dashboard"),
  actionButton("toggleSidebar", "Toggle Sidebar"),
  
  fluidRow(
    column(
      width = 3,
      id = "sidebar",
      h4("Data Mode:"),
      radioButtons("active_dataset", NULL, choices = c("Assay Data" = "assay", "Physical Endpoint Data" = "physical"), selected = "assay"),
      
      # Assay controls
      conditionalPanel(
        condition = "input.active_dataset == 'assay'",
        h4("Processing mode:"),
        radioButtons(
          "mode", NULL,
          choices = c("Baseline (filter only)" = "baseline",
                      "Well-level CV (wrangle_v3_avg)" = "wrangle",
                      "Inter-assay CV (inter_assay_CV)" = "inter"),
          selected = "baseline"
        ),
        br(),
        selectizeInput("assay", "Assay Type:", choices = filter_choices(config$assays) %||% unique(final_data$assay_type), multiple = TRUE),
        selectizeInput("fiber", "Fiber Type:", choices = filter_choices(config$fibers) %||% unique(final_data$fiber_type), multiple = TRUE),
        selectizeInput("sample_type", "Sample Type:", choices = filter_choices(config$samples) %||% unique(final_data$sample_type), multiple = TRUE),
        selectizeInput("week", "Week:", choices = week_choices_assay, multiple = TRUE),
        selectizeInput("treatment", "Treatment:", choices = c("Treated", "Untreated", "Control"), multiple = TRUE),
        selectizeInput("fiber_concentration", "Fiber Concentration:", choices = c("0", "100", "1000", "10000"), multiple = TRUE),
        conditionalPanel(condition = "input.mode != 'inter'",
                         selectizeInput("plate_replicate", "Plate Replicate:", choices = unique(final_data$plate_replicate), multiple = TRUE),
                         sliderInput("cv_threshold", "CV Threshold (%)", min = 0, max = 100, value = 15)
        ),
        selectizeInput("x_axis_assay", "X-axis:", choices = c(
          "Fiber Group" = "fiber_group",
          "Week" = "week",
          "Sample Type" = "sample_type",
          "Treatment" = "treatment",
          "Fiber Concentration" = "fiber_concentration"
        ), selected = "fiber_group"),
        actionButton("run_analysis", "Run Analysis"),
        br(),
        downloadButton("download_data", "Download Results")
      ),
      
      # Physical controls
      conditionalPanel(
        condition = "input.active_dataset == 'physical'",
        radioButtons("mode_phys", "Mode:", choices = c("Baseline" = "baseline", "Inter CV" = "inter"), selected = "baseline"),
        selectizeInput("fiber_type_phys", "Fiber Type", choices = unique(physical_master$fiber_type), multiple = TRUE),
        selectizeInput("week_phys", "Week", choices = sort(unique(physical_master$week)), multiple = TRUE),
        selectizeInput("endpoint", "Endpoint", choices = unique(physical_master$endpoint), multiple = TRUE),
        uiOutput("tissue_type_ui"),
        selectizeInput("fiber_concentration_phys", "Fiber Concentration", choices = c("0","100","1000","10000"), multiple = TRUE),
        selectizeInput("treatment_phys", "Treatment", choices = c("Control","Untreated","Treated"), multiple = TRUE),
        selectizeInput("x_axis_phys", "X Axis", choices = c(
          "Fiber Group"="fiber_group","Week"="week","Fiber Type"="fiber_type","Endpoint"="endpoint","Tissue Type"="tissue_type"
        ), selected = "week"),
        actionButton("run_phys_analysis", "Run Physical Analysis"),
        downloadButton("download_phys_data", "Download Physical Data")
      )
    ),
    
    column(
      9,
      
      tabsetPanel(
        id = "analysis_tabs",
        
        tabPanel(
          "Visualization",
          h4("Results"),
          DTOutput("table"),
          br(),
          tags$div(
            style = "width: 100%; max-width: 1200px; aspect-ratio: 4/3;",
            plotOutput("dist_plot", width = "100%", height = "100%")
          ),
          br(),
          tags$div(
            style = "width: 100%; max-width: 1200px; aspect-ratio: 4/3;",
            plotOutput("dot_plot", width = "100%", height = "100%")
          ),
          br(),
          verbatimTextOutput("info")
        ),
        
        tabPanel(
          "Regression Analysis",
          h4("Multiple Linear Regression"),
          fluidRow(
            column(6,
                   h5("Model Settings"),
                   selectInput("regression_endpoint", "Select Endpoint:",
                               choices = NULL),
                   selectInput("regression_dataset", "Dataset:",
                               choices = c("Assay Data" = "assay", 
                                           "Physical Data" = "physical"),
                               selected = "physical"),
                   checkboxInput("include_interaction", 
                                 "Include Treatment × Week Interaction", 
                                 value = TRUE),
                   checkboxInput("include_concentration", 
                                 "Include Fiber Concentration", 
                                 value = TRUE),
                   actionButton("run_regression", "Run Regression Model", 
                                class = "btn-primary")
            ),
            column(6,
                   h5("Post-hoc Comparisons"),
                   selectInput("emmeans_by", "Calculate EM Means by:",
                               choices = c("Treatment Category" = "treatment_category",
                                           "Week" = "week",
                                           "Treatment × Week" = "both")),
                   selectInput("comparison_type", "Comparison Type:",
                               choices = c("All Pairwise" = "pairwise",
                                           "vs. Control Only" = "control")),
                   selectInput("adjustment_method", "P-value Adjustment:",
                               choices = c("Tukey HSD" = "tukey",
                                           "Dunnett (vs Control)" = "dunnett",
                                           "Bonferroni" = "bonferroni",
                                           "Holm" = "holm",
                                           "FDR/Benjamini-Hochberg" = "fdr",
                                           "Scheffe" = "scheffe",
                                           "None (use with caution)" = "none"),
                               selected = "tukey"),
                   actionButton("run_emmeans", "Calculate EM Means")
            )
            
          ),
          hr(),
          h5("Model Summary"),
          verbatimTextOutput("model_summary"),
          h5("ANOVA Table"),
          DTOutput("anova_table"),
          h5("Diagnostic Plots"),
          plotOutput("diagnostic_plots", height = "600px"),
          h5("Estimated Marginal Means"),
          plotOutput("emmeans_plot"),
          DTOutput("emmeans_table"),
          h5("Pairwise Comparisons"),
          DTOutput("pairwise_table")
        ),
        
        # NEW TAB - Mixed Effects Analysis
        tabPanel(
          "Mixed Effects Analysis",
          h4("Linear Mixed Effects Regression (Hierarchical Modeling)"),
          
          fluidRow(
            column(6,
                   h5("Model Settings"),
                   selectInput("lmer_endpoint", "Select Endpoint:",
                               choices = NULL),
                   selectInput("lmer_dataset", "Dataset:",
                               choices = c("Assay Data" = "assay", 
                                           "Physical Data" = "physical"),
                               selected = "physical"),
                   checkboxInput("lmer_include_interaction",
                                 "Include Treatment × Week Interaction",
                                 value = TRUE),
                   checkboxInput("lmer_include_concentration",
                                 "Include Fiber Concentration",
                                 value = FALSE),  # Default off to avoid nesting issues
                   selectInput("random_structure", "Random Effects Structure:",
                               choices = c(
                                 "Tank/Sample (Random Intercept): (1|tank)" = "intercept",
                                 "Tank/Sample + Experiment Batch: (1|tank) + (1|experiment_batch)" = "batch",
                                 "Nested Tissue (mf_counts only): (1|individual/tissue)" = "nested",
                                 "Full Model: Batch + Tank + Tissue" = "full"
                               ),
                               selected = "intercept"),
                   helpText("Tank for assay data. Sample for physical data. Batch accounts for Cotton (April) vs PET (January) experiments.",
                            style = "color: #666;"),
                   helpText("⚠️ NOTE: Week effects and interactions require selecting 2+ weeks in the sidebar.",
                            style = "color: #d9534f;"),
                   helpText("⚠️ NOTE: Batch effects require BOTH Cotton AND PET data. If only one fiber type is selected, batch effect will be excluded.",
                            style = "color: #d9534f;"),
                   
                   actionButton("run_lmer", "Run Mixed Model",
                                class = "btn-primary")
            ),
            column(6,
                   h5("Post-hoc Comparisons"),
                   selectInput("lmer_emmeans_by", "Calculate EM Means by:",
                               choices = c("Treatment Category" = "treatment_category",
                                           "Week" = "week",
                                           "Treatment × Week" = "both")),
                   selectInput("lmer_comparison_type", "Comparison Type:",
                               choices = c("All Pairwise" = "pairwise",
                                           "vs. Control Only" = "control")),
                   selectInput("lmer_adjustment_method", "P-value Adjustment:",  # ADD THIS
                               choices = c("Tukey HSD" = "tukey",
                                           "Dunnett (vs Control)" = "dunnett",
                                           "Holm" = "holm",
                                           "Bonferroni" = "bonferroni",
                                           "FDR/Benjamini-Hochberg" = "fdr",
                                           "None (use with caution)" = "none"),
                               selected = "tukey"),
                   actionButton("run_lmer_emmeans", "Calculate EM Means"),
                   br(), br(),
                   downloadButton("download_lmer", "Download Results")
            )
          ),
          
          hr(),
          
          h5("Model Summary"),
          verbatimTextOutput("lmer_summary"),
          
          h5("Fixed Effects ANOVA"),
          DTOutput("lmer_anova_table"),
          
          h5("Random Effects Variance"),
          verbatimTextOutput("lmer_random_effects"),
          
          h5("Model Diagnostics"),
          plotOutput("lmer_diagnostics", height = "600px"),
          
          h5("Estimated Marginal Means"),
          plotOutput("lmer_emmeans_plot"),
          DTOutput("lmer_emmeans_table"),
          
          h5("Pairwise Comparisons"),
          DTOutput("lmer_pairwise_table")
        )
      ) 
      
    )
    
  )
)


# --- Server ---

server <- function(input, output, session) {
  
  observeEvent(input$toggleSidebar, {
    toggle("sidebar")
  })
  
  # Dynamic tissue filter UI for mf_counts endpoint
  output$tissue_type_ui <- renderUI({
    req(input$endpoint)
    if ("mf_counts" %in% input$endpoint) {
      selectizeInput("tissue_type_phys", "Tissue Type", choices = unique(na.omit(physical_master$tissue_type)), multiple = TRUE)
    } else {
      NULL
    }
  })
  
  # Assay data filter pipeline
  assay_filtered <- eventReactive(input$run_analysis, {
    df <- final_data
    
    if (length(input$assay)) df <- df %>% filter(assay_type %in% input$assay)
    if (length(input$fiber)) df <- df %>% filter(fiber_type %in% input$fiber)
    if (length(input$sample_type)) df <- df %>% filter(sample_type %in% input$sample_type)
    if (length(input$week)) df <- df %>% filter(week %in% as.integer(input$week))
    
    if (length(input$fiber_concentration) && "fiber_concentration" %in% names(df)) {
      df <- df %>% filter(fiber_concentration %in% input$fiber_concentration)
    }
    
    if (length(input$treatment) && "treatment" %in% names(df)) {
      if ("Control" %in% input$treatment) {
        df <- df %>% filter(fiber_concentration == "0" | treatment %in% setdiff(input$treatment, "Control"))
      } else {
        df <- df %>% filter(treatment %in% input$treatment)
      }
    }
    
    if (input$mode != "inter" && length(input$plate_replicate) && "plate_replicate" %in% names(df)) {
      df <- df %>% filter(plate_replicate %in% input$plate_replicate)
    }
    
    df
  })
  
  # Physical data filter pipeline
  phys_filtered <- eventReactive(input$run_phys_analysis, {
    df <- physical_master
    
    if (length(input$fiber_type_phys)) df <- df %>% filter(fiber_type %in% input$fiber_type_phys)
    if (length(input$week_phys)) df <- df %>% filter(week %in% input$week_phys)
    if (length(input$endpoint)) df <- df %>% filter(endpoint %in% input$endpoint)
    if ("mf_counts" %in% input$endpoint && length(input$tissue_type_phys)) df <- df %>% filter(tissue_type %in% input$tissue_type_phys)
    if (length(input$fiber_concentration_phys)) df <- df %>% filter(fiber_concentration %in% input$fiber_concentration_phys)
    if (length(input$treatment_phys)) df <- df %>% filter(treatment %in% input$treatment_phys)
    
    df
  })
  
  # Assay summarizations based on selected mode
  summarized <- reactive({
    req(input$active_dataset == "assay")
    df <- assay_filtered()
    if (nrow(df) == 0) return(NULL)
    threshold <- input$cv_threshold
    
    # Implementations assumed available in your code base
    if (input$mode == "baseline") {
      summarize_baseline(df, threshold)
    } else if (input$mode == "wrangle") {
      summarize_wrangle(df, threshold)
    } else {
      summarize_inter(df, threshold)
    }
  })
  
  # Physical summarization baseline raw data; inter mode compute CV
  phys_summarized <- reactive({
    df <- phys_filtered()
    req(nrow(df) > 0)
    
    # Compute blank correction value within each group (fiber_type, week, tissue_type, endpoint)
    df <- df %>%
      group_by(fiber_type, week, tissue_type, endpoint) %>%
      mutate(
        avg_control = if_else(
          endpoint == "mf_counts",
          mean(value[treatment == "Control"], na.rm = TRUE),
          0
        ),
        # Apply correction only to non-controls; controls show their raw count
        value_corr = if_else(
          endpoint == "mf_counts" & treatment != "Control",
          value - avg_control,
          value
        ),
        # Always round for counts (controls and non-controls)
        value_corr = if_else(endpoint == "mf_counts", round(value_corr), value_corr),
        value = if_else(endpoint == "mf_counts", round(value), value)
      ) %>%
      ungroup()
    
    if (input$mode_phys == "baseline") {
      # Return original values for all, but display rounded counts in table
      df
    } else {
      df %>%
        mutate(use_val = if_else(endpoint == "mf_counts", value_corr, value)) %>%
        group_by(fiber_group, week, endpoint, tissue_type, fiber_concentration) %>%
        summarise(
          mean_value = mean(use_val, na.rm = TRUE),
          sd_value = sd(use_val, na.rm = TRUE),
          n = n(),
          cv = ifelse(mean_value != 0, 100 * sd_value / mean_value, NA_real_),
          .groups = "drop"
        ) %>%
        mutate(
          # Round appropriately based on magnitude
          mean_value = if_else(
            abs(mean_value) < 1,
            round(mean_value, 4),
            round(mean_value, 2)  
          ),
          sd_value = signif(sd_value, 5),
          cv = round(cv, 2)
        )
    }
  })
  
  
  # Active table for display
  active_table <- reactive({
    if (input$active_dataset == "assay") summarized() else phys_summarized()
  })
  
  # Render data table
  output$table <- renderDT({
    df <- active_table()
    if (is.null(df) || nrow(df) == 0) return(datatable(data.frame(Message = "No data for selected filters / mode")))
    
    # Round only the displayed mf_counts value
    if("endpoint" %in% names(df) && "value" %in% names(df)) {
      df <- df %>% mutate(
        value = if_else(endpoint == "mf_counts", round(value), value)
      )
    }
    
    # Hide avg_control and value_corr for non-mf_counts endpoints (only in baseline mode)
    if("endpoint" %in% names(df) && 
       !all(df$endpoint == "mf_counts") && 
       "avg_control" %in% names(df) && 
       "value_corr" %in% names(df)) {
      df <- df %>% select(-avg_control, -value_corr)
    }
    
    
    datatable(df, filter = "top", options = list(pageLength = 15))
  })
  
  
  # Distribution plot
  output$dist_plot <- renderPlot({
    df <- active_table()
    req(df, nrow(df) > 0)
    
    if (input$active_dataset == "assay") {
      x_axis_var <- input$x_axis_assay %||% "fiber_group"
      # Convert week to factor if it's the x-axis
      if (x_axis_var == "week" && "week" %in% names(df)) {
        df$week <- factor(df$week)
      }
      ggplot(df, aes_string(x = x_axis_var, y = "mean_activity", fill = "fiber_concentration")) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
        facet_wrap(~ assay_type) +
        theme_minimal() +
        labs(title = glue("Assay Activity by {str_to_title(gsub('_', ' ', x_axis_var))}"),
             x = str_to_title(gsub("_", " ", x_axis_var)), y = "Mean Activity", fill = "Fiber Concentration") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } else {
      x_axis_var <- input$x_axis_phys %||% "week"
      x_axis_factor <- paste0("factor(", x_axis_var, ")")  # Make x-axis discrete factor
      
      if (input$mode_phys == "baseline") {
        ggplot(df, aes_string(x = x_axis_factor, y = "value", fill = "fiber_concentration")) +
          geom_boxplot(outlier.shape = NA, alpha = 0.7) +
          geom_jitter(width = 0.18, alpha = 0.6, size = 2) +
          facet_wrap(~ endpoint, scales = "free_y") +
          theme_minimal() +
          labs(title = glue("Physical Baseline by {str_to_title(gsub('_', ' ', x_axis_var))}"),
               x = str_to_title(gsub("_", " ", x_axis_var)), y = "Value", fill = "Fiber Concentration") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      } else {
        ggplot(df, aes_string(x = x_axis_factor, y = "mean_value", fill = "fiber_concentration")) +
          geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.7) +  # Bar plot for summarized means
          geom_point(aes_string(group = "fiber_concentration"),
                     position = position_dodge(width = 0.8), size = 2, color = "black", show.legend = FALSE) +
          facet_wrap(~ endpoint, scales = "free_y") +
          theme_minimal() +
          labs(title = glue("Physical Inter CV by {str_to_title(gsub('_', ' ', x_axis_var))}"),
               x = str_to_title(gsub("_", " ", x_axis_var)), y = "Mean Value", fill = "Fiber Concentration") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
    }
  })
  
  # Dot plot
  output$dot_plot <- renderPlot({
    df <- assay_filtered()
    req(!is.null(df), nrow(df) > 0)
    required_cols <- c("fiber_group", "calculated_concentration", "fiber_concentration", "assay_type")
    if (!all(required_cols %in% colnames(df))) return(NULL)
    ggplot(df, aes(x = fiber_group, y = calculated_concentration, color = fiber_concentration)) +
      geom_jitter(width = 0.3, alpha = 0.7, size = 2) +
      facet_wrap(~ assay_type) +
      theme_minimal() +
      labs(
        title = "Sample-Level Distribution by Group and Fiber Concentration",
        x = "Fiber Group",
        y = "Calculated Concentration",
        color = "Fiber Concentration"
      ) +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  })
  
  
  # Info text
  output$info <- renderText({
    df <- active_table()
    if (is.null(df)) return("No data for current filters / mode.")
    if (input$active_dataset == "assay" && "cv_flag" %in% colnames(df)) {
      flagged <- sum(dplyr::coalesce(df$cv_flag, FALSE), na.rm = TRUE)
      total <- nrow(df)
      glue("Assay mode: {input$mode}. Showing {total} summary rows. {flagged} rows exceed CV threshold ({input$cv_threshold}%).")
    } else {
      glue("Physical mode: {input$mode_phys}. Showing {nrow(df)} rows. Filtered endpoint: {paste(unique(df$endpoint), collapse = ', ')}")
    }
  })
  
  # Download handler for assay data
  output$download_data <- downloadHandler(
    filename = function() { paste0("ECOTOX_", input$mode, "_", Sys.Date(), ".xlsx") },
    content = function(file) {
      df <- summarized()
      writexl::write_xlsx(if (is.null(df)) data.frame(Message = "No data for selected filters / mode") else df, file)
    }
  )
  
  # Download handler for physical data
  output$download_phys_data <- downloadHandler(
    filename = function() { paste0("ECOTOX_physical_", Sys.Date(), ".xlsx") },
    content = function(file) {
      df <- phys_summarized()
      writexl::write_xlsx(if (is.null(df)) data.frame(Message = "No data for selected filters / mode") else df, file)
    }
  )
  
  # ============================================================================
  # REGRESSION ANALYSIS MODULE - ADD THIS ENTIRE SECTION
  # ============================================================================
  
  # Dynamic endpoint selection for regression
  observe({
    if (input$regression_dataset == "assay") {
      choices <- unique(assay_filtered()$assay_type)
    } else {
      choices <- unique(phys_filtered()$endpoint)
    }
    updateSelectInput(session, "regression_endpoint", choices = choices)
  })
  
  # Create five-category treatment variable for regression
  regression_data <- reactive({
    if (input$regression_dataset == "assay") {
      df <- assay_filtered()
      
      df <- df %>%
        mutate(
          # Convert to lowercase for matching
          fiber_type_lower = tolower(fiber_type),
          treatment_lower = tolower(treatment),
          
          treatment_category = case_when(
            treatment_lower == "control" ~ "Control",
            fiber_type_lower == "cotton" & treatment_lower == "untreated" ~ "Untreated Cotton",
            fiber_type_lower == "cotton" & treatment_lower == "treated" ~ "Treated Cotton",
            fiber_type_lower == "pet" & treatment_lower == "untreated" ~ "Untreated Polyester",
            fiber_type_lower == "pet" & treatment_lower == "treated" ~ "Treated Polyester",
            TRUE ~ NA_character_
          ),
          treatment_category = factor(
            treatment_category,
            levels = c("Control", "Untreated Cotton", "Treated Cotton",
                       "Untreated Polyester", "Treated Polyester")
          ),
          week = as.factor(week),
          fiber_concentration = as.factor(fiber_concentration),
          sample_type = as.factor(sample_type)
        ) %>%
        filter(assay_type == input$regression_endpoint) %>%
        droplevels()
      
      df$outcome <- df$calculated_concentration
      
    } else {
      df <- phys_filtered()
      
      df <- df %>%
        mutate(
          # Convert to lowercase for matching
          fiber_type_lower = tolower(fiber_type),
          treatment_lower = tolower(treatment),
          
          treatment_category = case_when(
            treatment_lower == "control" ~ "Control",
            fiber_type_lower == "cotton" & treatment_lower == "untreated" ~ "Untreated Cotton",
            fiber_type_lower == "cotton" & treatment_lower == "treated" ~ "Treated Cotton",
            fiber_type_lower == "pet" & treatment_lower == "untreated" ~ "Untreated Polyester",
            fiber_type_lower == "pet" & treatment_lower == "treated" ~ "Treated Polyester",
            TRUE ~ NA_character_
          ),
          treatment_category = factor(
            treatment_category,
            levels = c("Control", "Untreated Cotton", "Treated Cotton",
                       "Untreated Polyester", "Treated Polyester")
          ),
          week = as.factor(week),
          fiber_concentration = as.factor(fiber_concentration),
          tissue_type = as.factor(tissue_type)
        ) %>%
        filter(endpoint == input$regression_endpoint) %>%
        droplevels()
      
      df$outcome <- df$value
    }
    
    df
  })
  
  # Fit regression model
  regression_model <- eventReactive(input$run_regression, {
    df <- regression_data()
    req(nrow(df) > 0)
    
    # DEBUG: Print what we have
    cat("\n=== DEBUG INFO ===\n")
    cat("Number of rows:", nrow(df), "\n")
    cat("Treatment categories:", paste(unique(df$treatment_category), collapse=", "), "\n")
    cat("Weeks:", paste(unique(df$week), collapse=", "), "\n")
    cat("Fiber concentrations:", paste(unique(df$fiber_concentration), collapse=", "), "\n")
    
    # Count levels for each factor
    n_treatment <- length(unique(na.omit(df$treatment_category)))
    n_week <- length(unique(na.omit(df$week)))
    n_conc <- length(unique(na.omit(df$fiber_concentration)))
    
    # Validate minimum requirements
    if (n_treatment < 2) {
      stop("Cannot fit model: treatment_category has ", n_treatment, 
           " level(s), need at least 2. Please adjust your filters to include more groups.")
    }
    
    # Build formula based on user selections AND available data
    formula_parts <- "outcome ~ treatment_category"
    
    # Only include week if there are 2+ weeks
    if (n_week >= 2) {
      if (input$include_interaction) {
        formula_parts <- paste0(formula_parts, " * week")
      } else {
        formula_parts <- paste0(formula_parts, " + week")
      }
    } else {
      message("Only one week detected (", unique(df$week), "). Week effects and interactions excluded from model.")
    }
    
    # ONLY include concentration if it has 2+ levels
    if (input$include_concentration && n_conc >= 2) {
      formula_parts <- paste0(formula_parts, " + fiber_concentration")
      cat("Including fiber_concentration in model (", n_conc, " levels)\n")
    } else if (input$include_concentration && n_conc < 2) {
      cat("WARNING: Skipping fiber_concentration - only", n_conc, "level(s)\n")
    }
    
    # Include tissue/sample type if multiple are present
    if (input$regression_dataset == "assay" && length(unique(df$sample_type)) > 1) {
      formula_parts <- paste0(formula_parts, " + sample_type")
    } else if (input$regression_dataset == "physical" && length(unique(df$tissue_type)) > 1) {
      formula_parts <- paste0(formula_parts, " + tissue_type")
    }
    
    cat("Final formula:", formula_parts, "\n")
    cat("==================\n\n")
    
    formula_obj <- as.formula(formula_parts)
    
    # Fit model
    model <- lm(formula_obj, data = df)
    return(model)
  })
  
  
  
  
  # Model summary output
  output$model_summary <- renderPrint({
    model <- regression_model()
    summary(model)
  })
  
  # ANOVA table
  output$anova_table <- renderDT({
    model <- regression_model()
    anova_results <- Anova(model, type = "II")
    
    anova_df <- as.data.frame(anova_results) %>%
      tibble::rownames_to_column("Term") %>%
      mutate(across(where(is.numeric), ~round(.x, 4)))
    
    datatable(anova_df, options = list(pageLength = 10))
  })
  
  # Diagnostic plots
  output$diagnostic_plots <- renderPlot({
    model <- regression_model()
    par(mfrow = c(2, 2))
    plot(model)
  })
  
  # Calculate estimated marginal means
  emmeans_results <- eventReactive(list(input$run_emmeans, regression_model()), {
    model <- regression_model()
    df <- regression_data()
    
    has_tissue <- if (input$regression_dataset == "assay") {
      "sample_type" %in% names(df) && length(unique(df$sample_type)) > 1
    } else {
      "tissue_type" %in% names(df) && length(unique(df$tissue_type)) > 1
    }
    has_week <- "week" %in% names(df) && length(unique(na.omit(df$week))) > 1
    tissue_var <- if (input$regression_dataset == "assay") "sample_type" else "tissue_type"
    
    # Choose formula safely based on availability
    if (input$emmeans_by == "treatment_category") {
      emm <- emmeans(model, ~ treatment_category)
    } else if (input$emmeans_by == "week" && has_week) {
      emm <- emmeans(model, ~ week)
    } else if (input$emmeans_by == "both" && has_week) {
      emm <- emmeans(model, ~ treatment_category | week)
    } else if (input$emmeans_by == "tissue" && has_tissue) {
      emm <- emmeans(model, as.formula(paste0("~ treatment_category | ", tissue_var)))
    } else if (input$emmeans_by == "all" && has_tissue && has_week) {
      emm <- emmeans(model, as.formula(paste0("~ treatment_category | week + ", tissue_var)))
    } else {
      # Fallback for single-week or no-tissue cases
      emm <- emmeans(model, ~ treatment_category)
    }
    emm
  })
  
  
  # EM Means plot (REGRESSION)
  output$emmeans_plot <- renderPlot({
    emm <- emmeans_results()
    
    # Get confidence intervals with standard column names
    df_emm <- summary(emm, infer = TRUE)  # This gives emmean, SE, df, lower.CL, upper.CL
    df_emm <- as.data.frame(df_emm)
    
    # Debug: check what columns we actually have
    cat("EM means columns:", paste(names(df_emm), collapse=", "), "\n")
    
    # Detect available grouping columns
    has_week <- "week" %in% names(df_emm) && length(unique(na.omit(df_emm$week))) > 1
    has_treat <- "treatment_category" %in% names(df_emm)
    
    # Handle different possible column names for confidence intervals
    if ("lower.CL" %in% names(df_emm)) {
      lower_col <- "lower.CL"
      upper_col <- "upper.CL"
    } else if ("asymp.LCL" %in% names(df_emm)) {
      lower_col <- "asymp.LCL"
      upper_col <- "asymp.UCL"
    } else {
      # Fallback: compute CIs manually from SE
      df_emm$lower.CL <- df_emm$emmean - 1.96 * df_emm$SE
      df_emm$upper.CL <- df_emm$emmean + 1.96 * df_emm$SE
      lower_col <- "lower.CL"
      upper_col <- "upper.CL"
    }
    
    # Build a safe EM means plot
    p <- ggplot(df_emm, aes(
      x = if (has_treat) treatment_category else factor(1),
      y = emmean
    )) +
      geom_point(size = 3, color = "#2c7fb8") +
      geom_errorbar(aes(ymin = .data[[lower_col]], ymax = .data[[upper_col]]), 
                    width = 0.2, color = "#2c7fb8") +
      coord_flip() +
      theme_minimal() +
      labs(
        title = "Estimated Marginal Means with 95% Confidence Intervals",
        x = "Group",
        y = "Estimated Mean"
      ) +
      theme(axis.text.y = element_text(size = 10))
    
    # If multiple weeks are present, facet
    if (has_week) {
      p <- p + facet_wrap(~ week, scales = "free_y")
    }
    
    print(p)
  })
  
  output$emmeans_table <- renderDT({
    emm <- emmeans_results()
    emm_df <- as.data.frame(summary(emm, infer = TRUE)) %>%
      mutate(across(where(is.numeric), ~signif(.x, 5)))
    datatable(emm_df, options = list(pageLength = 10))
  })
  
  
  # Pairwise comparisons for REGULAR REGRESSION
  output$pairwise_table <- renderDT({
    emm <- emmeans_results()
    adj_method <- input$adjustment_method
    
    # Compute comparisons
    if (input$comparison_type == "pairwise") {
      pairs_result <- pairs(emm, adjust = adj_method)
    } else if (input$comparison_type == "control") {
      if (adj_method == "dunnett") {
        pairs_result <- contrast(emm, method = "trt.vs.ctrl", ref = "Control", 
                                 adjust = "dunnett")
      } else {
        pairs_result <- contrast(emm, method = "trt.vs.ctrl", ref = "Control", 
                                 adjust = adj_method)
      }
    } else {
      pairs_result <- pairs(emm, adjust = adj_method)
    }
    
    # Convert and add Significant column
    pairs_df <- as.data.frame(pairs_result) %>%
      mutate(
        across(where(is.numeric), ~signif(.x, 5)),
        Significant = if_else(p.value < 0.05, "Yes", "No"),  # ADD THIS
        Adjustment = adj_method
      )
    
    datatable(pairs_df,
              options = list(pageLength = 20, scrollX = TRUE),
              filter = "top") %>%
      formatStyle("Significant",
                  backgroundColor = styleEqual(c("Yes", "No"),
                                               c("lightgreen", "white")))
  })
  
  
  # Download regression results
  output$download_regression <- downloadHandler(
    filename = function() {
      paste0("regression_analysis_", input$regression_endpoint, "_", 
             Sys.Date(), ".xlsx")
    },
    content = function(file) {
      model <- regression_model()
      emm <- emmeans_results()
      
      # Create workbook with multiple sheets
      wb_list <- list(
        "Model_Summary" = broom::tidy(model),
        "ANOVA" = broom::tidy(Anova(model, type = "II")),
        "EM_Means" = as.data.frame(emm),
        "Pairwise_Comparisons" = as.data.frame(pairs(emm, adjust = "tukey"))
      )
      
      writexl::write_xlsx(wb_list, file)
    }
  )
  
  # END OF REGRESSION ANALYSIS MODULE
  # ============================================================================

  # ============================================================================
  # MIXED EFFECTS ANALYSIS MODULE
  # ============================================================================
  
  # Dynamic endpoint selection for lmer
  observe({
    if (input$lmer_dataset == "assay") {
      choices <- unique(assay_filtered()$assay_type)
    } else {
      choices <- unique(phys_filtered()$endpoint)
    }
    updateSelectInput(session, "lmer_endpoint", choices = choices)
  })
  
  # Prepare data for mixed model
  lmer_data <- reactive({
    if (input$lmer_dataset == "assay") {
      df <- assay_filtered()
      
      # Add experiment batch variable
      df <- create_experiment_batch(df)
      
      df <- df %>%
        mutate(
          fiber_type_lower = tolower(fiber_type),
          treatment_lower = tolower(treatment),
          
          treatment_category = case_when(
            treatment_lower == "control" ~ "Control",
            fiber_type_lower == "cotton" & treatment_lower == "untreated" ~ "Untreated Cotton",
            fiber_type_lower == "cotton" & treatment_lower == "treated" ~ "Treated Cotton",
            fiber_type_lower == "pet" & treatment_lower == "untreated" ~ "Untreated Polyester",
            fiber_type_lower == "pet" & treatment_lower == "treated" ~ "Treated Polyester",
            TRUE ~ NA_character_
          ),
          treatment_category = factor(
            treatment_category,
            levels = c("Control", "Untreated Cotton", "Treated Cotton",
                       "Untreated Polyester", "Treated Polyester")
          ),
          week = as.factor(week),
          fiber_concentration = as.factor(fiber_concentration),
          sample_type = as.factor(sample_type),
          tank = as.factor(tank)
        ) %>%
        filter(assay_type == input$lmer_endpoint) %>%
        droplevels()
      
      df$outcome <- df$calculated_concentration
      df$grouping_var <- df$tank  # Tank is the sample unit
      
    } else {
      # Physical data
      df <- phys_filtered() %>%
        filter(endpoint == input$lmer_endpoint)
      
      # Add experiment batch variable
      df <- create_experiment_batch(df)
      
      # For mf_counts, extract individual ID from sample (e.g., "1_1" -> "1")
      if (input$lmer_endpoint == "mf_counts") {
        df <- df %>%
          mutate(
            individual_id = str_extract(sample, "^[0-9]+"),
            tissue_replicate = sample  # Full identifier for nested effect
          )
      } else {
        # For other endpoints, sample = individual = tank
        df <- df %>%
          mutate(
            individual_id = as.character(sample),
            tank_id = as.character(sample)  # Sample IS the tank for physical data
          )
      }
      
      df <- df %>%
        mutate(
          fiber_type_lower = tolower(fiber_type),
          treatment_lower = tolower(treatment),
          
          treatment_category = case_when(
            treatment_lower == "control" ~ "Control",
            fiber_type_lower == "cotton" & treatment_lower == "untreated" ~ "Untreated Cotton",
            fiber_type_lower == "cotton" & treatment_lower == "treated" ~ "Treated Cotton",
            fiber_type_lower == "pet" & treatment_lower == "untreated" ~ "Untreated Polyester",
            fiber_type_lower == "pet" & treatment_lower == "treated" ~ "Treated Polyester",
            TRUE ~ NA_character_
          ),
          treatment_category = factor(
            treatment_category,
            levels = c("Control", "Untreated Cotton", "Treated Cotton",
                       "Untreated Polyester", "Treated Polyester")
          ),
          week = as.factor(week),
          fiber_concentration = as.factor(fiber_concentration),
          tissue_type = as.factor(tissue_type),
          individual_id = as.factor(individual_id)
        ) %>%
        droplevels()
      
      df$outcome <- df$value
      df$grouping_var <- df$individual_id  # Sample/individual/tank are the same
    }
    
    df
  })
  
  
  # Fit mixed effects model
  lmer_model <- eventReactive(input$run_lmer, {
    df <- lmer_data()
    req(nrow(df) > 0)
    
    # Validation
    n_treatment <- length(unique(na.omit(df$treatment_category)))
    n_week <- length(unique(na.omit(df$week)))
    n_groups <- length(unique(na.omit(df$grouping_var)))
    n_batch <- length(unique(na.omit(df$experiment_batch)))
    
    if (n_treatment < 2) stop("Need at least 2 treatment categories")
    if (n_groups < 3) stop("Need at least 3 grouping levels for random effects")
    
    # Build fixed effects formula
    fixed_formula <- "outcome ~ treatment_category"
    
    # Only include week if there are 2+ weeks
    if (n_week >= 2) {
      if (input$lmer_include_interaction) {
        fixed_formula <- paste0(fixed_formula, " * week")
      } else {
        fixed_formula <- paste0(fixed_formula, " + week")
      }
    } else {
      message("Only one week detected. Week effects and interactions excluded from model.")
    }
    
    if (input$lmer_include_concentration && length(unique(df$fiber_concentration)) >= 2) {
      fixed_formula <- paste0(fixed_formula, " + fiber_concentration")
    }
    
    # Build random effects formula based on selection
    if (input$random_structure == "intercept") {
      random_formula <- "(1|grouping_var)"
      
    } else if (input$random_structure == "batch") {
      if (n_batch < 2) {
        message("WARNING: experiment_batch has only ", n_batch, " level(s). Using simple random intercept instead.")
        random_formula <- "(1|grouping_var)"
      } else {
        random_formula <- "(1|grouping_var) + (1|experiment_batch)"
      }
      
    } else if (input$random_structure == "nested") {
      if (input$lmer_dataset == "physical" && input$lmer_endpoint == "mf_counts") {
        random_formula <- "(1|grouping_var/tissue_replicate)"
      } else {
        warning("Nested structure only applicable to mf_counts. Using simple random intercept.")
        random_formula <- "(1|grouping_var)"
      }
      
    } else if (input$random_structure == "full") {
      if (n_batch < 2) {
        message("WARNING: experiment_batch has only ", n_batch, " level(s). Excluding batch from model.")
        if (input$lmer_dataset == "physical" && input$lmer_endpoint == "mf_counts") {
          random_formula <- "(1|grouping_var/tissue_replicate)"
        } else {
          random_formula <- "(1|grouping_var)"
        }
      } else {
        if (input$lmer_dataset == "physical" && input$lmer_endpoint == "mf_counts") {
          random_formula <- "(1|experiment_batch) + (1|grouping_var/tissue_replicate)"
        } else {
          random_formula <- "(1|experiment_batch) + (1|grouping_var)"
        }
      }
    }
    
    # Complete formula
    full_formula <- paste(fixed_formula, "+", random_formula)
    message("Fitting model: ", full_formula)
    
    # Fit model
    model <- lmer(as.formula(full_formula), data = df)
    
    return(model)
  })
  
  
  
  # Model summary
  output$lmer_summary <- renderPrint({
    model <- lmer_model()
    summary(model)
  })
  
  # Fixed effects ANOVA
  output$lmer_anova_table <- renderDT({
    model <- lmer_model()
    anova_results <- anova(model)
    
    anova_df <- as.data.frame(anova_results) %>%
      tibble::rownames_to_column("Term") %>%
      mutate(across(where(is.numeric), ~signif(.x, 5)))
    
    datatable(anova_df, options = list(pageLength = 10))
  })
  
  # Random effects summary
  output$lmer_random_effects <- renderPrint({
    model <- lmer_model()
    cat("Random Effects Variance Components:\n\n")
    print(VarCorr(model))
    cat("\nNumber of observations:", nobs(model), "\n")
    cat("Number of grouping levels:", 
        length(unique(lmer_data()$grouping_var)), "\n")
  })
  
  # Diagnostic plots
  output$lmer_diagnostics <- renderPlot({
    model <- lmer_model()
    par(mfrow = c(2, 2))
    
    # Residuals vs Fitted
    plot(fitted(model), residuals(model),
         xlab = "Fitted Values", ylab = "Residuals",
         main = "Residuals vs Fitted")
    abline(h = 0, col = "red", lty = 2)
    
    # Q-Q plot
    qqnorm(residuals(model), main = "Normal Q-Q Plot")
    qqline(residuals(model), col = "red")
    
    # Scale-Location
    plot(fitted(model), sqrt(abs(residuals(model))),
         xlab = "Fitted Values", ylab = "Sqrt(|Residuals|)",
         main = "Scale-Location")
    
    # Random effects Q-Q
    re <- ranef(model)[[1]][,1]
    qqnorm(re, main = "Random Effects Q-Q")
    qqline(re, col = "red")
  })
  
  # EM Means
  lmer_emmeans_results <- eventReactive(input$run_lmer_emmeans, {
    model <- lmer_model()
    # Get names of fixed effects actually in the model
    fe_terms <- attr(terms(model), "term.labels")
    
    has_week <- "week" %in% fe_terms
    has_treat <- "treatment_category" %in% fe_terms
    
    # Choose the safest emmeans specification
    if (input$lmer_emmeans_by == "treatment_category" && has_treat) {
      emm <- emmeans(model, ~ treatment_category)
    } else if (input$lmer_emmeans_by == "week" && has_week) {
      emm <- emmeans(model, ~ week)
    } else if (input$lmer_emmeans_by == "both" && has_week && has_treat) {
      emm <- emmeans(model, ~ treatment_category | week)
    } else {
      # Fallback: whatever fixed effect exists
      target <- if (has_treat) "treatment_category" else if (has_week) "week" else fe_terms[1]
      message("EMMeans fallback using term: ", target)
      emm <- emmeans(model, as.formula(paste("~", target)))
    }
    emm
  })
  
  
  # EM Means plot for MIXED EFFECTS
  output$lmer_emmeans_plot <- renderPlot({
    emm <- lmer_emmeans_results()
    
    df_emm <- as.data.frame(summary(emm, infer = TRUE))
    if ("asymp.LCL" %in% names(df_emm)) {
      df_emm <- df_emm %>% rename(lower.CL = asymp.LCL, upper.CL = asymp.UCL)
    }
    df_emm <- df_emm %>% filter(!is.na(emmean))
    
    if (nrow(df_emm) == 0) {
      plot.new()
      text(0.5, 0.5, "EM means not available after simplifying model.\nTry deselecting fiber_concentration or include more data.", cex = 1)
      return()
    }
    
    # Check we have required columns
    if (!all(c("treatment_category", "emmean", "lower.CL", "upper.CL") %in% names(df_emm))) {
      plot.new()
      text(0.5, 0.5, "Error: Missing required columns in emmeans output", cex = 1.2)
      return()
    }
    
    # Build plot
    p <- ggplot(df_emm, aes(x = treatment_category, y = emmean)) +
      geom_point(size = 3, color = "#2c7fb8") +
      geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                    width = 0.25, linewidth = 1, color = "#2c7fb8") +
      coord_flip() +
      theme_minimal(base_size = 12) +
      labs(
        title = "Estimated Marginal Means (Mixed Model)",
        subtitle = "Error bars show 95% confidence intervals",
        x = NULL,
        y = "Estimated Mean"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 11)
      )
    
    # Add faceting if week exists with >1 level
    if ("week" %in% names(df_emm) && length(unique(df_emm$week)) > 1) {
      p <- p + facet_wrap(~ week, scales = "free_y")
    }
    
    print(p)
  })
  
  # EM Means table
  output$lmer_emmeans_table <- renderDT({
    emm <- lmer_emmeans_results()
    emm_df <- as.data.frame(summary(emm, infer = TRUE)) %>%
      mutate(across(where(is.numeric), ~signif(.x, 5)))
    datatable(emm_df, options = list(pageLength = 10))
  })
  
  
  # Pairwise comparisons for MIXED EFFECTS
  output$lmer_pairwise_table <- renderDT({
    emm <- lmer_emmeans_results()
    adj_method <- input$lmer_adjustment_method
    
    # Compute comparisons based on type and adjustment
    if (input$lmer_comparison_type == "pairwise") {
      pairs_result <- pairs(emm, adjust = adj_method)
    } else if (input$lmer_comparison_type == "control") {
      if (adj_method == "dunnett") {
        pairs_result <- contrast(emm, method = "trt.vs.ctrl", ref = "Control", 
                                 adjust = "dunnett")
      } else {
        pairs_result <- contrast(emm, method = "trt.vs.ctrl", ref = "Control", 
                                 adjust = adj_method)
      }
    } else {
      pairs_result <- pairs(emm, adjust = adj_method)
    }
    
    # Convert to data frame and add significance column
    pairs_df <- as.data.frame(pairs_result) %>%
      mutate(
        across(where(is.numeric), ~signif(.x, 5)),
        Significant = if_else(p.value < 0.05, "Yes", "No"),  # ADD THIS LINE
        Adjustment = adj_method
      )
    
    # Display table
    datatable(pairs_df,
              options = list(pageLength = 20, scrollX = TRUE),
              filter = "top") %>%
      formatStyle("Significant",
                  backgroundColor = styleEqual(c("Yes", "No"),
                                               c("lightgreen", "white")))
  })
  
  
  # Download handler
  output$download_lmer <- downloadHandler(
    filename = function() {
      paste0("mixed_model_", input$lmer_endpoint, "_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      model <- lmer_model()
      emm <- lmer_emmeans_results()
      
      wb_list <- list(
        "Fixed_Effects" = broom.mixed::tidy(model, effects = "fixed"),
        "Random_Effects" = broom.mixed::tidy(model, effects = "ran_pars"),
        "ANOVA" = broom::tidy(anova(model)),
        "EM_Means" = as.data.frame(emm),
        "Pairwise" = as.data.frame(pairs(emm, adjust = input$lmer_adjustment_method))
      )
      
      writexl::write_xlsx(wb_list, file)
    }
  )
  
} 
shinyApp(ui, server)
