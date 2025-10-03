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

`%||%` <- function(a, b) if (!is.null(a)) a else b

# --- Load config and data ---

config <- yaml::read_yaml(here("config.yml"))
filter_choices <- function(x) unlist(x[!startsWith(as.character(x), "#")])
tissue_weights <- read_excel(here("tissue_weights_clean.xlsx")) %>% clean_names()

# -- Assay Data --

final_data_raw <- read_excel(here("Assay_endpoint_DATA_cleaned.xlsx")) %>% clean_names()

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
    select(-week_start, -week_end, -well_number)
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
        paste0(round(mean_activity, 2), "*"),
        as.character(round(mean_activity, 2))
      )
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
      h4("Results"),
      DTOutput("table"),
      br(),
      plotOutput("dist_plot"),
      br(),
      plotOutput("dot_plot"),
      br(),
      verbatimTextOutput("info")
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
        value = if_else(endpoint == "mf_counts", round(value), value) # round visible mf_counts for table
      ) %>%
      ungroup()
    
    if (input$mode_phys == "baseline") {
      # Return original values for all, but display rounded counts in table
      df
    } else {
      # Use the *corrected* values only for downstream summary (not for controls in table)
      df %>%
        mutate(use_val = if_else(endpoint == "mf_counts", value_corr, value)) %>%
        group_by(fiber_group, week, endpoint, tissue_type, fiber_concentration) %>%
        summarise(
          mean_value = round(mean(use_val, na.rm = TRUE)),
          sd_value = sd(use_val, na.rm = TRUE),
          n = n(),
          cv = ifelse(mean_value != 0, 100 * sd_value / mean_value, NA_real_),
          .groups = "drop"
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
    datatable(df, filter = "top", options = list(pageLength = 15))
  })
  
  
  # Distribution plot
  output$dist_plot <- renderPlot({
    df <- active_table()
    req(df, nrow(df) > 0)
    
    if (input$active_dataset == "assay") {
      x_axis_var <- input$x_axis_assay %||% "fiber_group"
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
}

shinyApp(ui, server)
