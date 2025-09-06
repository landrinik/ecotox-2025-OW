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

filter_choices <- function(x) {
  unlist(x[!startsWith(as.character(x), "#")])
}

tissue_weights <- read_excel(here("tissue_weights_clean.xlsx")) %>% clean_names()
final_data <- read_excel(here("final_long_table_normalized.xlsx")) %>% clean_names() %>%
  mutate(calculated_concentration = as.numeric(calculated_concentration))

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
        well_number > 6 & well_letter %in% c("A", "B") ~ (well_number - 7) * 4 + 4,
        well_number > 6 & well_letter %in% c("C", "D") ~ (well_number - 7) * 4 + 5,
        well_number > 6 & well_letter %in% c("E", "F") ~ (well_number - 7) * 4 + 6,
        well_number > 6 & well_letter %in% c("G", "H") ~ (well_number - 7) * 4 + 7,
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
        tank <= 12 & !(tank %in% c(1, 2, 3)) ~ "Treated",
        tank > 12 ~ "Untreated",
        TRUE ~ NA_character_
      ),
      fiber_group = if_else(treatment == "Control" & fiber_concentration == "0", "Control",
                            paste(fiber_type, treatment)
      )
    )
}

final_data <- annotate_common(final_data)

# --- UI ---

ui <- fluidPage(
  useShinyjs(),
  titlePanel("ECOTOX Assay Dashboard"),
  actionButton("toggleSidebar", "Toggle Sidebar"),
  fluidRow(
    column(
      width = 3,
      id = "sidebar",
      h4("Processing mode:"),
      radioButtons("mode", NULL,
                   choices = c(
                     "Baseline (filter only)" = "baseline",
                     "Well-level CV (wrangle_v3_avg)" = "wrangle",
                     "Inter-assay CV (inter_assay_CV)" = "inter"
                   ), selected = "baseline"
      ),
      br(),
      selectizeInput("assay", "Assay Type:",
                     choices = filter_choices(config$assays) %||% unique(final_data$assay_type),
                     multiple = TRUE
      ),
      selectizeInput("fiber", "Fiber Type:",
                     choices = filter_choices(config$fibers) %||% unique(final_data$fiber_type),
                     multiple = TRUE
      ),
      selectizeInput("sample_type", "Sample Type:",
                     choices = filter_choices(config$samples) %||% unique(final_data$sample_type),
                     multiple = TRUE
      ),
      selectizeInput("week", "Week:",
                     choices = sort(unique(final_data$week)),
                     multiple = TRUE
      ),
      selectizeInput("treatment", "Treatment:",
                     choices = c("Treated", "Untreated", "Control"),
                     multiple = TRUE
      ),
      selectizeInput("fiber_concentration", "Fiber Concentration:",
                     choices = c("0", "100", "1000", "10000"),
                     multiple = TRUE
      ),
      conditionalPanel(
        condition = "input.mode != 'inter'",
        selectizeInput("plate_replicate", "Plate Replicate:",
                       choices = unique(final_data$plate_replicate),
                       multiple = TRUE
        )
      ),
      sliderInput("cv_threshold", "CV Threshold (%)", min = 0, max = 100, value = 15),
      actionButton("run_analysis", "Run Analysis"),
      br(),
      downloadButton("download_data", "Download Results")
    ),
    column(
      width = 9,
      DTOutput("table", width = "100%"),
      br(),
      plotOutput("dist_plot", width = "100%", height = "500px"),
      br(),
      plotOutput("dot_plot", width = "100%", height = "400px"),
      verbatimTextOutput("info_text")
    )
  )
)

# --- Server ---

server <- function(input, output, session) {
  
  observeEvent(input$toggleSidebar, {
    toggle("sidebar")
  })
  
  # Filtering
  base_filtered <- eventReactive(input$run_analysis, {
    df <- final_data
    if (!is.null(input$assay) && length(input$assay) > 0)
      df <- df %>% filter(assay_type %in% input$assay)
    if (!is.null(input$fiber) && length(input$fiber) > 0)
      df <- df %>% filter(fiber_type %in% input$fiber)
    if (!is.null(input$sample_type) && length(input$sample_type) > 0)
      df <- df %>% filter(sample_type %in% input$sample_type)
    if (!is.null(input$week) && length(input$week) > 0)
      df <- df %>% filter(week %in% as.numeric(input$week))
    if (!is.null(input$fiber_concentration) && length(input$fiber_concentration) > 0) {
      df <- df %>% filter(fiber_concentration %in% input$fiber_concentration)
    }
    if (!is.null(input$treatment) && length(input$treatment) > 0) {
      if ("Control" %in% input$treatment) {
        df <- df %>% filter(fiber_concentration == "0" | treatment %in% setdiff(input$treatment, "Control"))
      } else {
        df <- df %>% filter(treatment %in% input$treatment)
      }
    }
    if (input$mode != "inter" && !is.null(input$plate_replicate) && length(input$plate_replicate) > 0) {
      df <- df %>% filter(plate_replicate %in% input$plate_replicate)
    }
    df
  })
  
  #test summaries
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
        vals = list(calculated_concentration[is.finite(calculated_concentration)]),
        mean_activity = ifelse(length(vals[[1]]) == 0, NA_real_, mean(vals[[1]], na.rm = TRUE)),
        sd_activity_raw = ifelse(length(vals[[1]]) == 0, NA_real_, sd(vals[[1]], na.rm = TRUE)),
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
  
  # # Summaries
  # summarize_baseline <- function(df, threshold) {
  #   df %>%
  #     group_by(fiber_group, assay_type, week, sample_type, tank, well_letter,
  #              fiber_concentration, well_name, plate_replicate) %>%
  #     summarise(
  #       mean_activity = mean(calculated_concentration, na.rm = TRUE),
  #       sd_activity = sd(calculated_concentration, na.rm = TRUE),
  #       cv_percent = round(if_else(mean_activity != 0, (sd_activity / mean_activity) * 100, NA_real_), 2),
  #       cv_flag = cv_percent > threshold,
  #       .groups = "drop"
  #     ) %>%
  #     arrange(assay_type, fiber_group, week)
  # }
  # 
  # summarize_wrangle <- function(df, threshold) {
  #   df %>%
  #     left_join(tissue_weights,
  #               by = c("assay_type", "fiber_type", "sample_type", "sample_week", "well_name")) %>%
  #     group_by(fiber_group, assay_type, week, sample_type, tank, well_letter,
  #              fiber_concentration, well_name) %>%
  #     summarise(
  #       mean_activity = mean(calculated_concentration, na.rm = TRUE),
  #       sd_activity_raw = sd(calculated_concentration, na.rm = TRUE),
  #       .groups = "drop"
  #     ) %>%
  #     mutate(
  #       sd_activity = round(sd_activity_raw, 4),
  #       cv_activity = if_else(mean_activity != 0, (sd_activity / mean_activity) * 100, NA_real_),
  #       cv_percent = round(cv_activity, 2),
  #       cv_flag = cv_percent > threshold
  #     ) %>%
  #      dplyr::select(fiber_group, assay_type, week, sample_type, tank,
  #            well_letter, fiber_concentration, well_name,
  #            mean_activity, sd_activity, cv_percent, cv_flag)
  # }
  # 
  # summarize_inter <- function(df, threshold) {
  #   df %>%
  #     group_by(fiber_group, assay_type, week, sample_type, fiber_concentration, plate_replicate) %>%
  #     summarise(
  #       mean_activity = mean(calculated_concentration, na.rm = TRUE),
  #       sd_activity = sd(calculated_concentration, na.rm = TRUE),
  #       n_samples = n(),
  #       cv = if_else(mean_activity == 0, NA_real_, 100 * sd_activity / mean_activity),
  #       .groups = "drop"
  #     ) %>%
  #     mutate(
  #       cv = round(cv, 2),
  #       cv_flag = cv > threshold
  #     )
  # }
  # 
  summarized <- reactive({
    df <- base_filtered()
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
  
  # Outputs
  output$table <- renderDT({
    df <- summarized()
    if (is.null(df)) return(datatable(data.frame(Message = "No data for selected filters / mode")))
    
    # Optionally for debugging: 
    print(str(df))   
    print(class(df))
    
    # Remove problematic columns
    is_simple_col <- sapply(df, function(col) is.atomic(col) && !is.list(col))
    df <- df[, is_simple_col, drop = FALSE]
    
    # If mean_activity_display exists, always use that for display
    if ("mean_activity_display" %in% names(df)) {
      df$mean_activity <- df$mean_activity_display
      df$mean_activity_display <- NULL
    }
    
    datatable(df, filter = "top", options = list(pageLength = 15))
  })
  
 
  # Boxplot with dots
  output$dist_plot <- renderPlot({
    df <- base_filtered()
    if (nrow(df) == 0) return(NULL)
    ggplot(df, aes(x = fiber_group, y = calculated_concentration, fill = fiber_concentration)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
      facet_wrap(~ assay_type) +
      theme_minimal() +
      labs(
        title = "Comparison of Calculated Concentrations by Group and Fiber Concentration",
        x = "Fiber Group",
        y = "Calculated Concentration",
        fill = "Fiber Concentration"
      ) +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  })

  # Dots-only jitter plot (no box)
  output$dot_plot <- renderPlot({
    df <- base_filtered()
    if (nrow(df) == 0) return(NULL)
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
  
  output$info_text <- renderText({
    df <- summarized()
    if (is.null(df)) return("No data for current filters / mode.")
    flagged <- sum(dplyr::coalesce(df$cv_flag, FALSE), na.rm = TRUE)
    total <- nrow(df)
    glue("Mode: {input$mode}. Showing {total} summary rows. {flagged} rows exceed CV threshold ({input$cv_threshold}%).")
  })
  
  output$download_data <- downloadHandler(
    filename = function() { paste0("ECOTOX_", input$mode, "_", Sys.Date(), ".xlsx") },
    content = function(file) {
      df <- summarized()
      writexl::write_xlsx(if (is.null(df)) data.frame(Message = "No data for selected filters / mode") else df, file)
    }
  )
}

shinyApp(ui, server)



