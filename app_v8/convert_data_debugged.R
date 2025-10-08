# Data Conversion Script - DEBUGGED VERSION
# Run this ONCE locally to optimize data loading

library(readxl)
library(janitor)
library(tidyverse)

# Set working directory (update this path)
setwd("C:/Users/Kevin/Documents/ECOTOX/app_v8")

message("🔄 Converting Excel files to RDS for faster loading...")

# Create data directory if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
  message("✅ Created data/ directory")
}

# Check if Excel files exist
required_files <- c("Assay_endpoint_DATA_cleaned.xlsx", "physical_endpoint_master_sheet.xlsx", "tissue_weights_clean.xlsx")
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  message("❌ Missing files: ", paste(missing_files, collapse = ", "))
  message("Current directory: ", getwd())
  message("Files in directory: ", paste(list.files(), collapse = ", "))
  stop("Please ensure Excel files are in the correct directory")
}

# ============================================================================
# 1. CONVERT ASSAY DATA
# ============================================================================

message("🔄 Processing assay data...")

tryCatch({
  assay_raw <- read_excel("Assay_endpoint_DATA_cleaned.xlsx") %>% 
    clean_names() %>%
    mutate(calculated_concentration = suppressWarnings(as.numeric(calculated_concentration)))

  # Apply all static transformations
  assay_processed <- assay_raw %>%
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
        tank <= 12 & !(tank %in% c(1, 2, 3)) ~ "Untreated",
        tank > 12 ~ "Treated",
        TRUE ~ NA_character_
      ),
      fiber_group = if_else(
        treatment == "Control" & fiber_concentration == "0", 
        "Control", 
        paste(fiber_type, treatment)
      )
    ) %>%
    dplyr::select(-week_start, -week_end, -well_number) %>%
    arrange(fiber_type, week, tank) %>%
    as_tibble()

  # Save processed assay data
  saveRDS(assay_processed, "data/assay_data.rds")
  message("✅ Assay data saved to data/assay_data.rds")

}, error = function(e) {
  message("❌ Error processing assay data: ", e$message)
  stop("Failed to process assay data")
})

# ============================================================================
# 2. CONVERT PHYSICAL DATA
# ============================================================================

message("🔄 Processing physical data...")

tryCatch({
  physical_raw <- read_excel("physical_endpoint_master_sheet.xlsx") %>%
    clean_names()

  # Process with mf_counts averaging
  physical_processed <- physical_raw %>%
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
      fiber_group = if_else(
        treatment == "Control" & fiber_concentration == "0",
        "Control", paste(fiber_type, treatment)
      )
    ) %>%
    # Process mf_counts averaging efficiently
    split(.$endpoint) %>%
    map_dfr(function(.x) {
      ep <- unique(.x$endpoint)
      if (identical(ep, "mf_counts")) {
        .x %>%
          mutate(
            sample_base = str_extract(sample, "^[0-9]+"),
            replicate_id = str_extract(sample, "[0-9]+$")
          ) %>%
          group_by(fiber_type, week, tissue_type, endpoint, tank,
                  treatment, fiber_concentration, sample_base) %>%
          summarise(
            value = mean(value, na.rm = TRUE),
            n_reps = n(),
            fiber_group = first(fiber_group),
            .groups = "drop"
          ) %>%
          mutate(sample = sample_base) %>%
          dplyr::select(-sample_base)
      } else {
        .x
      }
    }) %>%
    as_tibble()

  # Save processed physical data
  saveRDS(physical_processed, "data/physical_data.rds")
  message("✅ Physical data saved to data/physical_data.rds")

}, error = function(e) {
  message("❌ Error processing physical data: ", e$message)
  stop("Failed to process physical data")
})

# ============================================================================
# 3. CONVERT TISSUE WEIGHTS
# ============================================================================

message("🔄 Processing tissue weights...")

tissue_weights <- tryCatch({
  read_excel("tissue_weights_clean.xlsx") %>% 
    clean_names()
}, error = function(e) {
  message("⚠️ Warning: tissue_weights_clean.xlsx not found, using empty dataframe")
  data.frame()
})

saveRDS(tissue_weights, "data/tissue_weights.rds")
message("✅ Tissue weights saved to data/tissue_weights.rds")

# ============================================================================
# 4. VERIFY FILES WERE CREATED
# ============================================================================

message("\n📁 Verifying created files...")

# Check if all RDS files exist
rds_files <- c("data/assay_data.rds", "data/physical_data.rds", "data/tissue_weights.rds")
created_files <- file.exists(rds_files)

if (all(created_files)) {
  message("✅ All RDS files created successfully!")
  
  # Show file sizes
  get_file_size <- function(file) {
    if (file.exists(file)) {
      size_mb <- round(file.info(file)$size / 1024 / 1024, 2)
      paste0(size_mb, " MB")
    } else {
      "File not found"
    }
  }
  
  message("📊 File sizes:")
  message("  Assay data: ", get_file_size("data/assay_data.rds"))
  message("  Physical data: ", get_file_size("data/physical_data.rds"))
  message("  Tissue weights: ", get_file_size("data/tissue_weights.rds"))
  
} else {
  missing_rds <- rds_files[!created_files]
  message("❌ Failed to create: ", paste(missing_rds, collapse = ", "))
}

# ============================================================================
# 5. CREATE OPTIMIZED GLOBAL.R
# ============================================================================

global_optimized <- '
# OPTIMIZED global.R - Uses RDS files for fast loading

# ===== LIBRARIES =====
library(shiny)
library(tidyverse)
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
library(MASS)
library(performance)

# ===== UTILITY FUNCTIONS =====
`%||%` <- function(a, b) if (!is.null(a)) a else b
filter_choices <- function(x) unlist(x[!startsWith(as.character(x), "#")])

# ===== LOAD OPTIMIZED DATA =====
message("Loading optimized data...")

# Load pre-processed data (much faster than Excel)
final_data <- readRDS("data/assay_data.rds")
physical_master <- readRDS("data/physical_data.rds")
tissue_weights <- readRDS("data/tissue_weights.rds")

# Extract UI choices
week_choices_assay <- sort(unique(na.omit(final_data$week)))
endpoint_choices_assay <- sort(unique(final_data$assay_type))
endpoint_choices_physical <- sort(unique(physical_master$endpoint))

# Create lookup tables
tank_fiber_lookup <- final_data %>%
  select(tank, fiber_type, fiber_concentration, treatment) %>%
  distinct() %>%
  arrange(tank)

# Load config
config <- tryCatch({
  yaml::read_yaml("config.yml")
}, error = function(e) {
  list(
    assays = endpoint_choices_assay,
    fibers = c("cotton", "pet"),
    samples = c("tissue", "gland", "gills")
  )
})

message("✅ Global data loading complete!")
message("📊 Data summary:")
message("  Assay data: ", nrow(final_data), " rows")
message("  Physical data: ", nrow(physical_master), " rows")
message("  Memory usage: ", round(as.numeric(object.size(final_data) + object.size(physical_master))/1024^2, 1), " MB")
'

writeLines(global_optimized, "global_optimized.R")
message("📝 Created global_optimized.R")

message("\n🎉 DATA CONVERSION COMPLETE!")
message("🚀 Next steps:")
message("1. Replace your current global.R with global_optimized.R")
message("2. Test your app with: runApp()")
message("3. Deploy to shinyapps.io when ready")