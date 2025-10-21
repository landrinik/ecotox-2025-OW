# ============================================================================
#                       convert_data_debugged.R
# ============================================================================
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
required_files <- c(
  "Assay_endpoint_DATA_cleaned.xlsx",
  "physical_endpoint_master_sheet.xlsx",
  "tissue_weights_clean.xlsx"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  message("❌ Missing files: ", paste(missing_files, collapse = ", "))
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
  
  saveRDS(assay_processed, "data/assay_data.rds")
  message("✅ Assay data saved to data/assay_data.rds")
  
}, error = function(e) {
  message("❌ Error processing assay data: ", e$message)
  stop("Failed to process assay data")
})

# ============================================================================
# 2. CONVERT PHYSICAL DATA - SIMPLIFIED WITH CONTROL CORRECTION
# ============================================================================

message("🔄 Processing physical data...")

tryCatch({
  physical_raw <- read_excel("physical_endpoint_master_sheet.xlsx") %>%
    clean_names()
  
  # Add treatment metadata
  physical_raw <- physical_raw %>%
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
    )
  
  # SEPARATE PROCESSING: mf_counts vs others
  # Process mf_counts separately (has tissue_type, needs averaging)
  mf_data <- physical_raw %>%
    filter(endpoint == "mf_counts")
  
  # Process other endpoints (no tissue_type, no averaging)
  other_data <- physical_raw %>%
    filter(endpoint != "mf_counts")
  
  # STEP 1: Control correction for mf_counts ONLY
  if (nrow(mf_data) > 0) {
    mf_corrected <- mf_data %>%
      group_by(fiber_type, week, tissue_type) %>%
      mutate(
        avg_control = mean(value[treatment == "Control"], na.rm = TRUE),
        n_controls = sum(treatment == "Control"),
        n_nonzero_controls = sum(treatment == "Control" & value > 0),
        value_corr = if_else(
          treatment != "Control",
          pmax(value - avg_control, 0),
          value
        )
      ) %>%
      ungroup()
    
    mf_corrected <- mf_corrected %>%
      # coerce any blanks like "" to NA, keep numeric zeros
      dplyr::mutate(
        value      = suppressWarnings(as.numeric(value)),
        value_corr = suppressWarnings(as.numeric(value_corr))
      ) %>%
      # drop only non-finite rows (NA/NaN/Inf); zeros are finite and kept
      dplyr::filter(is.finite(value) | is.finite(value_corr))
    
    # STEP 2: KEEP REPLICATES (no averaging) for mf_counts
    mf_noavg <- mf_corrected %>%
      # optional: parse a replicate_id if your sample encodes it (kept as-is if absent)
      dplyr::mutate(
        replicate_id = stringr::str_extract(sample, "(?<=_)\\d+"),
        value       = round(value, 1),
        value_corr  = round(value_corr, 1)
      )
    
  } else {
    mf_noavg <- dplyr::tibble()
  }
  
  # Combine back together (mf_counts with replicates + other endpoints)
  physical_processed <- dplyr::bind_rows(mf_noavg, other_data) %>%
    dplyr::arrange(fiber_type, week, endpoint)
  
  saveRDS(physical_processed, "data/physical_data.rds")
  message("✅ Physical data saved to data/physical_data.rds")
  message("   ✓ Control correction applied to mf_counts")
  
}, error = function(e) {
  message("❌ Error processing physical data: ", e$message)
  print(e)
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

rds_files <- c("data/assay_data.rds", "data/physical_data.rds", "data/tissue_weights.rds")
created_files <- file.exists(rds_files)

if (all(created_files)) {
  message("✅ All RDS files created successfully!")
  
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

message("\n✅ DATA CONVERSION COMPLETE!")