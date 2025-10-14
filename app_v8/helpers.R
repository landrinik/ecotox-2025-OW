# Data Processing Helper Functions
# Optimized for performance and clarity

# ===== EXPERIMENT BATCH CREATION =====
create_experiment_batch <- function(df) {
  df %>%
    dplyr::mutate(
      experiment_batch = dplyr::case_when(
        tolower(fiber_type) == "cotton" ~ "April_Cotton",
        tolower(fiber_type) == "pet" ~ "January_PET",
        TRUE ~ NA_character_
      ),
      experiment_batch = factor(experiment_batch)
    )
}

# ===== TREATMENT CATEGORIZATION =====
create_enhanced_treatment_categories <- function(df) {
  df %>%
    dplyr::mutate(
      fiber_type_lower = stringr::str_to_lower(fiber_type),
      treatment_lower = stringr::str_to_lower(treatment),
      # Create 4-category structure with controls as concentration 0
      treatment_category = dplyr::case_when(
        fiber_type_lower == "cotton" & treatment_lower %in% c("untreated", "control") ~ "Untreated Cotton",
        fiber_type_lower == "cotton" & treatment_lower == "treated" ~ "Treated Cotton",
        fiber_type_lower == "pet" & treatment_lower %in% c("untreated", "control") ~ "Untreated Polyester",
        fiber_type_lower == "pet" & treatment_lower == "treated" ~ "Treated Polyester",
        TRUE ~ NA_character_
      ),
      treatment_category = factor(treatment_category, levels = c(
        "Untreated Cotton", "Treated Cotton", "Untreated Polyester", "Treated Polyester"
      )),
      # Set fiber_concentration = 0 for controls
      fiber_concentration_numeric = dplyr::case_when(
        treatment_lower == "control" ~ 0,
        TRUE ~ as.numeric(fiber_concentration)
      )
    )
}

# ===== NORMALIZATION FUNCTIONS =====
normalize_controls_and_dose <- function(df) {
  # Standardize fiber_type
  df <- df %>%
    dplyr::mutate(
      fiber_type = stringr::str_to_lower(fiber_type),
      fiber_type = dplyr::if_else(
        fiber_type %in% c("cotton", "cot"), "cotton",
        dplyr::if_else(fiber_type %in% c("pet", "polyester"), "pet", fiber_type)
      ),
      fiber_type = factor(fiber_type, levels = c("cotton", "pet"))
    )
  
  # Map original treatment into chem_treatment (controls become Untreated)
  df <- df %>%
    dplyr::mutate(
      treatment_lower = stringr::str_to_lower(treatment),
      chem_treatment = dplyr::case_when(
        treatment_lower == "treated" ~ "treated",
        treatment_lower == "untreated" ~ "untreated",
        treatment_lower == "control" ~ "untreated", # controls are untreated by design
        TRUE ~ NA_character_
      ),
      chem_treatment = factor(chem_treatment, levels = c("untreated", "treated"))
    )
  
  # Build dose on log10 scale and control indicator
  df <- df %>%
    dplyr::mutate(
      # fiber_concentration must be numeric; coerce if needed
      fiber_concentration = suppressWarnings(as.numeric(fiber_concentration)),
      # Missing or non-numeric become 0 for safety
      fiber_concentration = dplyr::if_else(is.na(fiber_concentration), 0, fiber_concentration),
      is_control = as.integer(fiber_concentration == 0),
      dose_log10 = dplyr::if_else(
        fiber_concentration > 0,
        log10(fiber_concentration),
        0  # controls sit at 0 on log10 scale
      )
    )
  
  # Factor week for interactions
  if (!is.factor(df$week)) df$week <- factor(df$week)
  
  df
}

# ===== REPLICATE AVERAGING =====
# For assay data: average known replicate columns out of the grouping
average_assay_replicates <- function(df, outcome_col) {
  replicate_cols <- c("plate_replicate", "replicate", "rep", "plate")
  drop_cols <- c(replicate_cols, rlang::as_name(rlang::enquo(outcome_col)))
  
  df %>%
    dplyr::group_by(dplyr::across(tidyselect::all_of(setdiff(names(.), drop_cols)))) %>%
    dplyr::summarise(
      !!rlang::enquo(outcome_col) := mean({{ outcome_col }}, na.rm = TRUE),
      n_reps = dplyr::n(),
      .groups = "drop"
    )
}

# For physical endpoints: average "1_1, 1_2 → 1" by sample_base  
average_physical_replicates <- function(df) {
  df %>%
    dplyr::mutate(
      sample_base = stringr::str_extract(sample, "^[0-9]+"),
      replicate = stringr::str_extract(sample, "[0-9]+$")
    ) %>%
    dplyr::group_by(fiber_type, week, tissue_type, endpoint, sample_base) %>%
    dplyr::summarise(
      value = mean(value, na.rm = TRUE),
      n_reps = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(sample = sample_base) %>%
    dplyr::select(-sample_base)
}

# ===== DATA SUMMARIZATION FUNCTIONS =====
summarize_baseline <- function(df, threshold) {
  df %>%
    dplyr::group_by(fiber_group, assay_type, week, sample_type, tank, well_letter,
                   fiber_concentration, well_name, plate_replicate) %>%
    dplyr::summarise(
      vals = list(calculated_concentration[is.finite(calculated_concentration)]),
      mean_activity = ifelse(length(vals[[1]]) == 0, NA_real_, mean(vals[[1]], na.rm = TRUE)),
      sd_activity = ifelse(length(vals[[1]]) == 0, NA_real_, sd(vals[[1]], na.rm = TRUE)),
      has_inf = any(is.infinite(calculated_concentration)),
      has_zero = any(calculated_concentration == 0),
      has_na = any(is.na(calculated_concentration)),
      cv_percent = round(dplyr::if_else(mean_activity != 0 & !is.na(mean_activity),
                                       (sd_activity / mean_activity) * 100, NA_real_), 2),
      cv_flag = cv_percent > threshold,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # Create display version of mean_activity with asterisk if any flag true
      mean_activity_display = dplyr::if_else(
        has_inf | has_zero | has_na,
        paste0(round(mean_activity, 2), "*"),
        as.character(round(mean_activity, 2))
      )
    ) %>%
    # Remove columns you don't want to display after summarise
    dplyr::select(-has_inf, -has_zero, -has_na, -vals)
}

summarize_wrangle <- function(df, threshold) {
  df %>%
    dplyr::left_join(tissue_weights,
                    by = c("assay_type", "fiber_type", "sample_type", "sample_week", "well_name")) %>%
    dplyr::group_by(fiber_group, assay_type, week, sample_type, tank, well_letter,
                   fiber_concentration, well_name) %>%
    dplyr::summarise(
      vals = list(calculated_concentration[is.finite(calculated_concentration) & calculated_concentration != 0]),
      mean_activity = ifelse(length(vals[[1]]) == 0, NA_real_, mean(unlist(vals), na.rm = TRUE)),
      sd_activity_raw = ifelse(length(vals[[1]]) == 0, NA_real_, sd(unlist(vals), na.rm = TRUE)),
      has_inf = any(is.infinite(calculated_concentration)),
      has_zero = any(calculated_concentration == 0),
      has_na = any(is.na(calculated_concentration)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      sd_activity = round(sd_activity_raw, 4),
      cv_activity = dplyr::if_else(mean_activity != 0 & !is.na(mean_activity),
                                  (sd_activity / mean_activity) * 100,
                                  NA_real_),
      cv_percent = round(cv_activity, 2),
      cv_flag = cv_percent > threshold,
      mean_activity_display = dplyr::if_else(has_inf | has_na,
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
    dplyr::group_by(fiber_group, assay_type, week, sample_type, fiber_concentration, plate_replicate) %>%
    dplyr::summarise(
      vals = list(calculated_concentration[is.finite(calculated_concentration)]),
      mean_activity = ifelse(length(vals[[1]]) == 0, NA_real_, mean(vals[[1]], na.rm = TRUE)),
      sd_activity = ifelse(length(vals[[1]]) == 0, NA_real_, sd(vals[[1]], na.rm = TRUE)),
      n_samples = dplyr::n(),
      has_inf = any(is.infinite(calculated_concentration)),
      has_zero = any(calculated_concentration == 0),
      has_na = any(is.na(calculated_concentration)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      cv = dplyr::if_else(mean_activity == 0, NA_real_, 100 * sd_activity / mean_activity),
      cv = round(cv, 2),
      cv_flag = cv > threshold,
      mean_activity_display = dplyr::if_else(
        has_inf | has_zero | has_na,
        paste0(format(mean_activity, digits = 4, scientific = FALSE), "*"),
        format(mean_activity, digits = 4, scientific = FALSE)
      ),
      sd_activity = signif(sd_activity, 5)
    ) %>%
    dplyr::select(-has_inf, -has_zero, -has_na, -vals)
}

process_assay_data_for_recovery <- function(df) {
  df %>%
    dplyr::mutate(
      # Extract week components
      week_start = as.numeric(stringr::str_extract(sample_week, "^[0-9]+")),
      week_end = as.numeric(stringr::str_extract(sample_week, "[0-9]+$")),
      well_number = as.numeric(stringr::str_extract(well_name, "[0-9]+")),
      well_letter = stringr::str_extract(well_name, "[A-H]"),
      
      # Apply the EXACT same week splitting logic as convert_data_debugged.R
      week = dplyr::case_when(
        well_number <= 6 & well_letter %in% c("A", "B") ~ week_start,
        well_number <= 6 & well_letter %in% c("C", "D", "E", "F", "G", "H") ~ week_end,
        well_number > 6 ~ week_start,
        well_number > 6 ~ week_end,
        TRUE ~ NA_real_
      ),
      
      # Apply the EXACT same tank assignment logic as convert_data_debugged.R
      tank = dplyr::case_when(
        well_number <= 5 & well_letter %in% c("A", "B") ~ (well_number - 1) %/% 4 + 1,
        well_number <= 5 & well_letter %in% c("C", "D") ~ (well_number - 1) %/% 4 + 2,
        well_number <= 5 & well_letter %in% c("E", "F") ~ (well_number - 1) %/% 4 + 3,
        well_number <= 5 & well_letter %in% c("G", "H") ~ (well_number - 1) %/% 4 + 4,
        well_number > 6 & well_letter %in% c("A", "B") ~ 21,
        well_number > 6 & well_letter %in% c("C", "D") ~ 1,
        well_number > 6 & well_letter %in% c("E", "F") ~ 2,
        well_number > 6 & well_letter %in% c("G", "H") ~ 3,
        well_number > 6 & well_letter %in% c("A", "B") ~ (well_number - 7) %/% 4 + 4,
        well_number > 6 & well_letter %in% c("C", "D") ~ (well_number - 7) %/% 4 + 5,
        well_number > 6 & well_letter %in% c("E", "F") ~ (well_number - 7) %/% 4 + 6,
        well_number > 6 & well_letter %in% c("G", "H") ~ (well_number - 7) %/% 4 + 7,
        TRUE ~ NA_real_
      ),
      
      # Apply the EXACT same treatment logic as convert_data_debugged.R
      fiber_concentration = dplyr::case_when(
        tank %in% c(1, 2, 3) ~ "0",
        tank %in% c(4, 5, 6, 13, 14, 15) ~ "100",
        tank %in% c(7, 8, 9, 16, 17, 18) ~ "1000",
        tank %in% c(10, 11, 12, 19, 20, 21) ~ "10000",
        TRUE ~ NA_character_
      ),
      
      treatment = dplyr::case_when(
        tank %in% c(1, 2, 3) ~ "Control",
        tank > 12 & !tank %in% c(1, 2, 3) ~ "Untreated",
        tank <= 12 ~ "Treated",
        TRUE ~ NA_character_
      ),
      
      fiber_group = ifelse(treatment == "Control" | fiber_concentration == "0", 
                           "Control", 
                           paste(fiber_type, treatment))
    ) %>%
    # Remove temporary columns
    dplyr::select(-week_start, -week_end, -well_number, -well_letter) %>%
    # Remove rows with missing weeks or tanks
    dplyr::filter(!is.na(week), !is.na(tank)) %>%
    # Arrange like the original processing
    dplyr::arrange(fiber_type, week, tank) %>%
    tibble::as_tibble()
}

# Canonicalize fiber concentrations to {0,100,1000,10000} as character labels
canonicalize_concentration <- function(x) {
  v <- suppressWarnings(as.numeric(x))              # tolerate "100.0", factors, etc.
  allowed <- c(0, 100, 1000, 10000)
  # Snap to nearest allowed within a 1% (or 1 unit) tolerance
  snapped <- vapply(v, function(y) {
    if (!is.finite(y)) return(NA_real_)
    idx <- which.min(abs(allowed - y))
    tol <- max(1, 0.01 * allowed[idx])
    if (abs(allowed[idx] - y) <= tol) allowed[idx] else y
  }, numeric(1))
  as.character(as.integer(round(snapped)))          # canonical string labels
}
