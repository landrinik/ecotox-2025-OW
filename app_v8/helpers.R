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