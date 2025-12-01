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

# =============================================================================
# Create 14-level combined treatment variable (UNIFIED VERSION)
# =============================================================================
#' Combines fiber_type, chem_treatment, and concentration into one categorical variable
#' to avoid collinearity issues with control groups
#' 
#' @param data Data frame containing fiber_type, chem_treatment, concentration, is_control
#' @param reference_level Which level should be the reference (default: "Control Cotton")
#' @return Data frame with new 'treatment_combined' variable

create_combined_treatment <- function(data, reference_level = "Control Cotton") {
  df <- data
  
  # Map legacy snake_case if needed, without clobbering existing columns
  if (!"fiber_type" %in% names(df) && "fibertype" %in% names(df)) {
    df <- dplyr::rename(df, fiber_type = fibertype)
  }
  if (!"chem_treatment" %in% names(df) && "chemtreatment" %in% names(df)) {
    df <- dplyr::rename(df, chem_treatment = chemtreatment)
  }
  if (!"is_control" %in% names(df) && "iscontrol" %in% names(df)) {
    df <- dplyr::rename(df, is_control = iscontrol)
  }
  
  # Ensure a discrete 'concentration' column exists (snap to {0,100,1000,10000})
  if (!"concentration" %in% names(df)) {
    # Only coalesce from columns that actually exist and have length > 0
    candidates <- list()
    
    if ("fiber_concentration_numeric" %in% names(df)) {
      candidates[[length(candidates) + 1]] <- suppressWarnings(as.numeric(df[["fiber_concentration_numeric"]]))
    }
    if ("fiber_concentration" %in% names(df)) {
      candidates[[length(candidates) + 1]] <- suppressWarnings(as.numeric(df[["fiber_concentration"]]))
    }
    # Only add dose columns if they exist
    if ("dose_log10" %in% names(df)) {
      candidates[[length(candidates) + 1]] <- suppressWarnings(as.numeric(df[["dose_log10"]]))
    }
    if ("doselog10" %in% names(df)) {
      candidates[[length(candidates) + 1]] <- suppressWarnings(as.numeric(df[["doselog10"]]))
    }
    
    # Coalesce only if we have at least one candidate of correct length
    if (length(candidates) == 0 || all(vapply(candidates, length, integer(1)) == 0)) {
      # Fallback: create 0-vector of correct length
      conc_num <- rep(0, nrow(df))
    } else {
      # Filter out zero-length vectors before coalescing
      candidates <- Filter(function(x) length(x) == nrow(df), candidates)
      conc_num <- do.call(dplyr::coalesce, candidates)
      if (length(conc_num) == 0 || all(is.na(conc_num))) {
        conc_num <- rep(0, nrow(df))
      }
    }
    
    df[["concentration"]] <- canonicalize_concentration(conc_num)
  }
  
  # *** FIXED: Build with consistent casing - lowercase matching ***
  # Standardize fiber_type and chem_treatment to lowercase for robust matching
  df <- df %>%
    dplyr::mutate(
      fiber_type_lower = tolower(as.character(fiber_type)),
      chem_treatment_lower = tolower(as.character(chem_treatment)),
      concentration_char = as.character(concentration)
    )
  
  # Build the combined treatment variable with explicit case_when
  df <- df %>%
    dplyr::mutate(
      treatment_combined = dplyr::case_when(
        # Control groups - use "Control Cotton" and "Control Pet" (not PET)
        is_control == 1 & fiber_type_lower == "cotton" ~ "Control Cotton",
        is_control == 1 & fiber_type_lower == "pet" ~ "Control Pet",
        
        # Cotton treatments
        fiber_type_lower == "cotton" & chem_treatment_lower == "untreated" & concentration_char == "100" ~ "Untreated Cotton 100 mfL",
        fiber_type_lower == "cotton" & chem_treatment_lower == "untreated" & concentration_char == "1000" ~ "Untreated Cotton 1000 mfL",
        fiber_type_lower == "cotton" & chem_treatment_lower == "untreated" & concentration_char == "10000" ~ "Untreated Cotton 10000 mfL",
        fiber_type_lower == "cotton" & chem_treatment_lower == "treated" & concentration_char == "100" ~ "Treated Cotton 100 mfL",
        fiber_type_lower == "cotton" & chem_treatment_lower == "treated" & concentration_char == "1000" ~ "Treated Cotton 1000 mfL",
        fiber_type_lower == "cotton" & chem_treatment_lower == "treated" & concentration_char == "10000" ~ "Treated Cotton 10000 mfL",
        
        # PET treatments - use "Pet" (not "PET") to match str_to_title behavior
        fiber_type_lower == "pet" & chem_treatment_lower == "untreated" & concentration_char == "100" ~ "Untreated Pet 100 mfL",
        fiber_type_lower == "pet" & chem_treatment_lower == "untreated" & concentration_char == "1000" ~ "Untreated Pet 1000 mfL",
        fiber_type_lower == "pet" & chem_treatment_lower == "untreated" & concentration_char == "10000" ~ "Untreated Pet 10000 mfL",
        fiber_type_lower == "pet" & chem_treatment_lower == "treated" & concentration_char == "100" ~ "Treated Pet 100 mfL",
        fiber_type_lower == "pet" & chem_treatment_lower == "treated" & concentration_char == "1000" ~ "Treated Pet 1000 mfL",
        fiber_type_lower == "pet" & chem_treatment_lower == "treated" & concentration_char == "10000" ~ "Treated Pet 10000 mfL",
        
        TRUE ~ "Other"
      )
    ) %>%
    # Remove temporary columns
    dplyr::select(-fiber_type_lower, -chem_treatment_lower, -concentration_char)
  
  # Convert to factor with specified reference level first
  if (!is.null(reference_level) && reference_level %in% unique(df$treatment_combined)) {
    all_levels <- c(
      reference_level,
      setdiff(unique(df$treatment_combined), c(reference_level, "Other"))
    )
    df$treatment_combined <- factor(df$treatment_combined, levels = all_levels)
  } else {
    # If reference not found, just make it a factor
    df$treatment_combined <- factor(df$treatment_combined)
  }
  
  # Print summary
  cat("\n=== Combined Treatment Variable Created ===\n")
  cat("Reference level:", reference_level, "\n\n")
  cat("Level counts:\n")
  print(table(df$treatment_combined))
  cat("\n")
  
  return(df)
}



# ============================================================================
# EXPORT UTILITIES: high-resolution plots & tables
# ============================================================================

# ---- Plot styling helpers ----
#' Consistent report theme for ggplot objects
theme_report <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      plot.title = ggplot2::element_text(size = base_size + 2, face = "bold"),
      axis.title = ggplot2::element_text(size = base_size, face = "bold")
    )
}

#' Apply user overrides for title, axis labels, legend, palette
#' @param p ggplot object
#' @param title plot title (if non-empty)
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param legend_pos one of "right", "left", "top", "bottom", "none"
#' @param palette one of "brewer_set1", "viridis", "greys", "custom"
#' @return modified ggplot object
apply_plot_overrides <- function(p, title = NULL, xlab = NULL, ylab = NULL,
                                 legend_pos = "right", palette = "brewer_set1") {
  # Apply labels
  if (!is.null(title) && nzchar(title)) p <- p + ggplot2::ggtitle(title)
  if (!is.null(xlab)  && nzchar(xlab))  p <- p + ggplot2::xlab(xlab)
  if (!is.null(ylab)  && nzchar(ylab))  p <- p + ggplot2::ylab(ylab)
  
  # Legend position
  if (!is.null(legend_pos) && legend_pos %in% c("right","left","top","bottom","none")) {
    p <- p + ggplot2::theme(legend.position = ifelse(legend_pos == "none", "none", legend_pos))
  }
  
  # Color/fill palettes
  if (palette == "brewer_set1") {
    p <- p + ggplot2::scale_color_brewer(palette = "Set1") + ggplot2::scale_fill_brewer(palette = "Set1")
  } else if (palette == "viridis") {
    p <- p + ggplot2::scale_color_viridis_d() + ggplot2::scale_fill_viridis_d()
  } else if (palette == "greys") {
    p <- p + ggplot2::scale_color_grey(start = 0.2, end = 0.6) + ggplot2::scale_fill_grey(start = 0.2, end = 0.6)
  }
  p
}

# ---- Plot saving helpers (300 dpi raster + editable vector) ----
#' Save a ggplot in all formats: PNG, TIFF, SVG, PDF
#' @param p ggplot object
#' @param file_base path without extension
#' @param width_in width in inches
#' @param height_in height in inches
#' @param dpi dots per inch for raster
save_plot_all_formats <- function(p, file_base, width_in = 7, height_in = 5, dpi = 300) {
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  
  # Ensure consistent theme
  p <- p + theme_report()
  
  # High-res PNG via ragg (or fallback to png if ragg unavailable)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = paste0(file_base, ".png"), width = width_in, height = height_in, 
                  units = "in", res = dpi, scaling = 1)
    print(p)
    grDevices::dev.off()
  } else {
    grDevices::png(filename = paste0(file_base, ".png"), width = width_in, height = height_in, 
                   units = "in", res = dpi)
    print(p)
    grDevices::dev.off()
  }
  
  # High-res TIFF
  grDevices::tiff(filename = paste0(file_base, ".tiff"), width = width_in, height = height_in, 
                  units = "in", res = dpi, compression = "lzw")
  print(p)
  grDevices::dev.off()
  
  # Editable SVG (requires svglite)
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file = paste0(file_base, ".svg"), width = width_in, height = height_in)
    print(p)
    grDevices::dev.off()
  }
  
  # Editable PDF
  grDevices::pdf(file = paste0(file_base, ".pdf"), width = width_in, height = height_in, 
                 useDingbats = FALSE)
  print(p)
  grDevices::dev.off()
  
  invisible(paste0(file_base, c(".png",".tiff",".svg",".pdf")))
}

# ---- Table (gt) export from a data frame ----
#' Build a gt table from a data frame with optional title
#' @param df data frame
#' @param title optional title for the table
#' @return gt object
build_gt_from_df <- function(df, title = NULL) {
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required for table exports. Install via: install.packages('gt')")
  }
  gt_tbl <- gt::gt(df)
  if (!is.null(title) && nzchar(title)) {
    gt_tbl <- gt::tab_header(gt_tbl, title = gt::md(title))
  }
  gt_tbl
}

#' Save a gt table as PNG, PDF, HTML
#' @param gt_tbl gt object
#' @param file_path_no_ext path without extension
#' @param width_px width in pixels for PNG
#' @param height_px height in pixels for PNG
#' @param dpi DPI for PNG rendering
save_gt_table <- function(gt_tbl, file_path_no_ext, width_px = 2000, height_px = 1400, dpi = 300) {
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required. Install via: install.packages('gt')")
  }
  
  # PNG (high DPI)
  tryCatch({
    gt::gtsave(
      data = gt_tbl,
      filename = paste0(file_path_no_ext, ".png"),
      vwidth = width_px, vheight = height_px, dpi = dpi
    )
  }, error = function(e) warning("PNG export failed: ", e$message))
  
  # PDF (vector editable)
  tryCatch({
    gt::gtsave(
      data = gt_tbl,
      filename = paste0(file_path_no_ext, ".pdf")
    )
  }, error = function(e) warning("PDF export failed: ", e$message))
  
  # HTML
  tryCatch({
    gt::gtsave(
      data = gt_tbl,
      filename = paste0(file_path_no_ext, ".html")
    )
  }, error = function(e) warning("HTML export failed: ", e$message))
  
  invisible(file_path_no_ext)
}


# ==============================================================================
# GRANULAR ANOVA SUMMARY GENERATOR (ROBUST & REFIT-AWARE)
# ==============================================================================
generate_granular_anova <- function(model, data, return_format = "numeric") {
  require(emmeans)
  require(dplyr)
  require(tidyr)
  require(stringr)
  
  # 1. Identify predictors actually used in the model formula
  model_terms <- tryCatch(
    all.vars(stats::formula(model)), 
    error = function(e) character(0)
  )
  
  # 2. Detect Time Variable (Factor 'week' vs Spline 'week_numeric')
  time_var <- if ("week_numeric" %in% model_terms) "week_numeric" else "week"
  
  # 3. Define Stratification
  potential_by <- c(time_var, "fiber_type")
  by_vars <- intersect(potential_by, model_terms)
  
  tryCatch({
    # 4. Prepare 'at' list for Splines
    at_list <- list()
    if (time_var == "week_numeric" && "week_numeric" %in% names(data)) {
      at_list$week_numeric <- sort(unique(data$week_numeric))
    }
    
    # 5. Run Joint Tests
    jt <- emmeans::joint_tests(
      model, 
      by = by_vars, 
      data = data,
      at = if(length(at_list) > 0) at_list else NULL
    ) %>% as.data.frame()
    
    if (nrow(jt) == 0) return(data.frame(Message = "No joint tests results available."))
    
    # 6. Standardize 'week' column for downstream processing
    if (time_var == "week_numeric" && "week_numeric" %in% names(jt)) {
      jt$week <- jt$week_numeric
    }
    
    # 7. Add Context (Fiber/Sample) if missing
    if (!"fiber_type" %in% names(jt)) {
      val <- if("fiber_type" %in% names(data)) unique(as.character(data$fiber_type)) else "Mixed"
      jt$fiber_type <- if(length(val) == 1) val else "Mixed"
    }
    
    if (!"sample_type" %in% names(jt)) {
      val <- "Unknown"
      if ("sample_type" %in% names(data)) val <- unique(as.character(data$sample_type))
      else if ("tissue_type" %in% names(data)) val <- unique(as.character(data$tissue_type))
      jt$sample_type <- if(length(val) == 1) val else "Mixed"
    }
    
    # 8. Format Terms & Values
    jt <- jt %>%
      dplyr::mutate(
        term_clean = stringr::str_to_title(gsub("_", " ", `model term`))
      )
    
    if ("F.ratio" %in% names(jt)) {
      jt$stat_val <- jt$F.ratio; jt$stat_lab <- "F"
    } else if ("Chisq" %in% names(jt)) {
      jt$stat_val <- jt$Chisq; jt$stat_lab <- "Chi2"
    } else {
      jt$stat_val <- NA_real_; jt$stat_lab <- "?"
    }
    
    if (return_format == "composite") {
      jt <- jt %>%
        dplyr::mutate(
          cell_value = dplyr::case_when(
            is.na(p.value) | is.na(stat_val) ~ "NA",
            TRUE ~ sprintf("%s (%s=%.1f)", format.pval(p.value, digits = 3, eps = 0.001), stat_lab, stat_val)
          )
        )
    } else {
      jt <- jt %>% dplyr::mutate(cell_value = p.value)
    }
    
    # 9. Pivot
    if ("week" %in% names(jt)) {
      jt$week_label <- paste("Week", jt$week)
      jt_wide <- jt %>%
        dplyr::select(sample_type, fiber_type, term_clean, week_label, cell_value) %>%
        tidyr::pivot_wider(names_from = week_label, values_from = cell_value) %>%
        dplyr::arrange(sample_type, fiber_type, term_clean) %>%
        dplyr::rename(`Sample` = sample_type, `Fiber` = fiber_type, `Effect Tested` = term_clean)
      return(jt_wide)
    } else {
      return(jt %>% dplyr::select(sample_type, fiber_type, term_clean, stat_lab, stat_val, p.value))
    }
    
  }, error = function(e) {
    return(data.frame(Error = paste("Granular ANOVA Error:", e$message)))
  })
}

# ==============================================================================
# DOSE-RESPONSE VISUALIZATION (Refit-Aware + Shared Control Point)
# ==============================================================================

generate_dose_response_data <- function(model, data) {
  require(emmeans)
  require(dplyr)
  
  # Helper: Fix CI names
  fix_ci_names <- function(d) {
    if ("asymp.LCL" %in% names(d)) d <- dplyr::rename(d, lower.CL = asymp.LCL)
    if ("asymp.UCL" %in% names(d)) d <- dplyr::rename(d, upper.CL = asymp.UCL)
    d
  }
  
  vars <- all.vars(stats::formula(model))
  
  # Time variable check
  time_var <- if ("week_numeric" %in% vars) "week_numeric" else "week"
  
  # Determine Grouping Variables
  group_terms <- c("chem_treatment")
  if (time_var %in% vars) group_terms <- c(group_terms, time_var)
  if ("fiber_type" %in% vars) group_terms <- c(group_terms, "fiber_type")
  
  # --- SCENARIO A: Factor Dose ---
  if ("dose_factor" %in% vars) {
    f_str <- paste("~ dose_factor +", paste(group_terms, collapse = "+"))
    specs_formula <- stats::as.formula(f_str)
    
    emm <- emmeans::emmeans(model, specs = specs_formula, data = data)
    df <- as.data.frame(summary(emm))
    df <- fix_ci_names(df)
    
    if (time_var == "week_numeric" && "week_numeric" %in% names(df)) {
      df$week <- as.factor(df$week_numeric)
    }
    
    num_vals <- suppressWarnings(as.numeric(as.character(df$dose_factor)))
    df$dose_disp <- if(any(is.na(num_vals))) df$dose_factor else num_vals
    
    return(df)
  }
  
  # --- SCENARIO B: Numeric Log Dose (Continuous) ---
  
  # 1. Define Formulas
  f_base <- paste("~", paste(group_terms, collapse = "+"))
  specs_base <- stats::as.formula(f_base)
  f_dose <- paste("~", paste(c(group_terms, "dose_log10"), collapse = "+"))
  specs_dose <- stats::as.formula(f_dose)
  
  # 2. Control Points (Dose=0)
  at_ctrl <- list(dose_log10 = 0, is_control = 1)
  if (time_var == "week_numeric") {
    at_ctrl$week_numeric <- sort(unique(data$week_numeric))
  }
  
  emm_ctrl <- emmeans::emmeans(model, specs_base, at = at_ctrl, data = data)
  df_ctrl  <- as.data.frame(summary(emm_ctrl))
  df_ctrl  <- fix_ci_names(df_ctrl)
  df_ctrl$dose_disp <- 0 
  
  # --- VISUAL IMPROVEMENT: Force Shared Baseline at 0 ---
  # Overwrite the hypothetical "Treated @ 0" estimate with the "Untreated @ 0" estimate
  # so both lines visually originate from the Real Control data.
  if ("chem_treatment" %in% names(df_ctrl)) {
    # Identify Real Control rows (Untreated)
    is_untreated <- grepl("untreated", df_ctrl$chem_treatment, ignore.case = TRUE)
    
    if (any(is_untreated)) {
      real_controls <- df_ctrl[is_untreated, ]
      
      # Create a "Shadow" of these controls for the Treated group
      shadow_controls <- real_controls
      
      # Find the label for the Treated group (e.g. "Treated" or "treated")
      all_levels <- levels(df_ctrl$chem_treatment)
      if (is.null(all_levels)) all_levels <- unique(df_ctrl$chem_treatment)
      
      # Pick the level that is NOT untreated
      treated_lab <- all_levels[!grepl("untreated", all_levels, ignore.case = TRUE)][1]
      
      if (!is.na(treated_lab)) {
        shadow_controls$chem_treatment <- treated_lab
        # Replace original df_ctrl with the Shared Baseline version
        df_ctrl <- rbind(real_controls, shadow_controls)
      }
    }
  }
  # -----------------------------------------------------
  
  # 3. Treated Points
  at_dose <- list(dose_log10 = c(2, 3, 4), is_control = 0)
  if (time_var == "week_numeric") {
    at_dose$week_numeric <- sort(unique(data$week_numeric))
  }
  
  emm_dose <- emmeans::emmeans(model, specs_dose, at = at_dose, data = data)
  df_dose  <- as.data.frame(summary(emm_dose))
  df_dose  <- fix_ci_names(df_dose)
  df_dose$dose_disp <- 10^df_dose$dose_log10
  
  # 4. Standardize Time for Binding
  if (time_var == "week_numeric") {
    if ("week_numeric" %in% names(df_ctrl)) df_ctrl$week <- as.factor(df_ctrl$week_numeric)
    if ("week_numeric" %in% names(df_dose)) df_dose$week <- as.factor(df_dose$week_numeric)
  }
  
  # 5. Combine
  common_cols <- intersect(names(df_ctrl), names(df_dose))
  rbind(df_ctrl[, c(common_cols, "dose_disp"), drop = FALSE],
        df_dose[, c(common_cols, "dose_disp"), drop = FALSE])
}

# ==============================================================================
# SLOPE ANALYSIS (Refit-Aware)
# ==============================================================================
generate_dose_slopes <- function(model, data) {
  require(emmeans)
  
  vars <- all.vars(stats::formula(model))
  
  if (!"dose_log10" %in% vars) {
    return(list(error = "Model uses discrete 'Dose Factor'. Slopes cannot be calculated."))
  }
  
  # Dynamic Time Variable
  time_var <- if ("week_numeric" %in% vars) "week_numeric" else "week"
  
  # Construct Grid Variables
  grid_vars <- "chem_treatment"
  if ("fiber_type" %in% vars) grid_vars <- c(grid_vars, "fiber_type")
  
  # Formula construction
  lhs <- paste(grid_vars, collapse = "+")
  f_string <- if (time_var %in% vars) paste("~", lhs, "|", time_var) else paste("~", lhs)
  specs_formula <- stats::as.formula(f_string)
  
  # Prepare 'at' list
  at_list <- list()
  if (time_var == "week_numeric") {
    at_list$week_numeric <- sort(unique(data$week_numeric))
  }
  
  tryCatch({
    trends <- emmeans::emtrends(
      model, 
      specs = specs_formula, 
      var = "dose_log10",
      data = data,
      at = if(length(at_list) > 0) at_list else NULL
    )
    
    contrasts <- emmeans::contrast(trends, "pairwise", adjust = "none") 
    
    df_slopes <- as.data.frame(trends)
    df_contrasts <- as.data.frame(contrasts)
    
    # Standardize week
    if (time_var == "week_numeric") {
      if ("week_numeric" %in% names(df_slopes)) df_slopes$week <- as.factor(df_slopes$week_numeric)
      if ("week_numeric" %in% names(df_contrasts)) df_contrasts$week <- as.factor(df_contrasts$week_numeric)
    }
    
    list(slopes = df_slopes, comparisons = df_contrasts)
  }, error = function(e) {
    return(list(error = paste("Slope Calculation Error:", e$message)))
  })
}