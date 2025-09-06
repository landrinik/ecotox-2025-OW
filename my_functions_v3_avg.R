# Apply QC flags to a results dataframe using flag_plate_qc function
apply_qc_flags <- function(results_df, assay_type = NULL) 
{
  if (!"calculated_concentration" %in% colnames(results_df)) {
    stop("QC error: 'calculated_concentration' column missing in input data.")
  }
  if (nrow(results_df) == 0) return(results_df)
  
  flag_thresholds <- qc_thresholds_for_assay(assay_type)
  
  qc_list <- flag_plate_qc(results_df, flag_thresholds)
  
  if (!"qc_flag" %in% colnames(results_df)) {
    results_df$qc_flag <- FALSE
  }
  
  flagged_wells <- unique(c(
    if (!is.null(qc_list$out_of_bounds)) qc_list$out_of_bounds$well_name else character(0),
    if (!is.null(qc_list$negative)) qc_list$negative$well_name else character(0),
    if (!is.null(qc_list$na)) qc_list$na$well_name else character(0)
  ))
  
  results_df$qc_flag[results_df$well_name %in% flagged_wells] <- TRUE
  
  return(results_df)
}

# Compute replicate group statistics and add CV and QC flags
compute_qc_stats <- function(df, concentration_col = "calculated_concentration", group_cols) 
{
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      Mean = mean(.data[[concentration_col]], na.rm = TRUE),
      SD = sd(.data[[concentration_col]], na.rm = TRUE),
      N = sum(!is.na(.data[[concentration_col]])),
      CV = ifelse(Mean == 0 | is.na(Mean), NA_real_, 100 * SD / Mean),
      .groups = "drop"
    )
}

# MASTER ANALYSIS FUNCTION ------------------------------------------
analyze_assay_data <- function(full_dataset, 
                               assay_to_filter, 
                               concentrations_vector,
                               fiber_to_filter,
                               week_to_filter,
                               sample_type_to_filter,
                               replicate_to_filter,
                               global_standards,
                               tissue_weights) 
{
  base_data <- full_dataset %>%
    filter(
      assay_type == assay_to_filter,
      fiber_type == fiber_to_filter,
      sample_week == week_to_filter,
      sample_type == sample_type_to_filter
    )
  
  dilution_factors <- c(
    "SOD" = 50,
    "CAT" = 10 
  )
  
  if (assay_to_filter %in% c("SOD", "CAT", "ALP", "AChE", "Bradford")) {
    plate_data <- base_data %>%
      filter(plate_replicate == replicate_to_filter)
    
    if (nrow(plate_data) == 0) { return(NULL) }
    
    if (assay_to_filter %in% c("SOD", "CAT")) {
      results <- process_sod_cat(plate_data, 
                                 concentrations_vector, 
                                 assay_to_filter,
                                 global_standards,
                                 dilution_factor = dilution_factors[assay_to_filter])
    } else if (assay_to_filter == "Bradford") { 
      results <- process_bradford(plate_data, 
                                  concentrations_vector, 
                                  assay_to_filter,
                                  global_standards)
    } else if (assay_to_filter == "ALP") {
      results <- process_alp(plate_data, 
                             concentrations_vector,
                             global_standards,
                             tissue_weights)
    } else { 
      results <- process_ache(plate_data,
                              concentrations_vector,
                              global_standards)
    }
  } else if (assay_to_filter == "ACP") {
    plate_data <- base_data %>%
      filter(plate_replicate %in% c(replicate_to_filter, "BlankRep"))
    
    if (length(unique(plate_data$plate_replicate)) < 2) { return(NULL) }
    
    results <- process_acp(plate_data, 
                           concentrations_vector,
                           tissue_weights)
  } else {
    stop("Unknown assay type provided!")
  }
  
  # Calculate QC stats and attach flags
  group_columns <- c("fiber_type", "assay_type", "sample_week", "sample_type", "plate_replicate", "well_name")
  
  # Use the correct concentration column name in final data (calculated_concentration or normalized_activity)
  concentration_col <- ifelse("calculated_concentration" %in% names(results$final_data),
                              "calculated_concentration", "normalized_activity")
  
  qc_stats <- compute_qc_stats(results$final_data, concentration_col, group_columns)
  
  results$final_data <- results$final_data %>%
    left_join(qc_stats, by = group_columns)
  
  return(results)
}

# HELPER FUNCTION: SOD AND CAT processing -----------------------------------
process_sod_cat <- function(plate_data, 
                            concentrations_vector,
                            assay_to_filter,
                            global_standards,
                            dilution_factor = 1) 
{
  n_standards <- length(concentrations_vector)
  blank_well_letter <- LETTERS[n_standards]
  
  current_fiber <- unique(plate_data$fiber_type)
  current_sample <- unique(plate_data$sample_type)
  
  global_blank_abs <- global_standards %>%
    filter(
      assay_type == assay_to_filter,
      fiber_type == current_fiber,
      sample_type == current_sample,
      wells == blank_well_letter
    ) %>%
    pull(avg_value)
  
  blank_crtd_data <- plate_data %>%
    mutate(across(x1:x12, ~ (global_blank_abs - .x) / global_blank_abs * 100))
  
  wells_to_get <- LETTERS[1:n_standards]
  
  long_standards <- pivot_standards(blank_crtd_data, wells_to_get)
  renamed_long_standards <- conditionally_rename_absorbance(long_standards)
  std_curve_data <- renamed_long_standards %>%
    mutate(known_conc = concentrations_vector)
  
  pl_model <- drm(percent_inhibition ~ known_conc,
                  data = std_curve_data,
                  fct = LL.4(names = c("Slope", "Lower", "Upper", "ED50")))
  
  long_samples <- pivot_samples(blank_crtd_data)
  renamed_long_samples <- conditionally_rename_absorbance(long_samples)
  
  final_results <- renamed_long_samples %>%
    mutate(calculated_concentration = suppressWarnings(
      ED(pl_model, percent_inhibition, type = "absolute")[, 1]
    )) %>%
    mutate(calculated_concentration = if_else(is.nan(calculated_concentration), 0, calculated_concentration)) %>%
    mutate(calculated_concentration = calculated_concentration * dilution_factor)
  
  final_results <- apply_qc_flags(final_results, assay_type = assay_to_filter)
  
  # CREATE AND DISPLAY OUTPUTS
  
  # Extract the model's coefficient table
  model_coefs <- coef(summary(pl_model))
  
  # Get the four parameters from the "Estimate" column
  b_slope <- model_coefs["Slope:(Intercept)", "Estimate"]
  c_lower <- model_coefs["Lower:(Intercept)", "Estimate"]
  d_upper <- model_coefs["Upper:(Intercept)", "Estimate"]
  e_ed50 <- model_coefs["ED50:(Intercept)", "Estimate"]
  
  # Format the parameters into a clean text label
  # The "\n" creates a new line
  model_label <- paste(
    paste0("Slope: ", round(b_slope, 2)),
    paste0("Upper: ", round(d_upper, 2)),
    paste0("Lower: ", round(c_lower, 2)),
    paste0("ED50: ", round(e_ed50, 2)),
    sep = "\n"
  )
  # Get the unique identifiers for this specific plate to use in the title
  current_fiber <- unique(plate_data$fiber_type)
  current_week <- unique(plate_data$sample_week)
  current_sample <- unique(plate_data$sample_type)
  current_rep <- unique(plate_data$plate_replicate[plate_data$plate_replicate != "BlankRep"])
  
  # Create a dynamic title string
  dynamic_title <- paste(assay_to_filter," Standard Curve:", current_fiber, current_sample, current_week, current_rep)
  
  # --- Prepare data for plotting (handle the zero for the log scale) ---
  plot_data <- std_curve_data %>%
    mutate(known_conc_plot = if_else(known_conc == 0, 0.01, 
                                     known_conc))
  
  # --- Generate the 4PL curve data from the model ---
  x_curve <- seq(min(plot_data$known_conc_plot), max(plot_data$known_conc), length.out = 100)
  y_curve <- predict(pl_model, newdata = data.frame(known_conc = x_curve))
  curve_data <- data.frame(known_conc_plot = x_curve, absorbance = y_curve)
  
  print(names(std_curve_data))
  # --- Create the publication-ready plot ---
  publication_plot <- ggplot(data = plot_data, 
                             aes(x = known_conc_plot, 
                                 y = percent_inhibition)) +
    geom_point(size = 2) +
    geom_line(data = curve_data, 
              aes(x = known_conc_plot, y = absorbance),
              color = "black", linewidth = 1) +
    annotate(
      geom = "label",
      x = min(plot_data$known_conc_plot), 
      y = max(plot_data$percent_inhibition), 
      hjust = 0,  
      vjust = 1,
      label = model_label,
      fill = "white"
    ) +
    coord_cartesian(ylim = c(NA, max(plot_data$percent_inhibition) * 1.1)) +
    scale_x_log10(labels = scales::number_format(accuracy = 1, big.mark = ",")) +
    theme_classic() +
    labs(
      title = dynamic_title,
      x = "[Standards] (U/mL)",
      y = expression(paste("Corrected %Inhibition", {}[560*nm]))
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  model_summary <- summary(pl_model)
  print(publication_plot)
  print(model_summary)
  list(
    final_data = final_results,
    plot = publication_plot,
    model = pl_model,
    summary = model_summary
  )
}

# HELPER FUNCTION: Bradford processing -----------------------------------
process_bradford <- function(plate_data, 
                             concentrations_vector,
                             assay_to_filter,
                             global_standards) 
{
  n_standards <- length(concentrations_vector)
  blank_well_letter <- LETTERS[n_standards]
  current_fiber <- unique(plate_data$fiber_type)
  current_sample <- unique(plate_data$sample_type)
  
  dilution_factor_bf = 20
  
  global_blank_abs <- global_standards %>%
    filter(
      assay_type == assay_to_filter,
      fiber_type == current_fiber,
      sample_type == current_sample,
      wells == blank_well_letter
    ) %>%
    pull(avg_value)
  
  blank_crtd_data <- plate_data %>%
    mutate(across(x1:x12, ~ .x - global_blank_abs))
  
  wells_to_get <- LETTERS[1:n_standards]
  
  long_standards <- pivot_standards(blank_crtd_data, wells_to_get)
  std_curve_data <- long_standards %>%
    mutate(known_conc = concentrations_vector)
  
  pl_model <- drm(absorbance ~ known_conc,
                  data = std_curve_data,
                  fct = LL.4(names = c("Slope", "Lower", "Upper", "ED50")))
  
  long_samples <- pivot_samples(blank_crtd_data)
  
  final_results <- long_samples %>%
    mutate(calculated_concentration = suppressWarnings(
      ED(pl_model, absorbance, type = "absolute")[, 1]
    )) %>%
    mutate(calculated_concentration = if_else(is.nan(calculated_concentration), 0, calculated_concentration)) %>%
    mutate(calculated_concentration = calculated_concentration * dilution_factor_bf)
  
  final_results <- apply_qc_flags(final_results, assay_type = assay_to_filter)
  
  # CREATE AND DISPLAY OUTPUTS
  
  # Extract the model's coefficient table
  model_coefs <- coef(summary(pl_model))
  
  # Get the four parameters from the "Estimate" column
  b_slope <- model_coefs["Slope:(Intercept)", "Estimate"]
  c_lower <- model_coefs["Lower:(Intercept)", "Estimate"]
  d_upper <- model_coefs["Upper:(Intercept)", "Estimate"]
  e_ed50 <- model_coefs["ED50:(Intercept)", "Estimate"]
  
  # Format the parameters into a clean text label
  # The "\n" creates a new line
  model_label <- paste(
    paste0("Slope: ", round(b_slope, 2)),
    paste0("Upper: ", round(d_upper, 2)),
    paste0("Lower: ", round(c_lower, 2)),
    paste0("ED50: ", round(e_ed50, 2)),
    sep = "\n"
  )
  
  # --- Prepare data for plotting (handle the zero for the log scale) ---
  plot_data <- std_curve_data %>%
    mutate(known_conc_plot = if_else(known_conc == 0, 0.01, 
                                     known_conc))
  
  # --- Generate the 4PL curve data from the model ---
  x_curve <- seq(min(plot_data$known_conc_plot), max(plot_data$known_conc), length.out = 100)
  y_curve <- predict(pl_model, newdata = data.frame(known_conc = x_curve))
  curve_data <- data.frame(known_conc_plot = x_curve, absorbance = y_curve)
  
  # --- Create the publication-ready plot ---
  publication_plot <- ggplot(data = plot_data, 
                             aes(x = known_conc_plot, 
                                 y = absorbance)) +
    geom_point(size = 2) +
    geom_line(data = curve_data, 
              color = "black", linewidth = 1) +
    annotate(
      geom = "label",
      x = min(plot_data$known_conc_plot), 
      y = max(plot_data$absorbance), 
      hjust = 0,  
      vjust = 1,
      label = model_label,
      fill = "white"
    ) +
    coord_cartesian(ylim = c(NA, max(plot_data$absorbance) * 1.1)) +
    scale_x_log10() +
    theme_classic() +
    labs(
      title = paste(assay_to_filter, "Standard Curve"),
      x = "[Standards] (mg/mL)",
      y = expression(paste("Corrected Abs", {}[595*nm]))
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  model_summary <- summary(pl_model)
  print(publication_plot)
  print(model_summary)
  list(
    final_data = final_results,
    plot = publication_plot,
    model = pl_model
  )
}

# HELPER FUNCTION: ACP processing -----------------------------------
process_acp <- function(plate_data, 
                        concentrations_vector,
                        tissue_weights) 
{
  # Split into blanks and samples
  blanks_wide <- plate_data %>% filter(plate_replicate == "BlankRep")
  samples_wide <- plate_data %>% filter(plate_replicate != "BlankRep")
  
  # Prepare blank absorbance values in long format
  blanks_long <- pivot_samples(blanks_wide) %>%
    rename(blank_value = absorbance) %>%
    dplyr::select(sample_week, well_name, blank_value)
  
  # Prepare sample absorbance values in long format
  samples_long <- pivot_samples(samples_wide)
  
  # Subtract blanks from samples
  corrected_data <- samples_long %>%
    left_join(blanks_long, by = c("sample_week", "well_name")) %>%
    mutate(corrected_absorbance = absorbance - blank_value)
  
  # Prepare standards for calibration curve
  standards_wide <- plate_data %>%
    filter(plate_replicate != "BlankRep") %>%
    dplyr::select(wells, x12)
  
  # Extract blank wells
  blank_d <- standards_wide %>% filter(wells == "D") %>% pull(x12)
  blank_h <- standards_wide %>% filter(wells == "H") %>% pull(x12)
  
  # Average blank absorbance for zero standard
  zero_std_abs <- mean(c(blank_d, blank_h))
  
  # Generate corrected standard absorbance values
  acp_std_curve <- standards_wide %>%
    mutate(corrected_abs = case_when(
      wells %in% c("A", "B", "C") ~ x12 - blank_d,
      wells %in% c("E", "F", "G") ~ x12 - blank_h,
      TRUE ~ 0
    )) %>%
    filter(wells %in% c("A", "B", "C", "E", "F", "G")) %>%
    mutate(known_conc = rep(concentrations_vector[1:3], 2))
  
  # Fit linear model to standard curve data
  lm_model <- lm(corrected_abs ~ known_conc, data = acp_std_curve)
  slope <- coef(lm_model)[2]
  
  # Dilution factor
  dilution_factor <- 3.75
  
  # Calculate initial concentrations for samples
  initial_data <- corrected_data %>%
    mutate(calculated_concentration = ((corrected_absorbance - zero_std_abs) / (30 * slope)) * dilution_factor)
  
  # Normalize by tissue weight
  final_data <- initial_data %>%
    left_join(tissue_weights, by = c("assay_type", "fiber_type", "sample_week", "sample_type", "well_name")) %>%
    mutate(tissue_modifier = 80 / tissue_weight_mg,
           calculated_concentration = calculated_concentration * tissue_modifier)
  
  # Calculate replicate statistics
  group_cols <- c("fiber_type", "assay_type", "sample_week", "sample_type", "well_name")
  replicate_stats <- final_data %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(Mean = mean(calculated_concentration, na.rm = TRUE),
              SD = sd(calculated_concentration, na.rm = TRUE),
              N = n(),
              CV = ifelse(Mean == 0, NA_real_, 100 * SD / Mean),
              .groups = "drop")
  
  # Join stats back into data
  final_data <- final_data %>%
    left_join(replicate_stats, by = group_cols)
  
  # Apply QC flags with designated thresholds for ACP assay
  final_data <- apply_qc_flags(final_data, assay_type = "ACP")
  
  # Get the R-squared value from the model summary
  r_squared <- summary(lm_model)$r.squared
  
  # Format the value into a clean text label
  # \u00B2 is ^2
  r_squared_label <- paste0("R\u00B2 = ", round(r_squared, 3))
  
  # Get the unique identifiers for this specific plate to use in the title
  current_fiber <- unique(plate_data$fiber_type)
  current_week <- unique(plate_data$sample_week)
  current_sample <- unique(plate_data$sample_type)
  current_rep <- unique(plate_data$plate_replicate[plate_data$plate_replicate != "BlankRep"])
  
  # Create a dynamic title string
  dynamic_title <- paste("ACP Standard Curve:", current_fiber, current_sample, current_week, current_rep)
  
  # Create plot
  publication_plot <- ggplot(data = acp_std_curve, aes(x = known_conc, y = corrected_abs)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
    annotate(
      geom = "text",
      x = Inf,    
      y = -Inf,       
      hjust = 1.1,    
      vjust = -0.5,     
      label = r_squared_label,
      size = 4
    ) +
    theme_classic() +
    labs(
      title = dynamic_title,
      x = "[Nitrophenol] (\u03BCM)",
      y = expression(paste("Corrected OD", {}[405*nm]))
    ) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12)
    )
  # Print plot and model summary
  print(publication_plot)
  print(summary(lm_model))
  
  # Return results
  list(
    final_data = final_data,
    plot = publication_plot,
    model = lm_model
  )
}

# HELPER FUNCTION: ALP processing -----------------------------------
process_alp <- function(plate_data, 
                        concentrations_vector = NULL,
                        global_standards,
                        tissue_weights) 
{
  current_assay <- unique(plate_data$assay_type)
  current_fiber <- unique(plate_data$fiber_type)
  current_week <- unique(plate_data$sample_week)
  current_sample_type <- unique(plate_data$sample_type)
  
  summary_values <- global_standards %>%
    filter(assay_type == current_assay,
           fiber_type == current_fiber,
           sample_week == current_week,
           sample_type == current_sample_type) %>%
    group_by(reading) %>%
    summarise(
      mean_calibrator = mean(avg_value[wells %in% c("A","B")]),
      mean_blank = mean(avg_value[wells %in% c("C","D")])
    )
  
  blank_corrected_data <- plate_data %>%
    left_join(summary_values, by = "reading") %>%
    mutate(across(x1:x12, ~ .x - mean_blank))
  
  long_corrected_data <- pivot_samples(blank_corrected_data)
  
  delta_abs_data <- long_corrected_data %>%
    pivot_wider(
      id_cols = c(well_name, fiber_type, assay_type, sample_week, sample_type, plate_replicate),
      names_from = reading,
      values_from = absorbance) %>%
    mutate(Delta_Abs = SecondRead - FirstRead)
  
  dilution_factor <- 3.75
  final_alp_data <- delta_abs_data %>%
    mutate(Calculated_Activity = (Delta_Abs * 35.32 / 
                                    (summary_values$mean_calibrator[1] - summary_values$mean_blank[1])) * dilution_factor) %>%
    filter(!well_name %in% c("A12", "B12", "C12", "D12"))
  
  final_normalized_data <- final_alp_data %>%
    left_join(tissue_weights, by = c("assay_type", "fiber_type", "sample_week", "sample_type", "well_name")) %>%
    mutate(
      tissue_modifier = 80 / tissue_weight_mg,
      normalized_activity = Calculated_Activity * tissue_modifier
    ) %>%
    dplyr::select(fiber_type, assay_type, sample_week, sample_type, plate_replicate, well_name, normalized_activity) %>%
    rename(calculated_concentration = normalized_activity)
  
  final_results <- apply_qc_flags(final_normalized_data, assay_type = "ALP")
  
  list(final_data = final_results, plot = NULL, model = NULL)
}

# HELPER FUNCTION: AChE processing -----------------------------------
process_ache <- function(plate_data,
                         concentrations_vector = NULL,
                         global_standards) 
{
  dilution_factor = 12.5
  current_fiber <- unique(plate_data$fiber_type)
  current_assay <- unique(plate_data$assay_type)
  
  if (current_fiber == "Cotton") {
    summary_values <- global_standards %>%
      filter(assay_type == current_assay, fiber_type == current_fiber) %>%
      group_by(reading) %>%
      summarise(
        mean_blank = mean(avg_value[wells %in% c("A", "B", "C", "D")]),
        mean_calibrator = mean(avg_value[wells %in% c("E", "F", "G", "H")])
      )
  } else if (current_fiber == "PET") {
    summary_values <- global_standards %>%
      filter(assay_type == current_assay, fiber_type == current_fiber) %>%
      group_by(reading) %>%
      summarise(
        mean_blank = mean(avg_value[wells %in% c("E", "F", "G", "H")]),
        mean_calibrator = mean(avg_value[wells %in% c("A", "B", "C", "D")])
      )
  } else {
    stop("Unknown fiber_type found in AChE assay!")
  }
  
  blank_corrected_data <- plate_data %>%
    left_join(summary_values, by = "reading") %>%
    mutate(across(x1:x11, ~ .x - mean_blank))
  
  final_ache_data <- blank_corrected_data %>%
    pivot_samples() %>%
    pivot_wider(
      id_cols = c(well_name, fiber_type, assay_type, sample_week, sample_type, plate_replicate),
      names_from = reading,
      values_from = absorbance) %>%
    mutate(Calculated_Activity = (((SecondRead - FirstRead) * dilution_factor * 200) / 
             (summary_values$mean_calibrator[2] - summary_values$mean_blank[2]))/1000) %>%
    rename(calculated_concentration = Calculated_Activity)
  
  final_results <- apply_qc_flags(final_ache_data, assay_type = "AChE")
  
  list(final_data = final_results, plot = NULL, model = NULL)
}

# HELPER FUNCTIONS: CONDITIONALLY RENAME AND PIVOT -----------------------------------

conditionally_rename_absorbance <- function(data) 
{
  assay_type <- unique(data$assay_type)
  if (assay_type %in% c("SOD", "CAT")) {
    renamed_data <- data %>%
      rename(percent_inhibition = absorbance)
    return(renamed_data)
  } else {
    return(data)
  }
}

pivot_standards <- function(data_to_pivot, wells_to_filter) 
{
  data_to_pivot %>%
    filter(wells %in% wells_to_filter) %>%
    dplyr::select(fiber_type, assay_type, sample_week, sample_type, plate_replicate, reading, wells, x12) %>%
    pivot_longer(cols = x12, names_to = "well_column", values_to = "absorbance") %>%
    mutate(well_column = stringr::str_remove(well_column, "x")) %>%
    unite(col = "well_name", wells, well_column, sep = "")
}

pivot_samples <- function(data_to_pivot) 
{
  data_to_pivot %>%
    pivot_longer(cols = x1:x11, names_to = "well_column", values_to = "absorbance") %>%
    mutate(well_column = stringr::str_remove(well_column, "x")) %>%
    unite(col = "well_name", wells, well_column, sep = "") %>%
    dplyr::select(-x12) %>%
    filter(!well_name %in% c("E11", "F11", "G11", "H11"))
}
