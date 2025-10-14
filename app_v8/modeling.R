# Statistical Modeling Functions
# Optimized for count data and mixed effects

# ===== ENHANCED MODEL FITTING =====
# Enhanced modeling for count data (mf_counts) with automatic family selection
fit_endpoint_model <- function(df, endpoint_name, model_type = "lmer") {
  if (endpoint_name == "mf_counts" && model_type %in% c("glmer", "glm")) {
    # Test for overdispersion in count data
    mean_val <- mean(df$outcome, na.rm = TRUE)
    var_val <- var(df$outcome, na.rm = TRUE)
    dispersion_ratio <- ifelse(mean_val > 0, var_val / mean_val, 1)
    
    if (dispersion_ratio > 2) {
      # Use negative binomial for overdispersed counts
      if (model_type == "glm") {
        model <- MASS::glm.nb(
          outcome ~ treatment_category * week * log(fiber_concentration_numeric + 1),
          data = df
        )
      } else {
        model <- lme4::glmer.nb(
          outcome ~ treatment_category * week * log(fiber_concentration_numeric + 1) + 
            (1 | tank) + (1 | experiment_batch),
          data = df
        )
      }
      family_used <- "negative_binomial"
    } else {
      # Use Poisson for equidispersed counts
      if (model_type == "glm") {
        model <- glm(
          outcome ~ treatment_category * week * log(fiber_concentration_numeric + 1),
          family = poisson(), data = df
        )
      } else {
        model <- lme4::glmer(
          outcome ~ treatment_category * week * log(fiber_concentration_numeric + 1) + 
            (1 | tank) + (1 | experiment_batch),
          family = poisson(), data = df
        )
      }
      family_used <- "poisson"
    }
  } else {
    # Use linear models for continuous data
    if (model_type == "lm") {
      model <- lm(outcome ~ treatment_category * week * fiber_concentration_numeric, data = df)
    } else {
      model <- lme4::lmer(
        outcome ~ treatment_category * week * fiber_concentration_numeric + 
          (1 | tank) + (1 | experiment_batch), 
        data = df
      )
    }
    family_used <- "gaussian"
    dispersion_ratio <- NA
  }
  
  return(list(
    model = model, 
    family = family_used, 
    dispersion_ratio = dispersion_ratio
  ))
}

# ===== RECOVERY ANALYSIS =====
# Recovery period analysis - separate analysis for week 6 vs week 5
analyze_recovery <- function(df, endpoint_name) {
  recovery_data <- df %>%
    dplyr::filter(week %in% c(5, 6)) %>%
    dplyr::mutate(
      is_recovery = factor(ifelse(week == 6, "Recovery", "Exposure"),
                          levels = c("Exposure", "Recovery")),
      week_comparison = "Week5_vs_Week6"
    )
  
  if (nrow(recovery_data) == 0) {
    return(list(model = NULL, data = recovery_data, message = "No recovery data available"))
  }
  
  # Recovery model: compare week 6 to week 5
  tryCatch({
    if (endpoint_name == "mf_counts") {
      recovery_model <- lme4::glmer(
        outcome ~ treatment_category * is_recovery + 
          (1 | tank) + (1 | experiment_batch),
        family = poisson(),
        data = recovery_data
      )
    } else {
      recovery_model <- lme4::lmer(
        outcome ~ treatment_category * is_recovery + 
          (1 | tank) + (1 | experiment_batch), 
        data = recovery_data
      )
    }
    
    return(list(model = recovery_model, data = recovery_data))
  }, error = function(e) {
    return(list(model = NULL, data = recovery_data, error = e$message))
  })
}

# ===== TRANSLOCATION ANALYSIS =====
# Translocation analysis for mf_counts across tissues
analyze_translocation <- function(df) {
  df %>%
    dplyr::filter(endpoint == "mf_counts") %>%
    dplyr::group_by(fiber_type, treatment_category, week, tissue_type) %>%
    dplyr::summarise(
      mean_count = mean(outcome, na.rm = TRUE),
      median_count = median(outcome, na.rm = TRUE),
      se_count = sd(outcome, na.rm = TRUE) / sqrt(dplyr::n()),
      n_samples = dplyr::n(),
      n_nonzero = sum(outcome > 0, na.rm = TRUE),
      prop_nonzero = n_nonzero / n_samples,
      .groups = "drop"
    ) %>%
    dplyr::arrange(fiber_type, treatment_category, week, tissue_type)
}

# ===== MIXED EFFECTS MODEL BUILDING =====
build_lmer_formula <- function(include_three_way = TRUE, dose_by_fiber = FALSE, 
                              dose_by_treat = FALSE, include_recovery = FALSE,
                              random_terms = c("tank", "experiment_batch")) {
  
  # Base formula
  base_formula <- "outcome ~ fiber_type * chem_treatment + week + fiber_concentration + 
                    fiber_type:week + chem_treatment:week + is_control * fiber_type + dose_log10"
  
  # Add optional interactions
  if (include_three_way) {
    base_formula <- paste(base_formula, "+ fiber_type:chem_treatment:week")
  }
  
  if (dose_by_fiber) {
    base_formula <- paste(base_formula, "+ dose_log10:fiber_type")
  }
  
  if (dose_by_treat) {
    base_formula <- paste(base_formula, "+ dose_log10:chem_treatment")
  }
  
  if (include_recovery) {
    base_formula <- paste(base_formula, "+ is_recovery + is_recovery:chem_treatment + is_recovery:dose_log10")
  }
  
  # Add random effects
  random_effects <- paste(sprintf("(1|%s)", random_terms), collapse = " + ")
  final_formula <- paste(base_formula, "+", random_effects)
  
  return(as.formula(final_formula))
}

# ===== MODEL DIAGNOSTICS =====
check_model_assumptions <- function(model, endpoint_name) {
  diagnostics <- list()
  
  if (endpoint_name == "mf_counts") {
    # Count model diagnostics
    fitted_vals <- fitted(model)
    residuals_vals <- residuals(model, type = "pearson")
    
    diagnostics$dispersion <- var(residuals_vals) / mean(residuals_vals)
    diagnostics$zero_inflation <- sum(model$model$outcome == 0) / length(model$model$outcome)
    
  } else {
    # Linear model diagnostics
    diagnostics$normality_test <- shapiro.test(residuals(model)[1:min(5000, length(residuals(model)))])
    diagnostics$residual_plots <- list(
      fitted_vs_residuals = plot(fitted(model), residuals(model)),
      qq_plot = qqnorm(residuals(model))
    )
  }
  
  return(diagnostics)
}

# ===== CACHING HELPERS FOR PERFORMANCE =====
# Cache expensive model fitting operations
cached_model_fit <- function(df, endpoint_name, model_type, cache_key = NULL) {
  if (is.null(cache_key)) {
    # Generate cache key from data characteristics
    cache_key <- digest::digest(list(
      nrow(df), 
      names(df), 
      endpoint_name, 
      model_type,
      range(df$outcome, na.rm = TRUE)
    ))
  }
  
  # Check if model already cached (in production, use proper caching system)
  fit_endpoint_model(df, endpoint_name, model_type)
}