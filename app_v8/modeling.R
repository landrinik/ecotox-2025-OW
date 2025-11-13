# Statistical Modeling Functions
# Optimized for count data and mixed effects

# Always create a centered numeric time covariate from 'week' (or 'time_wk' if present)
prepare_time_wk <- function(df) {
  # prefer time_wk if present, else coerce week
  time_raw <- if ("time_wk" %in% names(df)) df$time_wk else df$week
  
  # coerce to numeric if needed
  if (!is.numeric(time_raw)) {
    time_raw <- suppressWarnings(as.numeric(as.character(time_raw)))
  }
  
  # if everything failed but df$week is numeric, use it
  if (all(is.na(time_raw)) && is.numeric(df$week)) time_raw <- df$week
  
  # center without scaling (keep slope in per-week units)
  df$time_wk_z <- as.numeric(scale(time_raw, center = TRUE, scale = FALSE))
  
  df
}


# ===== ENHANCED MODEL FITTING =====
# Main fitter gains a switch for the random slope
fit_endpoint_model <- function(df,
                               endpoint_name,
                               model_type = "lmer",
                               include_tank_time_slope = TRUE) {
  df <- prepare_time_wk(df)
  
  # Only tank random effects; optional time slope in tank
  rand_tank <- if (isTRUE(include_tank_time_slope)) {
    "(1 + time_wk_z | tank)"     # random intercept + week slope per tank
  } else {
    "(1 | tank)"                 # random intercept per tank only
  }
  
  if (endpoint_name == "mf_counts" && model_type %in% c("glmer", "glm")) {
    mean_val <- mean(df$outcome, na.rm = TRUE)
    var_val  <- var(df$outcome, na.rm = TRUE)
    dispersion_ratio <- ifelse(mean_val > 0, var_val / mean_val, 1)
    
    if (dispersion_ratio > 2) {
      # Negative binomial counts
      if (model_type == "glm") {
        model <- MASS::glm.nb(
          outcome ~ treatment_category * week * log(fiber_concentration_numeric + 1),
          data = df
        )
        family_used <- "negative_binomial"
      } else {
        model <- lme4::glmer.nb(
          stats::as.formula(
            paste(
              "outcome ~ treatment_category * week * log(fiber_concentration_numeric + 1) +",
              rand_tank
            )
          ),
          data = df,
          control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
        )
        family_used <- "negative_binomial"
      }
    } else {
      # Poisson counts
      if (model_type == "glm") {
        model <- stats::glm(
          outcome ~ treatment_category * week * log(fiber_concentration_numeric + 1),
          family = poisson(), data = df
        )
        family_used <- "poisson"
      } else {
        model <- lme4::glmer(
          stats::as.formula(
            paste(
              "outcome ~ treatment_category * week * log(fiber_concentration_numeric + 1) +",
              rand_tank
            )
          ),
          family = poisson(), data = df,
          control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
        )
        family_used <- "poisson"
      }
    }
  } else {
    # Continuous endpoints
    if (model_type == "lm") {
      model <- stats::lm(
        outcome ~ treatment_category * week * fiber_concentration_numeric,
        data = df
      )
      family_used <- "gaussian"
      dispersion_ratio <- NA_real_
    } else {
      model <- lme4::lmer(
        stats::as.formula(
          paste(
            "outcome ~ treatment_category * week * fiber_concentration_numeric +",
            rand_tank
          )
        ),
        data = df,
        control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
      )
      family_used <- "gaussian"
      dispersion_ratio <- NA_real_
    }
  }
  
  list(
    model = model,
    family = family_used,
    dispersion_ratio = if (exists("dispersion_ratio")) dispersion_ratio else NA_real_
  )
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
# Optional: formula builder with slope switch (used wherever you construct formulas)
build_lmer_formula <- function(include_three_way = TRUE,
                               dose_by_fiber = FALSE,
                               dose_by_treat = FALSE,
                               include_recovery = FALSE,
                               random_intercepts = character(0),
                               include_time_slope = TRUE,
                               slope_var = "time_wk_z",
                               slope_group = "tank") {
  base_formula <- "outcome ~ fiber_type * chem_treatment + week + fiber_concentration +
                   fiber_type:week + chem_treatment:week + is_control * fiber_type + dose_log10"
  if (isTRUE(include_three_way)) base_formula <- paste(base_formula, "+ fiber_type:chem_treatment:week")
  if (isTRUE(dose_by_fiber))    base_formula <- paste(base_formula, "+ dose_log10:fiber_type")
  if (isTRUE(dose_by_treat))    base_formula <- paste(base_formula, "+ dose_log10:chem_treatment")
  if (isTRUE(include_recovery)) base_formula <- paste(base_formula, "+ is_recovery + is_recovery:chem_treatment + is_recovery:dose_log10")
  
  rand_parts <- c(
    if (isTRUE(include_time_slope)) sprintf("(1 + %s | %s)", slope_var, slope_group) else "(1 | tank)",
    if (length(random_intercepts)) sprintf("(1|%s)", random_intercepts) else NULL
  )
  as.formula(paste(base_formula, "+", paste(rand_parts, collapse = " + ")))
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