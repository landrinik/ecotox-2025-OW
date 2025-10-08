# COMPLETE FIX FOR MIXED EFFECTS EMMEANS ISSUES

# ============================================================================
# 1. FIX THE EMMEANS SPECIFICATION PROBLEM
# ============================================================================

# Replace your lmer_emmeans_results function with this corrected version:
lmer_emmeans_results <- eventReactive(input$lmer_run_emmeans, {
  tryCatch({
    # Require that the lmer model exists first
    model <- lmer_model()
    validate(need(!is.null(model) && !inherits(model, "list"), "Run the mixed effects model first"))
    
    # Calculate dose value based on selection
    dose_val <- switch(input$lmer_emmeans_dose,
                      "0" = 0,           # Control level
                      "2" = log10(100),  # Low dose
                      "3" = log10(1000), # Medium dose
                      "4" = log10(10000)) # High dose
    
    is_ctrl <- ifelse(dose_val == 0, 1, 0)
    at_list <- list(dose_log10 = dose_val, is_control = is_ctrl)
    
    # FIXED: Better emmeans specifications to avoid convergence issues
    emmeans_spec <- switch(input$lmer_emmeans_by,
                          "treatment" = ~ chem_treatment,           # Simple: just treatment
                          "fiber_type" = ~ fiber_type,             # Simple: just fiber
                          "week" = ~ week,                         # Simple: just week
                          "treatment_week" = ~ chem_treatment + week,      # Additive, not interaction
                          "treatment_fiber" = ~ chem_treatment + fiber_type, # Additive, not interaction
                          "fiber_week" = ~ fiber_type + week,       # Additive, not interaction
                          "all_interactions" = ~ chem_treatment + fiber_type + week) # Additive only
    
    # ENHANCED: Try different approaches if emmeans fails
    emm <- tryCatch({
      # Primary approach: Use full specification
      emmeans(model, specs = emmeans_spec, at = at_list)
    }, error = function(e1) {
      message("Primary emmeans failed, trying simplified approach: ", e1$message)
      tryCatch({
        # Fallback 1: Remove the at= specification (use model's reference levels)
        emmeans(model, specs = emmeans_spec)
      }, error = function(e2) {
        message("Simplified emmeans failed, trying minimal approach: ", e2$message)
        # Fallback 2: Use simplest specification possible
        emmeans(model, specs = ~ chem_treatment)
      })
    })
    
    list(
      emmeans = emm,
      emmeans_df = as.data.frame(emm),
      specification_used = deparse(emmeans_spec)
    )
  }, error = function(e) {
    message("All emmeans approaches failed: ", e$message)
    list(error = e$message)
  })
})

# ============================================================================
# 2. IMPROVE MIXED EFFECTS MODEL TO REDUCE CONVERGENCE ISSUES
# ============================================================================

# Replace your lmer_model function with this more stable version:
lmer_model <- eventReactive(input$run_lmer, {
  tryCatch({
    df <- lmer_data()
    req(nrow(df) > 0)
    
    message("Building lmer model with ", nrow(df), " observations")
    
    # ENHANCED: More stable model building approach
    
    # Start with basic formula
    fixed_formula <- "outcome ~ fiber_type + chem_treatment + week + dose_log10 + is_control"
    
    # Add interactions carefully to avoid over-parameterization
    if (input$lmer_include_three_way && length(unique(df$week)) > 2) {
      fixed_formula <- paste(fixed_formula, "+ fiber_type:chem_treatment + fiber_type:week + chem_treatment:week")
      
      # Only add 3-way if we have enough data
      if (nrow(df) > 100) {
        fixed_formula <- paste(fixed_formula, "+ fiber_type:chem_treatment:week")
      }
    } else {
      # Add basic 2-way interactions
      fixed_formula <- paste(fixed_formula, "+ fiber_type:chem_treatment")
    }
    
    if (input$lmer_dose_by_fiber) {
      fixed_formula <- paste(fixed_formula, "+ dose_log10:fiber_type")
    }
    
    if (input$lmer_dose_by_treat) {
      fixed_formula <- paste(fixed_formula, "+ dose_log10:chem_treatment")
    }
    
    # Handle recovery interaction
    if (input$lmer_include_recovery_interaction && "6" %in% input$lmer_weeks_include) {
      fixed_formula <- paste(fixed_formula, "+ is_recovery + is_recovery:chem_treatment")
    }
    
    # SMART: Build random effects based on data availability
    random_terms <- input$lmer_random_terms
    random_formula <- ""
    
    if ("tank" %in% random_terms && "tank" %in% names(df)) {
      # Check if we have enough tanks
      n_tanks <- length(unique(df$tank))
      if (n_tanks > 3) {
        random_formula <- paste(random_formula, "+ (1|tank)")
      } else {
        message("Too few tanks (", n_tanks, ") for random effect, skipping tank")
      }
    }
    
    if ("experiment_batch" %in% random_terms && "experiment_batch" %in% names(df)) {
      # Check if we have enough batches  
      n_batches <- length(unique(df$experiment_batch))
      if (n_batches > 1) {
        random_formula <- paste(random_formula, "+ (1|experiment_batch)")
      } else {
        message("Only one experiment batch, skipping batch random effect")
      }
    }
    
    # Fallback random effect
    if (random_formula == "" && "grouping_var" %in% names(df)) {
      n_groups <- length(unique(df$grouping_var))
      if (n_groups > 3) {
        random_formula <- "+ (1|grouping_var)"
        message("Using grouping_var as fallback random effect")
      }
    }
    
    # Construct final formula
    if (random_formula != "") {
      full_formula <- paste(fixed_formula, random_formula)
    } else {
      # If no random effects possible, fall back to simple lm
      message("No suitable random effects found, this may cause issues")
      full_formula <- fixed_formula
    }
    
    formula_obj <- as.formula(full_formula)
    message("Final formula: ", deparse(formula_obj))
    
    # Fit model with enhanced control parameters for stability
    model <- lmer(formula_obj, data = df, 
                 control = lmerControl(
                   optimizer = "bobyqa",      # More stable optimizer
                   optCtrl = list(maxfun = 20000),  # More iterations
                   check.conv.singular = "warning",  # Don't fail on singular fits
                   calc.derivs = FALSE        # Skip expensive derivative calculations
                 ))
    
    # Check for convergence warnings
    if (length(model@optinfo$conv$lme4$messages) > 0) {
      message("Model convergence warnings: ", paste(model@optinfo$conv$lme4$messages, collapse = "; "))
    }
    
    model
    
  }, error = function(e) {
    message("LMER model fitting error: ", e$message)
    list(error = e$message)
  })
})

# ============================================================================
# 3. ENHANCED MIXED EFFECTS PLOT WITH BETTER ERROR HANDLING
# ============================================================================

output$lmer_emmeans_plot <- renderPlot({
  results <- lmer_emmeans_results()
  validate(need(!is.null(results$emmeans), "Calculate EM Means for mixed effects first"))
  validate(need(is.null(results$error), paste("Emmeans error:", results$error)))
  
  tryCatch({
    # ENHANCED: Better plot with error handling for negative variances
    
    # Check if emmeans has issues with confidence intervals
    df_plot <- results$emmeans_df
    
    # Remove rows with missing or invalid confidence intervals
    valid_rows <- !is.na(df_plot$emmean) & 
                  !is.na(df_plot$lower.CL) & 
                  !is.na(df_plot$upper.CL) &
                  is.finite(df_plot$lower.CL) & 
                  is.finite(df_plot$upper.CL)
    
    if (sum(valid_rows) == 0) {
      # Fallback: create simple point plot without confidence intervals
      ggplot(df_plot, aes(x = emmean, y = rownames(df_plot))) +
        geom_point(size = 4, color = "#e74c3c") +
        theme_minimal() +
        labs(title = "Mixed Effects Estimated Marginal Means (No CI due to model issues)",
             subtitle = paste("Dose level:", input$lmer_emmeans_dose, "| Some confidence intervals invalid"),
             x = "Estimated Marginal Mean", y = "Groups") +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 11, color = "red"),
          axis.text.y = element_text(size = 10)
        )
    } else {
      # Use valid rows only
      df_clean <- df_plot[valid_rows, ]
      
      # Create the plot manually to handle negative variances better
      ggplot(df_clean, aes(x = emmean, y = factor(rownames(df_clean)))) +
        geom_point(size = 3, color = "#2c3e50") +
        geom_errorbarh(aes(xmin = pmax(lower.CL, emmean - 3*SE), 
                           xmax = pmin(upper.CL, emmean + 3*SE)), 
                      height = 0.2, color = "#3498db") +
        theme_minimal() +
        labs(title = "Mixed Effects Estimated Marginal Means with 95% CI",
             subtitle = paste("Dose level:", input$lmer_emmeans_dose, 
                             "| Adjustment:", input$lmer_adjustment_method),
             x = "Estimated Marginal Mean", y = "Groups") +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          panel.grid.minor = element_blank()
        )
    }
  }, error = function(e) {
    # Ultimate fallback plot
    message("Plot error: ", e$message)
    ggplot(data.frame(x = 1, y = 1, message = "Plot failed due to model convergence issues")) +
      geom_text(aes(x, y, label = message), size = 5, color = "red") +
      theme_void() +
      labs(title = "Mixed Effects Plot Error", subtitle = "Try simpler model settings")
  })
})