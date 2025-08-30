# Flagging function to identify problematic wells
flag_plate_qc <- function(results_df, flag_thresholds = list()) 
  
{
  if (!"calculated_concentration" %in% colnames(results_df)) {
    stop("QC error: 'calculated_concentration' column missing in input data.")
  }
  if (!"well_name" %in% colnames(results_df)) {
    stop("QC error: 'well_name' column missing in input data.")
  }
  
  # Set default thresholds if missing
  if (is.null(flag_thresholds$min_concentration)) flag_thresholds$min_concentration <- 0
  if (is.null(flag_thresholds$max_concentration)) flag_thresholds$max_concentration <- Inf
  
  qc_flags <- list()
  
  # Detect out-of-bound values
  out_of_bounds <- results_df %>%
    dplyr::filter(calculated_concentration < flag_thresholds$min_concentration |
                    calculated_concentration > flag_thresholds$max_concentration)
  if (nrow(out_of_bounds) > 0) {
    qc_flags$out_of_bounds <- dplyr::select(out_of_bounds, well_name, calculated_concentration)
  }
  
  # Detect negative concentrations
  negative_vals <- results_df %>%
    dplyr::filter(calculated_concentration < 0)
  if (nrow(negative_vals) > 0) {
    qc_flags$negative <- dplyr::select(negative_vals, well_name, calculated_concentration)
  }
  
  # Detect NA concentrations
  na_vals <- results_df %>%
    dplyr::filter(is.na(calculated_concentration))
  if (nrow(na_vals) > 0) {
    qc_flags$na <- dplyr::select(na_vals, well_name)
  }
  
  return(qc_flags)
}

# Helper to provide assay-specific QC thresholds
qc_thresholds_for_assay <- function(assay_type) 
{
  switch(
    assay_type,
    "SOD" = list(min_concentration = 0, max_concentration = Inf),
    "CAT" = list(min_concentration = 0, max_concentration = Inf),
    "ACP" = list(min_concentration = 0, max_concentration = 1000),
    "ALP" = list(min_concentration = 0, max_concentration = 1000),
    "AChE" = list(min_concentration = 0, max_concentration = 1000),
    "Bradford" = list(min_concentration = 0, max_concentration = 10),
    # Default
    list(min_concentration = 0, max_concentration = Inf)
  )
}

# Wrapper to add QC flags to results based on assay type
apply_qc_flags <- function(results_df, assay_type = NULL) 
{
  if (!"calculated_concentration" %in% colnames(results_df)) {
    stop("QC error: 'calculated_concentration' column missing in input data.")
  }
  if (nrow(results_df) == 0) return(results_df)
  
  flag_thresholds <- qc_thresholds_for_assay(assay_type)
  
  qc_list <- flag_plate_qc(results_df, flag_thresholds)
  
  # Initialize the qc_flag column if missing
  if (!"qc_flag" %in% colnames(results_df)) {
    results_df$qc_flag <- FALSE
  }
  
  # Mark wells with any QC issues
  flagged_wells <- unique(c(
    if (!is.null(qc_list$out_of_bounds)) qc_list$out_of_bounds$well_name else character(0),
    if (!is.null(qc_list$negative)) qc_list$negative$well_name else character(0),
    if (!is.null(qc_list$na)) qc_list$na$well_name else character(0)
  ))
  
  results_df$qc_flag[results_df$well_name %in% flagged_wells] <- TRUE
  
  return(results_df)
}
