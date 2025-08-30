# LOAD LIBRARIES & DATA ---------------------------------------------
library(tidyverse)
library(readxl)
library(writexl)

# Load the long-format results from your main analysis script
final_long_table <- read_excel("final_long_table_avgstd.xlsx")
tissue_weights <- read_excel("tissue_weights_clean.xlsx") %>% clean_names()

# CV threshold for QC flagging
cv_threshold <- 15  # percent

# ANNOTATE DATA -----------------------------------------------------
# This first part creates all the necessary descriptive columns
annotated_data <- final_long_table %>%
  mutate(
    well_number = as.numeric(str_extract(well_name, "\\d+")),
    well_letter = str_extract(well_name, "[A-H]"),
    week_start = as.numeric(str_extract(sample_week, "^[0-9]+")),
    week_end = as.numeric(str_extract(sample_week, "[0-9]+$")),
    week = if_else(well_number <= 6, week_start, week_end),
    tank = case_when(
      well_number <= 5 & well_letter %in% c("A","B") ~ (well_number - 1)*4 + 1,
      well_number <= 5 & well_letter %in% c("C","D") ~ (well_number - 1)*4 + 2,
      well_number <= 5 & well_letter %in% c("E","F") ~ (well_number - 1)*4 + 3,
      well_number <= 5 & well_letter %in% c("G","H") ~ (well_number - 1)*4 + 4,
      well_number == 6 & well_letter %in% c("A","B") ~ 21,
      well_number == 6 & well_letter %in% c("C","D") ~ 1,
      well_number == 6 & well_letter %in% c("E","F") ~ 2,
      well_number == 6 & well_letter %in% c("G","H") ~ 3,
      well_number > 6 & well_letter %in% c("A","B") ~ (well_number - 7)*4 + 4,
      well_number > 6 & well_letter %in% c("C","D") ~ (well_number - 7)*4 + 5,
      well_number > 6 & well_letter %in% c("E","F") ~ (well_number - 7)*4 + 6,
      well_number > 6 & well_letter %in% c("G","H") ~ (well_number - 7)*4 + 7
    ),
    fiber_concentration = case_when(
      tank %in% c(1,2,3) ~ "Control",
      tank %in% c(4,5,6,13,14,15) ~ "100",
      tank %in% c(7,8,9,16,17,18) ~ "1000",
      tank %in% c(10,11,12,19,20,21) ~ "10000"
    ),
    treatment = if_else(tank <= 12, "Treated", "Untreated"),
    fiber_group = if_else(
      fiber_concentration == "Control",
      "Control",
      paste(fiber_type, treatment)
    )
  )

# Summarize activity and calculate CV and QC flag
summary_with_qc <- annotated_data %>%
  group_by(fiber_group, assay_type, week, sample_type, fiber_concentration, plate_replicate) %>%
  summarise(
    mean_activity = mean(calculated_concentration, na.rm = TRUE),
    sd_activity = sd(calculated_concentration, na.rm = TRUE),
    n_samples = n(),
    cv = if_else(mean_activity == 0, NA_real_, 100 * sd_activity / mean_activity),
    cv_flag = cv > cv_threshold,
    .groups = "drop"
  ) %>%
  mutate(cv = round(cv, 2))

# Output flagged samples
flagged <- filter(summary_with_qc, cv_flag == TRUE)
if (nrow(flagged) > 0) {
  message(glue::glue("Warning: {nrow(flagged)} groups exceed CV threshold of {cv_threshold}%"))
  write_xlsx(flagged, "flagged_cv_samples.xlsx")
  View(flagged)
}

# Save full summary
write_xlsx(summary_with_qc, "summary_with_qc_flags.xlsx")
View(summary_with_qc)

# # SUMMARIZE BY FIBER CONCENTRATION ----------------------------------
# summary_by_concentration <- annotated_data %>%
#   group_by(
#     fiber_group, 
#     assay_type, 
#     week, 
#     sample_type, 
#     fiber_concentration, 
#     plate_replicate
#   ) %>%
#   summarise(
#     mean_activity = mean(calculated_concentration, na.rm = TRUE),
#     sd_activity = sd(calculated_concentration, na.rm = TRUE),
#     n_samples = n(),
#     .groups = "drop"
#   )
# 
# # View and save the new summary report
# View(summary_by_concentration)
# write_xlsx(summary_by_concentration, "summary_by_concentration.xlsx")