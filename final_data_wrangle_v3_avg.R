library(tidyverse)
library(readxl)
library(writexl)
library(dplyr)
library(stringr)

# Load the new long-format results from the updated pipeline
final_long_table <- read_excel("final_long_table_avgstd.xlsx")

# Check if tissue_weights_clean exists, else create empty placeholder
tissue_weights_clean <- read_excel("tissue_weights_clean.xlsx")

# CV threshold for flagging
cv_threshold <- 15

# Perform final summarization and CV flagging
final_summary_cv <- final_long_table %>%
  
  mutate(well_name = stringr::str_remove(well_name, "x")) %>%
  
  # Join tissue weight data
  left_join(tissue_weights_clean, by = c("assay_type", "fiber_type", "sample_type", "sample_week", "well_name")) %>%
  
  # Add annotations
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
  ) %>%
  
  # Group and summarize mean activity, SD and calculate CV
  group_by(fiber_group, assay_type, week, sample_type, tank, well_letter, fiber_concentration, well_name) %>%
  summarise(
    mean_activity = mean(calculated_concentration, na.rm=TRUE),
    sd_activity_raw = sd(calculated_concentration, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sd_activity = round(sd_activity_raw, 4),
    cv_activity = if_else(mean_activity != 0, (sd_activity / mean_activity) * 100, NA_real_),
    cv_percent = round(cv_activity, 2),
    cv_flag = cv_percent > cv_threshold
  ) %>%
  dplyr::select(fiber_group, assay_type, week, sample_type, tank, well_letter, fiber_concentration, well_name,
         mean_activity, sd_activity, cv_percent, cv_flag)

# Show flagged results
flagged <- final_summary_cv %>% filter(cv_flag)

if(nrow(flagged) > 0){
  message(glue::glue("Warning: {nrow(flagged)} groups exceed CV threshold of {cv_threshold}%"))
  View(flagged)
  write_xlsx(flagged, "high_cv_flags.xlsx")
}

# Save full report
write_xlsx(final_summary_cv, "final_summary_report_cv_avg.xlsx")
View(final_summary_cv)

# # Optional ordering by column plates instead of rows
# ordered_summary_cv <- final_summary_cv %>%
#   mutate(
#     well_letter = factor(str_extract(well_name, "[A-H]"), levels = LETTERS[1:8]),
#     well_number = as.numeric(str_extract(well_name, "\\d+"))
#   ) %>%
#   arrange(well_number, well_letter)
# 
# write_xlsx(ordered_summary_cv, "final_summary_report_cv_avg.xlsx")
# View(ordered_summary_cv)
