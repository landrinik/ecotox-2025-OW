library(tidyverse)
library(readxl)
library(writexl)
library(ggplot2)
library(yaml)
library(drc)
library(diffdf)
library(janitor)
library(dplyr)
library(stringr)

# Source your QC function first
source("flag_plate_qc.R")

# Source analyzed assay functions with QC enhancements
source("my_functions_v3_avg.R")

# Load full cleaned data table
mydata <- read_excel("final_data_clean.xlsx") %>%
  clean_names()

# Load tissue weights for activity normalization as needed
tissue_weights_clean <- read_excel("tissue_weights_clean.xlsx")  %>%
  clean_names()

# Define standard concentrations for each assay type
sod_concentrations <- c(92, 46, 23, 11.5, 5.8, 2.9, 1.4, 0)
cat_concentrations <- c(5, 2.5, 1.25, 0.625, 0.313, 0.156, 0)
acp_concentrations <- c(1000, 600, 300, 0)
bradford_concentrations <- c(2, 1.5, 1, 0.75, 0.5, 0.25, 0.125, 0)

# Load config YAML specifying filters for assays, fibers, weeks, etc.
config <- read_yaml("config.yml")

# Compute global summary averages for blanks and standards from raw data column 12
global_standards <- mydata %>%
  dplyr::select(assay_type, fiber_type, sample_week, sample_type, plate_replicate, reading, wells, x12) %>%
  group_by(assay_type, fiber_type, sample_week, sample_type, reading, wells) %>%
  summarise(avg_value = mean(x12, na.rm = TRUE), .groups = "drop")

# Build analysis job list of unique plates excluding blanks
full_job_list <- mydata %>%
  dplyr::select(
    assay_to_filter = assay_type,
    fiber_to_filter = fiber_type,
    week_to_filter = sample_week,
    sample_type_to_filter = sample_type,
    replicate_to_filter = plate_replicate
  ) %>%
  distinct() %>%
  filter(replicate_to_filter != "BlankRep") %>%
  mutate(
    concentrations_vector = case_when(
      assay_to_filter == "SOD" ~ list(sod_concentrations),
      assay_to_filter == "CAT" ~ list(cat_concentrations),
      assay_to_filter == "ACP" ~ list(acp_concentrations),
      assay_to_filter == "Bradford" ~ list(bradford_concentrations),
      assay_to_filter %in% c("ALP", "AChE") ~ list(NULL)
    )
  )

# Filter by config selections if they exist
job_list <- full_job_list %>%
  filter(
    is.null(config$assays) | assay_to_filter %in% config$assays,
    is.null(config$fibers) | fiber_to_filter %in% config$fibers,
    is.null(config$weeks) | week_to_filter %in% config$weeks,
    is.null(config$samples) | sample_type_to_filter %in% config$samples,
    is.null(config$replicates) | replicate_to_filter %in% config$replicates
  )

# Run analysis on all filtered jobs, returning full analysis_output list
all_results <- job_list %>%
  mutate(
    analysis_output = pmap(
      .l = list(
        assay_to_filter = assay_to_filter,
        concentrations_vector = concentrations_vector,
        fiber_to_filter = fiber_to_filter,
        week_to_filter = week_to_filter,
        sample_type_to_filter = sample_type_to_filter,
        replicate_to_filter = replicate_to_filter
      ),
      .f = analyze_assay_data,
      full_dataset = mydata,
      global_standards = global_standards,
      tissue_weights = tissue_weights_clean
    )
  )

# Filter to keep only successful (non-null) results
successful_results <- all_results %>%
  filter(!sapply(analysis_output, is.null))

# Pull out long format final data tables with QC flags
final_data_list <- successful_results %>%
  pull(analysis_output) %>%
  map(~ .x$final_data)

# Combine all results into a single long table with replicate CV and QC flags
final_long_table_avgstd <- bind_rows(final_data_list) %>%
  clean_names() 

# Create Bradford lookup (protein mg/ml) per sample well
bradford_lookup <- final_long_table_avgstd %>%
  filter(assay_type == "Bradford") %>%
  dplyr::select(fiber_type, sample_week, sample_type, plate_replicate, well_name,
         protein_conc_mg_ml = calculated_concentration)

# Normalize CAT, SOD, ACP, ALP assays
assays_to_normalize <- c("CAT", "SOD", "ACP", "ALP")

final_long_table_normalized <- final_long_table_avgstd %>%
  left_join(
    bradford_lookup,
    by = c("fiber_type", "sample_week", "sample_type", "plate_replicate", "well_name")
  ) %>%
  mutate(
    normalized_concentration = if_else(
      assay_type %in% assays_to_normalize & !is.na(protein_conc_mg_ml) & protein_conc_mg_ml != 0,
      calculated_concentration / protein_conc_mg_ml,
      calculated_concentration
    )
  ) %>%
  mutate(calculated_concentration = normalized_concentration) %>%
  dplyr::select(fiber_type,
                assay_type,
                sample_week,
                sample_type,
                plate_replicate,
                well_name,
                calculated_concentration,
                qc_flag)

# Save to Excel file for downstream use
write_xlsx(final_long_table_normalized, "final_long_table_normalized.xlsx")
write_xlsx(final_long_table_avgstd, "final_long_table_avgstd.xlsx")
# View(bradford_lookup)
# View(final_long_table_avgstd)
# View the final long table to confirm it's correct
View(final_long_table_normalized)