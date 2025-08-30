library(tidyverse)
library(writexl)
library(readxl)

# 1. Get the names of all the sheets (plates) in your tissue weight Excel file
tissue_sheets <- readxl::excel_sheets("deepwells_weights.xlsx")

# 2. Read and tidy every sheet, then bind them together
tissue_weights_long <- map_df(tissue_sheets, ~{
  
  # Read one sheet
  raw_weights <- read_excel("deepwells_weights.xlsx", sheet = .x)
  
  # Tidy the data
  tidy_weights <- raw_weights %>%
    pivot_longer(
      cols = -Wells, # Pivot all columns except the 'Well' identifier
      names_to = "Well_Column",
      values_to = "Tissue_Weight_mg"
    ) %>%
    unite(col = "well_name", Wells, Well_Column, sep = "") %>%
    # Add a column for the plate name from the sheet name
    mutate(plate_id = .x) 
  
  return(tidy_weights)
})

# 3. Separate the plate_id into your descriptive columns
# This now creates a simpler table without a replicate column
tissue_weights_final <- tissue_weights_long %>%
  separate(
    plate_id, 
    into = c("AssayType", "FiberType", "SampleWeek", "SampleType"), 
    sep = "_"
  ) %>%
  
  # Add this filter to remove the empty wells
  filter(!well_name %in% c("E11", "F11", "G11", "H11")) %>%
  
  # Keep only the columns needed for the join
  dplyr::select(assay_type = AssayType,
                fiber_type = FiberType,
                sample_week = SampleWeek,
                sample_type = SampleType,
                well_name,
                tissue_weight_mg = Tissue_Weight_mg)

# 4. Save the final, clean data to a new file
write_xlsx(tissue_weights_final, "tissue_weights_clean.xlsx")