library(tidyverse)
library(readxl)
library(writexl)

# Create a list of all files to read (exclude the output of the script in case 
# changes need to be made)
all_files <- list.files(
  path = "Plates data/",
  pattern = ".xlsx", 
  recursive = TRUE,
  full.names = TRUE)

files_to_process <- all_files[!all_files %in% "final_data_clean.xlsx"]

# Separate the first file as a standard for all other files(used for names)
first_file <- files_to_process[1]
other_files <- files_to_process[-1]

# Reading the first plate data to make a tibble
first_plate_data <- read_excel(
  first_file,
  col_types = c("text", rep("numeric", 12))
  )

# Gather names from first tibble
correct_names <- colnames(first_plate_data)

# Use map function to read all files in other_files list, skipping the first 
# line and avoid using new names
other_plates_data <- map(
  other_files,
  ~ read_excel(
    .x,
    skip = 1,
    col_names = FALSE,
    col_types = c("text", rep("numeric", 12))
  )
)


other_plates_data_named <- map(other_plates_data,  set_names, correct_names)

all_plates_list <- c(list(first_plate_data), other_plates_data_named)

names(all_plates_list) <- files_to_process

# Connect all files with bind_rows and add a source file column
final_data <- bind_rows(all_plates_list, .id = "source_file")

# Use mutate function with case_when to add columns based on the file name
# using str_detect to look for specific keyword in the file name
final_data_cleaned <- final_data %>%
  mutate(
    FiberType = case_when(
      str_detect(source_file, "Cotton") ~ "Cotton",
      str_detect(source_file, "PET")    ~ "PET",
      TRUE ~ "Unknown"
    ),
    
    AssayType = case_when(
      str_detect(source_file, "AChE") ~ "AChE",
      str_detect(source_file, "ACP")  ~ "ACP",
      str_detect(source_file, "ALP")  ~ "ALP",
      str_detect(source_file, "Bradford") ~ "Bradford",
      str_detect(source_file, "CAT")  ~ "CAT",
      str_detect(source_file, "SOD")  ~ "SOD",
      TRUE ~ "Unknown"
    ),
    
    SampleWeek = case_when(
      str_detect(source_file, "1-3") ~ "1-3",
      str_detect(source_file, "5-6") ~ "5-6",
      TRUE ~ "Unknown"
    ),
    
    SampleType = case_when(
      str_detect(source_file, "Gills") ~ "Gills",
      str_detect(source_file, "Gland") ~ "Gland",
      TRUE ~ "Hemolymph"
    ),
    
    PlateReplicate = case_when(
      str_detect(source_file, "P1|Plate1") ~ "Rep1",
      str_detect(source_file, "P2|Plate2") ~ "Rep2",
      str_detect(source_file, "Blank") ~ "BlankRep",
      TRUE ~ "Blankplate"
    ),
    
    Reading = case_when(
      str_detect(source_file, "FirstRead") ~ "FirstRead",
      str_detect(source_file, "SecondRead") ~ "SecondRead",
      TRUE ~ "SingleRead"
    )
  ) %>%

# Remove the source file column which was used only to help
dplyr::select(-source_file) %>% 

relocate(FiberType, AssayType, SampleWeek, SampleType, PlateReplicate, Reading)

write_xlsx(final_data_cleaned, "final_data_clean.xlsx")