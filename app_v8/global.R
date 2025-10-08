
# OPTIMIZED global.R - Uses RDS files for fast loading

# ===== LIBRARIES =====
library(shiny)
library(tidyverse)
library(DT)
library(glue)
library(yaml)
library(stringr)
library(shinyjs)
library(emmeans)
library(broom)
library(broom.mixed)
library(car)
library(lme4)
library(MASS)
library(performance)

# ===== UTILITY FUNCTIONS =====
`%||%` <- function(a, b) if (!is.null(a)) a else b
filter_choices <- function(x) unlist(x[!startsWith(as.character(x), "#")])

# ===== LOAD OPTIMIZED DATA =====
message("Loading optimized data...")

# Load pre-processed data (much faster than Excel)
final_data <- readRDS("data/assay_data.rds")
physical_master <- readRDS("data/physical_data.rds")
tissue_weights <- readRDS("data/tissue_weights.rds")

# Extract UI choices
week_choices_assay <- sort(unique(na.omit(final_data$week)))
endpoint_choices_assay <- sort(unique(final_data$assay_type))
endpoint_choices_physical <- sort(unique(physical_master$endpoint))

# Create lookup tables
tank_fiber_lookup <- final_data %>%
  dplyr::select(tank, fiber_type, fiber_concentration, treatment) %>%
  distinct() %>%
  arrange(tank)

# Load config
config <- tryCatch({
  yaml::read_yaml("config.yml")
}, error = function(e) {
  list(
    assays = endpoint_choices_assay,
    fibers = c("cotton", "pet"),
    samples = c("tissue", "gland", "gills")
  )
})

message("✅ Global data loading complete!")
message("📊 Data summary:")
message("  Assay data: ", nrow(final_data), " rows")
message("  Physical data: ", nrow(physical_master), " rows")
message("  Memory usage: ", round(as.numeric(object.size(final_data) + object.size(physical_master))/1024^2, 1), " MB")

