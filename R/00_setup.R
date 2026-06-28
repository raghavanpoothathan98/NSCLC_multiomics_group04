############################################################
# 00_setup.R
# Project setup for NSCLC multi-omics project
############################################################

# -----------------------------
# 1. Install and load packages
# -----------------------------

packages <- c(
  "here",
  "data.table",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "readr",
  "matrixStats",
  "stringr"
)

installed <- packages %in% rownames(installed.packages())

if (any(!installed)) {
  install.packages(packages[!installed])
}

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(matrixStats)
  library(stringr)
})

# -----------------------------
# 2. Define project folders
# -----------------------------

project_dir <- here::here()

data_dir <- here::here("data")
raw_dir <- here::here("data", "raw")
processed_dir <- here::here("data", "processed")

luad_raw_dir <- file.path(raw_dir, "LUAD", "luad_tcga_pan_can_atlas_2018")
lusc_raw_dir <- file.path(raw_dir, "LUSC", "lusc_tcga_pan_can_atlas_2018")

results_dir <- here::here("results")
figures_dir <- here::here("results", "figures")
tables_dir <- here::here("results", "tables")

notebooks_dir <- here::here("notebooks")
report_dir <- here::here("report")
presentation_dir <- here::here("presentation")

# -----------------------------
# 3. Create output folders
# -----------------------------

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 4. Define raw input folders
# -----------------------------

raw_paths <- list(
  luad_dir = luad_raw_dir,
  lusc_dir = lusc_raw_dir
)

# -----------------------------
# 5. Define generated output paths
# -----------------------------
# These files are created by 01_load_data.R.
# They are NOT required before running the workflow.

paths <- list(
  loaded_matched_rds = file.path(processed_dir, "01_loaded_matched_data.rds"),
  labels = file.path(processed_dir, "sample_labels.tsv"),
  mrna = file.path(processed_dir, "matched_mrna.tsv"),
  methylation = file.path(processed_dir, "matched_methylation.tsv"),
  rppa = file.path(processed_dir, "matched_rppa.tsv"),
  qc_model_input = file.path(processed_dir, "02_qc_model_input_unimputed.rds"),
  eda_results = file.path(processed_dir, "03_eda_results.rds")
)

# -----------------------------
# 6. Print project information
# -----------------------------

message("Project root:")
message(project_dir)

message("\nRaw data folder:")
message(raw_dir)

message("\nLUAD raw folder:")
message(luad_raw_dir)

message("\nLUSC raw folder:")
message(lusc_raw_dir)

message("\nProcessed data folder:")
message(processed_dir)

message("\nResults folder:")
message(results_dir)

# -----------------------------
# 7. Check raw folders only
# -----------------------------

raw_folder_check <- c(
  LUAD_raw_folder = dir.exists(luad_raw_dir),
  LUSC_raw_folder = dir.exists(lusc_raw_dir)
)

message("\nRaw folder check:")
print(raw_folder_check)

if (!all(raw_folder_check)) {
  stop(
    "One or more raw data folders are missing. Expected:\n",
    luad_raw_dir, "\n",
    lusc_raw_dir, "\n",
    "Check data/raw/ folder names."
  )
}

message("\nSetup completed successfully.")