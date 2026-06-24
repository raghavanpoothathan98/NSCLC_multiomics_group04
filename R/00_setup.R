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
  "matrixStats"
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
})

# -----------------------------
# 2. Define project folders
# -----------------------------

project_dir <- here::here()

data_dir <- here::here("data")
raw_dir <- here::here("data", "raw")
processed_dir <- here::here("data", "processed")

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
# 4. Define matched data paths
# -----------------------------
# Matched files are stored in data/processed/

paths <- list(
  labels = file.path(processed_dir, "sample_labels.tsv"),
  mrna = file.path(processed_dir, "matched_mrna.tsv"),
  methylation = file.path(processed_dir, "matched_methylation.tsv"),
  rppa = file.path(processed_dir, "matched_rppa.tsv"),
  integrated = file.path(processed_dir, "integrated_feature_matrix.tsv"),
  missing_summary = file.path(processed_dir, "missing_value_summary.xlsx"),
  cv_results = file.path(processed_dir, "cv_results.tsv")
)

# -----------------------------
# 5. Print project information
# -----------------------------

message("Project root:")
message(project_dir)

message("\nProcessed data folder:")
message(processed_dir)

message("\nFiles found in processed folder:")
print(list.files(processed_dir, all.files = TRUE))

# -----------------------------
# 6. Check required files
# -----------------------------

required_paths <- paths[c(
  "labels",
  "mrna",
  "methylation",
  "rppa"
)]

file_check <- sapply(required_paths, file.exists)

message("\nRequired file check:")
print(file_check)

if (!all(file_check)) {
  missing_files <- names(file_check)[!file_check]
  
  stop(
    "Missing required file(s): ",
    paste(missing_files, collapse = ", "),
    "\nCheck filenames and location inside data/processed/."
  )
}

# -----------------------------
# 7. Optional file check
# -----------------------------

optional_paths <- paths[c(
  "integrated",
  "missing_summary",
  "cv_results"
)]

optional_check <- sapply(optional_paths, file.exists)

message("\nOptional file check:")
print(optional_check)

message("\nSetup completed successfully.")