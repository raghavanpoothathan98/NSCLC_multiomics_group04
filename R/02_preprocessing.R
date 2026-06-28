############################################################
# 02_preprocessing.R
# Leakage-aware preprocessing and basic QC
# NSCLC multi-omics project: LUAD vs LUSC
#
# Input:
#   data/processed/01_loaded_matched_data.rds
#
# Output:
#   data/processed/02_qc_model_input_unimputed.rds
#
# This script does:
#   - loads matched labels, metadata, mRNA, methylation, RPPA
#   - checks sample alignment
#   - summarizes missingness and variance
#   - removes high-missingness features
#   - removes zero/invalid variance features
#   - identifies candidate batch/source variables from metadata
#
# This script does NOT do:
#   - raw file loading
#   - sample matching
#   - global imputation
#   - global scaling
#   - global top-variable feature selection
#   - supervised feature selection
#   - batch correction
#
# Final modelling scripts must perform:
#   - imputation
#   - scaling
#   - top-variable filtering
#   - supervised feature selection
# inside each cross-validation training fold.
############################################################

# -----------------------------
# 1. Setup and load matched data
# -----------------------------

source("R/00_setup.R")

loaded_rds_path <- file.path(processed_dir, "01_loaded_matched_data.rds")

if (!file.exists(loaded_rds_path)) {
  stop(
    "01_loaded_matched_data.rds not found.\n",
    "Run source('R/01_load_data.R') first."
  )
}

loaded_data <- readRDS(loaded_rds_path)

labels <- loaded_data$labels
metadata <- loaded_data$metadata
mrna <- loaded_data$mrna
methylation <- loaded_data$methylation
rppa <- loaded_data$rppa

labels$sample_id <- as.character(labels$sample_id)
labels$subtype <- factor(labels$subtype, levels = c("LUAD", "LUSC"))

if (is.null(metadata)) {
  stop(
    "No metadata found inside 01_loaded_matched_data.rds.\n",
    "Update and rerun R/01_load_data.R so it saves matched metadata."
  )
}

metadata$sample_id <- as.character(metadata$sample_id)

# Ensure metadata has subtype
metadata <- metadata %>%
  dplyr::select(-dplyr::any_of("subtype")) %>%
  dplyr::left_join(labels, by = "sample_id") %>%
  dplyr::arrange(match(sample_id, labels$sample_id))

# Check sample order
stopifnot(identical(labels$sample_id, metadata$sample_id))
stopifnot(identical(labels$sample_id, rownames(mrna)))
stopifnot(identical(labels$sample_id, rownames(methylation)))
stopifnot(identical(labels$sample_id, rownames(rppa)))

message("Loaded matched data successfully.")
message("Samples: ", nrow(labels))
message("mRNA features: ", ncol(mrna))
message("DNA methylation features: ", ncol(methylation))
message("RPPA features: ", ncol(rppa))

# -----------------------------
# 2. Helper functions
# -----------------------------

summarize_matrix <- function(mat, layer_name) {
  feature_missing_percent <- colMeans(is.na(mat)) * 100
  sample_missing_percent <- rowMeans(is.na(mat)) * 100
  feature_var <- matrixStats::colVars(mat, na.rm = TRUE)
  
  tibble::tibble(
    layer = layer_name,
    samples = nrow(mat),
    features = ncol(mat),
    total_values = as.numeric(nrow(mat)) * as.numeric(ncol(mat)),
    missing_values = sum(is.na(mat)),
    missing_percent = round(mean(is.na(mat)) * 100, 4),
    features_with_missing = sum(feature_missing_percent > 0),
    max_feature_missing_percent = round(max(feature_missing_percent, na.rm = TRUE), 4),
    max_sample_missing_percent = round(max(sample_missing_percent, na.rm = TRUE), 4),
    zero_or_invalid_variance_features = sum(is.na(feature_var) | feature_var == 0)
  )
}

basic_qc_filter <- function(mat, layer_name, max_missing = 0.20) {
  message("\nBasic QC filtering: ", layer_name)
  
  feature_missing_prop <- colMeans(is.na(mat))
  feature_var <- matrixStats::colVars(mat, na.rm = TRUE)
  
  keep <- feature_missing_prop <= max_missing &
    !is.na(feature_var) &
    feature_var > 0
  
  filtered_mat <- mat[, keep, drop = FALSE]
  
  removed_high_missing <- sum(feature_missing_prop > max_missing, na.rm = TRUE)
  removed_zero_or_invalid_var <- sum(is.na(feature_var) | feature_var == 0)
  
  summary <- tibble::tibble(
    layer = layer_name,
    features_before = ncol(mat),
    features_after_basic_qc = ncol(filtered_mat),
    removed_features_total = ncol(mat) - ncol(filtered_mat),
    removed_percent_total = round(
      (ncol(mat) - ncol(filtered_mat)) / ncol(mat) * 100,
      4
    ),
    removed_high_missing_features = removed_high_missing,
    removed_zero_or_invalid_variance_features = removed_zero_or_invalid_var,
    max_missing_allowed_percent = max_missing * 100
  )
  
  print(summary)
  
  list(
    matrix = filtered_mat,
    summary = summary
  )
}

# -----------------------------
# 3. Summary before QC
# -----------------------------

preprocessing_summary_before <- dplyr::bind_rows(
  summarize_matrix(mrna, "mRNA"),
  summarize_matrix(methylation, "DNA methylation"),
  summarize_matrix(rppa, "RPPA")
)

message("\nSummary before basic QC:")
print(preprocessing_summary_before)

readr::write_tsv(
  preprocessing_summary_before,
  file.path(tables_dir, "preprocessing_summary_before_basic_qc.tsv")
)

# -----------------------------
# 4. Basic technical QC only
# -----------------------------
# This removes technically unusable features only.
# No imputation, no scaling, no feature selection, no batch correction.

mrna_qc <- basic_qc_filter(
  mat = mrna,
  layer_name = "mRNA",
  max_missing = 0.20
)

methylation_qc <- basic_qc_filter(
  mat = methylation,
  layer_name = "DNA methylation",
  max_missing = 0.20
)

rppa_qc <- basic_qc_filter(
  mat = rppa,
  layer_name = "RPPA",
  max_missing = 0.20
)

mrna_qc_mat <- mrna_qc$matrix
methylation_qc_mat <- methylation_qc$matrix
rppa_qc_mat <- rppa_qc$matrix

# Check sample order after feature filtering
stopifnot(identical(labels$sample_id, rownames(mrna_qc_mat)))
stopifnot(identical(labels$sample_id, rownames(methylation_qc_mat)))
stopifnot(identical(labels$sample_id, rownames(rppa_qc_mat)))

basic_qc_summary <- dplyr::bind_rows(
  mrna_qc$summary,
  methylation_qc$summary,
  rppa_qc$summary
)

readr::write_tsv(
  basic_qc_summary,
  file.path(tables_dir, "basic_qc_filtering_summary.tsv")
)

# -----------------------------
# 5. Summary after basic QC
# -----------------------------

preprocessing_summary_after_basic_qc <- dplyr::bind_rows(
  summarize_matrix(mrna_qc_mat, "mRNA"),
  summarize_matrix(methylation_qc_mat, "DNA methylation"),
  summarize_matrix(rppa_qc_mat, "RPPA")
)

message("\nSummary after basic QC:")
print(preprocessing_summary_after_basic_qc)

readr::write_tsv(
  preprocessing_summary_after_basic_qc,
  file.path(tables_dir, "preprocessing_summary_after_basic_qc.tsv")
)

# -----------------------------
# 6. Class balance
# -----------------------------

class_balance <- labels %>%
  dplyr::count(subtype, name = "n") %>%
  dplyr::mutate(percent = round(n / sum(n) * 100, 2))

message("\nClass balance:")
print(class_balance)

readr::write_tsv(
  class_balance,
  file.path(tables_dir, "class_balance_basic_qc.tsv")
)

# -----------------------------
# 7. Candidate batch/source columns for later EDA
# -----------------------------
# This only identifies variables for EDA coloring.
# It does NOT apply batch correction.

diagnostic_metadata_cols <- grep(
  "batch|plate|center|lab|platform|source|site|tss|ship|date|study|file|slide",
  colnames(metadata),
  ignore.case = TRUE,
  value = TRUE
)

# Exclude biological/label-like columns from correction candidates.
# They may still appear in metadata, but they should not be treated
# as technical batch variables for correction.
exclude_as_batch <- c(
  "sample_id",
  "subtype",
  "CANCER_TYPE",
  "CANCER_TYPE_DETAILED",
  "ONCOTREE_CODE",
  "SAMPLE_TYPE",
  "SAMPLE_CLASS",
  "TUMOR_TYPE",
  "clinical_source"
)

candidate_batch_cols <- setdiff(
  diagnostic_metadata_cols,
  exclude_as_batch
)

candidate_batch_cols <- candidate_batch_cols[
  candidate_batch_cols %in% colnames(metadata)
]

# Keep only columns with more than one non-missing value.
candidate_batch_cols <- candidate_batch_cols[
  vapply(
    candidate_batch_cols,
    function(col) length(unique(na.omit(metadata[[col]]))) > 1,
    logical(1)
  )
]

message("\nMetadata columns available:")
print(colnames(metadata))

message("\nDiagnostic metadata columns:")
print(diagnostic_metadata_cols)

message("\nCandidate batch/source columns for later EDA:")
print(candidate_batch_cols)

readr::write_tsv(
  metadata,
  file.path(tables_dir, "sample_metadata_matched_for_batch_check.tsv")
)

readr::write_tsv(
  tibble::tibble(
    diagnostic_metadata_column = diagnostic_metadata_cols
  ),
  file.path(tables_dir, "diagnostic_metadata_columns.tsv")
)

readr::write_tsv(
  tibble::tibble(
    candidate_batch_column = candidate_batch_cols
  ),
  file.path(tables_dir, "candidate_batch_columns.tsv")
)

if (length(candidate_batch_cols) > 0) {
  batch_counts <- lapply(candidate_batch_cols, function(col) {
    metadata %>%
      dplyr::mutate(batch_value = as.character(.data[[col]])) %>%
      dplyr::count(batch_value, subtype, name = "n") %>%
      dplyr::mutate(batch_variable = col) %>%
      dplyr::select(batch_variable, batch_value, subtype, n)
  }) %>%
    dplyr::bind_rows()
  
  readr::write_tsv(
    batch_counts,
    file.path(tables_dir, "batch_candidate_counts_by_subtype.tsv")
  )
  
  message("\nBatch/source candidate counts by subtype:")
  print(batch_counts)
}

# -----------------------------
# 8. Save one official QC modelling object
# -----------------------------

qc_model_input <- list(
  labels = labels,
  metadata = metadata,
  mrna = mrna_qc_mat,
  methylation = methylation_qc_mat,
  rppa = rppa_qc_mat,
  diagnostic_metadata_cols = diagnostic_metadata_cols,
  candidate_batch_cols = candidate_batch_cols,
  preprocessing_summary_before = preprocessing_summary_before,
  basic_qc_summary = basic_qc_summary,
  preprocessing_summary_after_basic_qc = preprocessing_summary_after_basic_qc,
  class_balance = class_balance,
  note = paste(
    "This object contains matched, basic-QC-filtered data only.",
    "The omics layers are stored separately as mrna, methylation, and rppa.",
    "No global imputation, scaling, top-variable filtering, supervised feature selection, or batch correction was applied.",
    "For final modelling, imputation, scaling, feature selection, and any batch correction used for model input must be fitted inside cross-validation training folds."
  )
)

qc_rds_path <- file.path(
  processed_dir,
  "02_qc_model_input_unimputed.rds"
)

saveRDS(
  qc_model_input,
  qc_rds_path,
  compress = FALSE
)

message("\nSaved leakage-aware QC modelling object to:")
message(qc_rds_path)

# -----------------------------
# 9. Save dimension table
# -----------------------------

qc_dimension_summary <- tibble::tibble(
  layer = c("mRNA", "DNA methylation", "RPPA"),
  samples = c(
    nrow(mrna_qc_mat),
    nrow(methylation_qc_mat),
    nrow(rppa_qc_mat)
  ),
  features = c(
    ncol(mrna_qc_mat),
    ncol(methylation_qc_mat),
    ncol(rppa_qc_mat)
  )
)

readr::write_tsv(
  qc_dimension_summary,
  file.path(tables_dir, "qc_dimension_summary.tsv")
)

message("\nQC dimension summary:")
print(qc_dimension_summary)

message("\n02_preprocessing.R completed successfully.")