############################################################
# 02_preprocessing.R
# QC filtering, missingness summary, imputation and feature filtering
############################################################

# -----------------------------
# 1. Setup and load data
# -----------------------------

source("R/00_setup.R")

loaded_data <- readRDS(
  file.path(processed_dir, "01_loaded_matched_data.rds")
)

labels <- loaded_data$labels
mrna <- loaded_data$mrna
methylation <- loaded_data$methylation
rppa <- loaded_data$rppa

# -----------------------------
# 2. Helper functions
# -----------------------------

summarize_matrix <- function(mat, layer_name) {
  feature_missing <- colMeans(is.na(mat)) * 100
  sample_missing <- rowMeans(is.na(mat)) * 100
  feature_var <- matrixStats::colVars(mat, na.rm = TRUE)
  
  tibble::tibble(
    layer = layer_name,
    samples = nrow(mat),
    features = ncol(mat),
    total_values = nrow(mat) * ncol(mat),
    missing_values = sum(is.na(mat)),
    missing_percent = round(mean(is.na(mat)) * 100, 4),
    features_with_missing = sum(feature_missing > 0),
    max_feature_missing_percent = round(max(feature_missing, na.rm = TRUE), 4),
    max_sample_missing_percent = round(max(sample_missing, na.rm = TRUE), 4),
    zero_or_invalid_variance_features = sum(is.na(feature_var) | feature_var == 0)
  )
}

filter_qc_features <- function(mat, layer_name, max_missing = 0.20) {
  message("\nQC filtering: ", layer_name)
  
  missing_prop <- colMeans(is.na(mat))
  feature_var <- matrixStats::colVars(mat, na.rm = TRUE)
  
  keep <- missing_prop <= max_missing &
    !is.na(feature_var) &
    feature_var > 0
  
  filtered <- mat[, keep, drop = FALSE]
  
  summary <- tibble::tibble(
    layer = layer_name,
    features_before = ncol(mat),
    features_after_qc = ncol(filtered),
    removed_features = ncol(mat) - ncol(filtered),
    removed_percent = round((ncol(mat) - ncol(filtered)) / ncol(mat) * 100, 4),
    max_missing_allowed_percent = max_missing * 100
  )
  
  print(summary)
  
  list(matrix = filtered, summary = summary)
}

median_impute <- function(mat, layer_name) {
  message("\nMedian imputation: ", layer_name)
  
  missing_before <- sum(is.na(mat))
  
  if (missing_before == 0) {
    message("No missing values found in ", layer_name)
    return(mat)
  }
  
  for (j in seq_len(ncol(mat))) {
    idx <- is.na(mat[, j])
    
    if (any(idx)) {
      med <- median(mat[, j], na.rm = TRUE)
      mat[idx, j] <- med
    }
  }
  
  missing_after <- sum(is.na(mat))
  
  message("Missing before: ", missing_before)
  message("Missing after: ", missing_after)
  
  return(mat)
}

select_top_variable_features <- function(mat, layer_name, top_n) {
  message("\nVariance filtering: ", layer_name)
  
  if (ncol(mat) <= top_n) {
    message(layer_name, " has <= ", top_n, " features. Keeping all features.")
    
    summary <- tibble::tibble(
      layer = layer_name,
      features_before_variance_filter = ncol(mat),
      features_after_variance_filter = ncol(mat),
      selected_top_n = ncol(mat)
    )
    
    return(list(matrix = mat, summary = summary))
  }
  
  vars <- matrixStats::colVars(mat, na.rm = TRUE)
  names(vars) <- colnames(mat)
  
  top_features <- names(sort(vars, decreasing = TRUE))[seq_len(top_n)]
  
  filtered <- mat[, top_features, drop = FALSE]
  
  summary <- tibble::tibble(
    layer = layer_name,
    features_before_variance_filter = ncol(mat),
    features_after_variance_filter = ncol(filtered),
    selected_top_n = top_n
  )
  
  print(summary)
  
  list(matrix = filtered, summary = summary)
}

# -----------------------------
# 3. Summary before preprocessing
# -----------------------------

preprocessing_summary_before <- dplyr::bind_rows(
  summarize_matrix(mrna, "mRNA"),
  summarize_matrix(methylation, "DNA methylation"),
  summarize_matrix(rppa, "RPPA")
)

message("\nSummary before preprocessing:")
print(preprocessing_summary_before)

readr::write_tsv(
  preprocessing_summary_before,
  file.path(tables_dir, "preprocessing_summary_before.tsv")
)

# -----------------------------
# 4. Remove high-missingness and zero-variance features
# -----------------------------
# General QC only.
# No subtype-based feature selection is done here.

mrna_qc <- filter_qc_features(mrna, "mRNA", max_missing = 0.20)
methylation_qc <- filter_qc_features(methylation, "DNA methylation", max_missing = 0.20)
rppa_qc <- filter_qc_features(rppa, "RPPA", max_missing = 0.20)

mrna_qc_mat <- mrna_qc$matrix
methylation_qc_mat <- methylation_qc$matrix
rppa_qc_mat <- rppa_qc$matrix

qc_filtering_summary <- dplyr::bind_rows(
  mrna_qc$summary,
  methylation_qc$summary,
  rppa_qc$summary
)

readr::write_tsv(
  qc_filtering_summary,
  file.path(tables_dir, "qc_filtering_summary.tsv")
)

# -----------------------------
# 5. Median imputation
# -----------------------------
# This is acceptable here for EDA and baseline modelling.
# In the final nested CV pipeline, imputation should be fitted inside training folds.

mrna_imputed <- median_impute(mrna_qc_mat, "mRNA")
methylation_imputed <- median_impute(methylation_qc_mat, "DNA methylation")
rppa_imputed <- median_impute(rppa_qc_mat, "RPPA")

# -----------------------------
# 6. Variance-based feature filtering
# -----------------------------
# Rationale:
# mRNA and methylation are very high-dimensional.
# We keep the most variable features for practical baseline modelling.
#
# Conservative but feasible:
# mRNA: top 5000
# methylation: top 10000
# RPPA: all features

mrna_var <- select_top_variable_features(
  mrna_imputed,
  "mRNA",
  top_n = 5000
)

methylation_var <- select_top_variable_features(
  methylation_imputed,
  "DNA methylation",
  top_n = 10000
)

rppa_var <- select_top_variable_features(
  rppa_imputed,
  "RPPA",
  top_n = ncol(rppa_imputed)
)

mrna_processed <- mrna_var$matrix
methylation_processed <- methylation_var$matrix
rppa_processed <- rppa_var$matrix

variance_filtering_summary <- dplyr::bind_rows(
  mrna_var$summary,
  methylation_var$summary,
  rppa_var$summary
)

readr::write_tsv(
  variance_filtering_summary,
  file.path(tables_dir, "variance_filtering_summary.tsv")
)

# -----------------------------
# 7. Summary after preprocessing
# -----------------------------

preprocessing_summary_after <- dplyr::bind_rows(
  summarize_matrix(mrna_processed, "mRNA"),
  summarize_matrix(methylation_processed, "DNA methylation"),
  summarize_matrix(rppa_processed, "RPPA")
)

message("\nSummary after preprocessing:")
print(preprocessing_summary_after)

readr::write_tsv(
  preprocessing_summary_after,
  file.path(tables_dir, "preprocessing_summary_after.tsv")
)

# -----------------------------
# 8. Create class-balance plot
# -----------------------------

class_balance <- labels %>%
  dplyr::count(subtype, name = "n") %>%
  dplyr::mutate(percent = round(n / sum(n) * 100, 2))

p_class <- ggplot2::ggplot(class_balance, ggplot2::aes(x = subtype, y = n)) +
  ggplot2::geom_col() +
  ggplot2::labs(
    title = "Class balance after sample matching",
    x = "NSCLC subtype",
    y = "Number of samples"
  ) +
  ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(
  file.path(figures_dir, "class_balance.png"),
  p_class,
  width = 5,
  height = 4,
  dpi = 300
)

# -----------------------------
# 9. Create feature-count plot
# -----------------------------

feature_counts <- tibble::tibble(
  layer = rep(c("mRNA", "DNA methylation", "RPPA"), each = 3),
  stage = rep(c("Original", "After QC", "After variance filtering"), times = 3),
  features = c(
    ncol(mrna), ncol(mrna_qc_mat), ncol(mrna_processed),
    ncol(methylation), ncol(methylation_qc_mat), ncol(methylation_processed),
    ncol(rppa), ncol(rppa_qc_mat), ncol(rppa_processed)
  )
)

readr::write_tsv(
  feature_counts,
  file.path(tables_dir, "feature_counts_by_stage.tsv")
)

p_features <- ggplot2::ggplot(
  feature_counts,
  ggplot2::aes(x = stage, y = features, fill = layer)
) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::facet_wrap(~ layer, scales = "free_y") +
  ggplot2::labs(
    title = "Feature counts across preprocessing steps",
    x = "Preprocessing stage",
    y = "Number of features"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
  )

ggplot2::ggsave(
  file.path(figures_dir, "feature_counts_by_stage.png"),
  p_features,
  width = 9,
  height = 5,
  dpi = 300
)

# -----------------------------
# 10. Save processed data
# -----------------------------

processed_data <- list(
  labels = labels,
  mrna = mrna_processed,
  methylation = methylation_processed,
  rppa = rppa_processed,
  preprocessing_summary_before = preprocessing_summary_before,
  qc_filtering_summary = qc_filtering_summary,
  variance_filtering_summary = variance_filtering_summary,
  preprocessing_summary_after = preprocessing_summary_after,
  feature_counts = feature_counts
)

saveRDS(
  processed_data,
  file.path(processed_dir, "02_preprocessed_data.rds"),
  compress = FALSE
)

message("\nSaved preprocessed data to:")
message(file.path(processed_dir, "02_preprocessed_data.rds"))

message("\n02_preprocessing.R completed successfully.")