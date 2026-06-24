############################################################
# 03_eda.R
# Exploratory data analysis for preprocessed NSCLC multi-omics data
############################################################

# -----------------------------
# 1. Setup and load preprocessed data
# -----------------------------

source("R/00_setup.R")

processed_data <- readRDS(
  file.path(processed_dir, "02_preprocessed_data.rds")
)

labels <- processed_data$labels
mrna <- processed_data$mrna
methylation <- processed_data$methylation
rppa <- processed_data$rppa

# Final sample-order checks
stopifnot(identical(labels$sample_id, rownames(mrna)))
stopifnot(identical(labels$sample_id, rownames(methylation)))
stopifnot(identical(labels$sample_id, rownames(rppa)))

message("EDA data loaded successfully.")

# -----------------------------
# 2. Plot colors
# -----------------------------

subtype_colors <- c(
  "LUAD" = "#1F77B4",  # blue
  "LUSC" = "#D55E00"   # orange
)

layer_colors <- c(
  "mRNA" = "#1F77B4",
  "DNA methylation" = "#009E73",
  "RPPA" = "#CC79A7"
)

# -----------------------------
# 3. Class balance plot
# -----------------------------

class_balance <- labels %>%
  dplyr::count(subtype, name = "n") %>%
  dplyr::mutate(percent = round(n / sum(n) * 100, 2))

readr::write_tsv(
  class_balance,
  file.path(tables_dir, "eda_class_balance.tsv")
)

p_class <- ggplot2::ggplot(
  class_balance,
  ggplot2::aes(x = subtype, y = n, fill = subtype)
) +
  ggplot2::geom_col(width = 0.65, alpha = 0.85) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0(n, " (", percent, "%)")),
    vjust = -0.4,
    size = 4
  ) +
  ggplot2::scale_fill_manual(values = subtype_colors) +
  ggplot2::labs(
    title = "Class balance after sample matching",
    x = "NSCLC subtype",
    y = "Number of samples",
    fill = "Subtype"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "none"
  )

ggplot2::ggsave(
  file.path(figures_dir, "eda_class_balance.png"),
  p_class,
  width = 5.5,
  height = 4.2,
  dpi = 300
)

# -----------------------------
# 4. Feature-count plot
# -----------------------------

if ("feature_counts" %in% names(processed_data)) {
  
  feature_counts <- processed_data$feature_counts
  
  readr::write_tsv(
    feature_counts,
    file.path(tables_dir, "eda_feature_counts_by_stage.tsv")
  )
  
  p_features <- ggplot2::ggplot(
    feature_counts,
    ggplot2::aes(x = stage, y = features, fill = layer)
  ) +
    ggplot2::geom_col(width = 0.65, alpha = 0.85) +
    ggplot2::facet_wrap(~ layer, scales = "free_y") +
    ggplot2::scale_fill_manual(values = layer_colors) +
    ggplot2::labs(
      title = "Feature counts across preprocessing steps",
      x = "Preprocessing stage",
      y = "Number of features",
      fill = "Omics layer"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position = "none"
    )
  
  ggplot2::ggsave(
    file.path(figures_dir, "eda_feature_counts_by_stage.png"),
    p_features,
    width = 9,
    height = 5,
    dpi = 300
  )
}

# -----------------------------
# 5. Sample-level summary function
# -----------------------------

get_sample_summary <- function(mat, layer_name) {
  tibble::tibble(
    sample_id = rownames(mat),
    layer = layer_name,
    sample_mean = rowMeans(mat, na.rm = TRUE),
    sample_sd = matrixStats::rowSds(mat, na.rm = TRUE),
    sample_missing_percent = rowMeans(is.na(mat)) * 100
  ) %>%
    dplyr::left_join(labels, by = "sample_id")
}

sample_summary <- dplyr::bind_rows(
  get_sample_summary(mrna, "mRNA"),
  get_sample_summary(methylation, "DNA methylation"),
  get_sample_summary(rppa, "RPPA")
)

readr::write_tsv(
  sample_summary,
  file.path(tables_dir, "eda_sample_level_summary.tsv")
)

message("\nSample-level summary:")
print(head(sample_summary))

# -----------------------------
# 6. Sample mean distribution plot
# -----------------------------

p_sample_mean <- ggplot2::ggplot(
  sample_summary,
  ggplot2::aes(x = subtype, y = sample_mean, fill = subtype)
) +
  ggplot2::geom_boxplot(outlier.alpha = 0.35, alpha = 0.75) +
  ggplot2::facet_wrap(~ layer, scales = "free_y") +
  ggplot2::scale_fill_manual(values = subtype_colors) +
  ggplot2::labs(
    title = "Sample-level mean signal by omics layer",
    x = "NSCLC subtype",
    y = "Mean feature value per sample",
    fill = "Subtype"
  ) +
  ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(
  file.path(figures_dir, "eda_sample_mean_distribution_by_layer.png"),
  p_sample_mean,
  width = 9,
  height = 5,
  dpi = 300
)

# -----------------------------
# 7. PCA helper function
# -----------------------------
# Important correction:
# Variance explained is calculated using total variance of the full scaled matrix,
# not only the variance of the PCs returned by prcomp(rank.).

run_pca_eda <- function(mat, layer_name, n_pcs_to_show = 20) {
  message("\nRunning PCA for ", layer_name)
  
  # Scale each feature before PCA
  x <- scale(mat)
  x <- as.matrix(x)
  
  # Remove invalid columns if any remain after scaling
  valid_cols <- colSums(is.na(x)) == 0 &
    colSums(is.infinite(x)) == 0
  
  x <- x[, valid_cols, drop = FALSE]
  
  message(layer_name, " PCA input:")
  message(nrow(x), " samples x ", ncol(x), " valid features")
  
  # Total variance of the full scaled matrix
  total_variance <- sum(matrixStats::colVars(x, na.rm = TRUE))
  
  # Safe number of PCs
  n_pcs <- min(n_pcs_to_show, nrow(x) - 1, ncol(x))
  
  pca <- prcomp(
    x,
    center = FALSE,
    scale. = FALSE,
    rank. = n_pcs
  )
  
  # Correct PVE denominator
  pve <- (pca$sdev^2) / total_variance
  
  pca_df <- tibble::tibble(
    sample_id = rownames(x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3]
  ) %>%
    dplyr::left_join(labels, by = "sample_id")
  
  variance_df <- tibble::tibble(
    layer = layer_name,
    PC = paste0("PC", seq_along(pve)),
    PC_number = seq_along(pve),
    variance_explained = pve,
    variance_explained_percent = round(pve * 100, 3)
  )
  
  clean_layer_name <- gsub(" ", "_", layer_name)
  
  readr::write_tsv(
    pca_df,
    file.path(tables_dir, paste0("eda_pca_coordinates_", clean_layer_name, ".tsv"))
  )
  
  readr::write_tsv(
    variance_df,
    file.path(tables_dir, paste0("eda_pca_variance_", clean_layer_name, ".tsv"))
  )
  
  p_pca <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = PC1, y = PC2, color = subtype, shape = subtype)
  ) +
    ggplot2::geom_point(size = 2.5, alpha = 0.85) +
    ggplot2::scale_color_manual(values = subtype_colors) +
    ggplot2::labs(
      title = paste0("PCA of ", layer_name),
      x = paste0("PC1: ", round(pve[1] * 100, 2), "% variance"),
      y = paste0("PC2: ", round(pve[2] * 100, 2), "% variance"),
      color = "Subtype",
      shape = "Subtype"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  
  ggplot2::ggsave(
    file.path(figures_dir, paste0("eda_pca_", clean_layer_name, ".png")),
    p_pca,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  return(
    list(
      pca = pca,
      pca_df = pca_df,
      variance_df = variance_df,
      plot = p_pca
    )
  )
}

# -----------------------------
# 8. Run PCA for each omics layer
# -----------------------------

pca_mrna <- run_pca_eda(
  mrna,
  "mRNA",
  n_pcs_to_show = 20
)

pca_methylation <- run_pca_eda(
  methylation,
  "DNA methylation",
  n_pcs_to_show = 20
)

pca_rppa <- run_pca_eda(
  rppa,
  "RPPA",
  n_pcs_to_show = 20
)

# -----------------------------
# 9. Combine and save PCA variance table
# -----------------------------

pca_variance <- dplyr::bind_rows(
  pca_mrna$variance_df,
  pca_methylation$variance_df,
  pca_rppa$variance_df
)

readr::write_tsv(
  pca_variance,
  file.path(tables_dir, "eda_pca_variance_explained.tsv")
)

message("\nPCA variance explained:")
print(pca_variance)

# -----------------------------
# 10. PCA variance explained plot
# -----------------------------

p_pve <- ggplot2::ggplot(
  pca_variance,
  ggplot2::aes(
    x = PC_number,
    y = variance_explained_percent,
    color = layer,
    group = layer
  )
) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_color_manual(values = layer_colors) +
  ggplot2::facet_wrap(~ layer, scales = "free_y") +
  ggplot2::labs(
    title = "Variance explained by first principal components",
    x = "Principal component",
    y = "Variance explained (%)",
    color = "Omics layer"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "none"
  )

ggplot2::ggsave(
  file.path(figures_dir, "eda_pca_variance_explained.png"),
  p_pve,
  width = 9,
  height = 5,
  dpi = 300
)

# -----------------------------
# 11. Save EDA object
# -----------------------------

eda_results <- list(
  class_balance = class_balance,
  sample_summary = sample_summary,
  pca_mrna = pca_mrna,
  pca_methylation = pca_methylation,
  pca_rppa = pca_rppa,
  pca_variance = pca_variance
)

saveRDS(
  eda_results,
  file.path(processed_dir, "03_eda_results.rds"),
  compress = FALSE
)

# -----------------------------
# 12. Print plots in RStudio
# -----------------------------

print(p_class)

if (exists("p_features")) {
  print(p_features)
}

print(p_sample_mean)
print(pca_mrna$plot)
print(pca_methylation$plot)
print(pca_rppa$plot)
print(p_pve)

message("\nSaved EDA results to:")
message(file.path(processed_dir, "03_eda_results.rds"))

message("\n03_eda.R completed successfully.")