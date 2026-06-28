############################################################
# 03_eda.R
# Exploratory data analysis and source-site/batch inspection
# NSCLC multi-omics project: LUAD vs LUSC
#
# Input:
#   data/processed/02_qc_model_input_unimputed.rds
#
# Output:
#   results/figures/eda/*.png
#   results/tables/*.tsv
#   data/processed/03_eda_results.rds
#
# This script does:
#   - basic visual EDA
#   - PCA for each omics layer
#   - PCA colored by subtype
#   - PCA colored by metadata/source-site variables
#   - source-site/subtype confounding analysis
#   - visual explanation for why source-site correction is not applied
#
# This script does NOT do:
#   - model training
#   - batch correction
#   - save globally imputed/scaled data for modelling
#
# Important:
# PCA uses temporary imputation, top-variable selection and scaling
# only for visualization.
############################################################

# -----------------------------
# 1. Setup
# -----------------------------

source("R/00_setup.R")

qc_rds_path <- file.path(processed_dir, "02_qc_model_input_unimputed.rds")

if (!file.exists(qc_rds_path)) {
  stop(
    "02_qc_model_input_unimputed.rds not found.\n",
    "Run source('R/02_preprocessing.R') first."
  )
}

qc_data <- readRDS(qc_rds_path)

labels <- qc_data$labels
metadata <- qc_data$metadata
mrna <- qc_data$mrna
methylation <- qc_data$methylation
rppa <- qc_data$rppa

labels$sample_id <- as.character(labels$sample_id)
labels$subtype <- factor(labels$subtype, levels = c("LUAD", "LUSC"))

metadata$sample_id <- as.character(metadata$sample_id)

metadata <- metadata %>%
  dplyr::select(-dplyr::any_of("subtype")) %>%
  dplyr::left_join(labels, by = "sample_id") %>%
  dplyr::arrange(match(sample_id, labels$sample_id))

stopifnot(identical(labels$sample_id, metadata$sample_id))
stopifnot(identical(labels$sample_id, rownames(mrna)))
stopifnot(identical(labels$sample_id, rownames(methylation)))
stopifnot(identical(labels$sample_id, rownames(rppa)))

figures_eda_dir <- file.path(figures_dir, "eda")
dir.create(figures_eda_dir, recursive = TRUE, showWarnings = FALSE)

message("Loaded QC data for EDA.")
message("Samples: ", nrow(labels))
message("mRNA features: ", ncol(mrna))
message("DNA methylation features: ", ncol(methylation))
message("RPPA features: ", ncol(rppa))

# -----------------------------
# 2. Plot theme and colors
# -----------------------------

subtype_colors <- c(
  "LUAD" = "#1F77B4",
  "LUSC" = "#D62728"
)

omics_colors <- c(
  "mRNA" = "#1F77B4",
  "DNA methylation" = "#2CA02C",
  "RPPA" = "#9467BD"
)

qc_colors <- c(
  "Kept" = "#2CA02C",
  "Removed" = "#D62728"
)

theme_eda <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 3),
      plot.subtitle = ggplot2::element_text(size = base_size, color = "grey30"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "grey20"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      plot.caption = ggplot2::element_text(color = "grey35", size = base_size - 2)
    )
}

save_plot <- function(plot, filename, width = 8, height = 6, dpi = 320) {
  ggplot2::ggsave(
    filename = file.path(figures_eda_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

# -----------------------------
# 3. Basic EDA summary tables
# -----------------------------

class_balance <- labels %>%
  dplyr::count(subtype, name = "n") %>%
  dplyr::mutate(percent = round(n / sum(n) * 100, 2))

dimension_summary <- tibble::tibble(
  layer = c("mRNA", "DNA methylation", "RPPA"),
  samples = c(nrow(mrna), nrow(methylation), nrow(rppa)),
  features = c(ncol(mrna), ncol(methylation), ncol(rppa))
)

basic_qc_summary <- qc_data$basic_qc_summary
missing_summary_after_qc <- qc_data$preprocessing_summary_after_basic_qc

readr::write_tsv(class_balance, file.path(tables_dir, "eda_class_balance.tsv"))
readr::write_tsv(dimension_summary, file.path(tables_dir, "eda_dimension_summary_after_qc.tsv"))
readr::write_tsv(missing_summary_after_qc, file.path(tables_dir, "eda_missing_summary_after_qc.tsv"))

# -----------------------------
# 4. Basic EDA visualizations
# -----------------------------

# 4.1 Class balance
p_class <- ggplot2::ggplot(
  class_balance,
  ggplot2::aes(x = subtype, y = n, fill = subtype)
) +
  ggplot2::geom_col(width = 0.65, alpha = 0.92) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0(n, " samples\n", percent, "%")),
    vjust = -0.25,
    fontface = "bold",
    size = 4.2
  ) +
  ggplot2::scale_fill_manual(values = subtype_colors) +
  ggplot2::ylim(0, max(class_balance$n) * 1.18) +
  ggplot2::labs(
    title = "Class balance after sample matching",
    subtitle = "Matched LUAD and LUSC samples used across all omics layers",
    x = "Subtype",
    y = "Number of samples",
    caption = "Samples are matched across mRNA, DNA methylation and RPPA."
  ) +
  theme_eda() +
  ggplot2::theme(legend.position = "none")

save_plot(p_class, "01_class_balance.png", width = 7, height = 5)

# 4.2 Feature counts after QC
p_features <- ggplot2::ggplot(
  dimension_summary,
  ggplot2::aes(x = reorder(layer, features), y = features, fill = layer)
) +
  ggplot2::geom_col(width = 0.7, alpha = 0.92) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = omics_colors) +
  ggplot2::scale_y_continuous(labels = scales::comma) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::comma(features)),
    hjust = -0.08,
    fontface = "bold",
    size = 4
  ) +
  ggplot2::labs(
    title = "Feature counts after basic QC",
    subtitle = "Each omics layer is kept separate for downstream modelling",
    x = "Omics layer",
    y = "Number of features",
    caption = "Only high-missingness and zero/invalid variance features were removed."
  ) +
  theme_eda() +
  ggplot2::theme(legend.position = "none")

save_plot(p_features, "02_feature_counts_after_qc.png", width = 8.5, height = 5)

# 4.3 Kept vs removed features
qc_long <- basic_qc_summary %>%
  dplyr::transmute(
    layer,
    Kept = features_after_basic_qc,
    Removed = removed_features_total
  ) %>%
  tidyr::pivot_longer(
    cols = c("Kept", "Removed"),
    names_to = "status",
    values_to = "features"
  )

p_qc <- ggplot2::ggplot(
  qc_long,
  ggplot2::aes(x = layer, y = features, fill = status)
) +
  ggplot2::geom_col(width = 0.7, alpha = 0.92) +
  ggplot2::scale_fill_manual(values = qc_colors) +
  ggplot2::scale_y_continuous(labels = scales::comma) +
  ggplot2::labs(
    title = "Basic QC retained most features",
    subtitle = "Removed features had high missingness or zero/invalid variance",
    x = "Omics layer",
    y = "Number of features",
    fill = "QC status"
  ) +
  theme_eda()

save_plot(p_qc, "03_qc_kept_vs_removed_features.png", width = 8, height = 5)

# 4.4 Missingness after QC
missing_plot_df <- missing_summary_after_qc %>%
  dplyr::select(layer, missing_percent, features_with_missing, max_feature_missing_percent)

p_missing <- ggplot2::ggplot(
  missing_plot_df,
  ggplot2::aes(x = layer, y = missing_percent, fill = layer)
) +
  ggplot2::geom_col(width = 0.65, alpha = 0.92) +
  ggplot2::scale_fill_manual(values = omics_colors) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0(missing_percent, "%")),
    vjust = -0.35,
    fontface = "bold",
    size = 4
  ) +
  ggplot2::ylim(0, max(missing_plot_df$missing_percent) * 1.4 + 0.01) +
  ggplot2::labs(
    title = "Missingness after basic QC",
    subtitle = "Remaining methylation missingness is low and will be imputed inside CV folds",
    x = "Omics layer",
    y = "Missing values (%)",
    caption = "No global imputation was saved as modelling input."
  ) +
  theme_eda() +
  ggplot2::theme(legend.position = "none")

save_plot(p_missing, "04_missingness_after_qc.png", width = 8, height = 5)

# -----------------------------
# 5. Sample-level missingness
# -----------------------------

sample_missing_df <- dplyr::bind_rows(
  tibble::tibble(
    sample_id = rownames(mrna),
    layer = "mRNA",
    missing_percent = rowMeans(is.na(mrna)) * 100
  ),
  tibble::tibble(
    sample_id = rownames(methylation),
    layer = "DNA methylation",
    missing_percent = rowMeans(is.na(methylation)) * 100
  ),
  tibble::tibble(
    sample_id = rownames(rppa),
    layer = "RPPA",
    missing_percent = rowMeans(is.na(rppa)) * 100
  )
) %>%
  dplyr::left_join(labels, by = "sample_id")

readr::write_tsv(
  sample_missing_df,
  file.path(tables_dir, "eda_sample_level_missingness.tsv")
)

p_sample_missing <- ggplot2::ggplot(
  sample_missing_df,
  ggplot2::aes(x = layer, y = missing_percent, fill = subtype)
) +
  ggplot2::geom_boxplot(alpha = 0.78, outlier.alpha = 0.55) +
  ggplot2::scale_fill_manual(values = subtype_colors) +
  ggplot2::labs(
    title = "Sample-level missingness after basic QC",
    subtitle = "Checks whether one subtype has systematically higher missingness",
    x = "Omics layer",
    y = "Missing values per sample (%)",
    fill = "Subtype"
  ) +
  theme_eda()

save_plot(p_sample_missing, "05_sample_level_missingness_by_subtype.png", width = 8, height = 5)

# -----------------------------
# 6. Sampled value distributions
# -----------------------------

sample_values_for_distribution <- function(mat, layer_name, n_values = 50000, seed = 42) {
  set.seed(seed)
  
  total_values <- length(mat)
  n_take <- min(n_values, total_values)
  
  sampled_values <- as.numeric(mat[sample.int(total_values, n_take)])
  sampled_values <- sampled_values[is.finite(sampled_values)]
  
  tibble::tibble(
    layer = layer_name,
    value = sampled_values
  )
}

value_distribution_df <- dplyr::bind_rows(
  sample_values_for_distribution(mrna, "mRNA", n_values = 50000, seed = 1),
  sample_values_for_distribution(methylation, "DNA methylation", n_values = 50000, seed = 2),
  sample_values_for_distribution(rppa, "RPPA", n_values = 50000, seed = 3)
)

readr::write_tsv(
  value_distribution_df,
  file.path(tables_dir, "eda_sampled_value_distribution.tsv")
)

p_value_dist <- ggplot2::ggplot(
  value_distribution_df,
  ggplot2::aes(x = value, fill = layer)
) +
  ggplot2::geom_density(alpha = 0.45, linewidth = 0.75) +
  ggplot2::facet_wrap(~layer, scales = "free", ncol = 1) +
  ggplot2::scale_fill_manual(values = omics_colors) +
  ggplot2::labs(
    title = "Sampled value distributions by omics layer",
    subtitle = "Different omics layers have different value ranges",
    x = "Feature value",
    y = "Density",
    fill = "Omics layer",
    caption = "This motivates fold-wise scaling before modelling."
  ) +
  theme_eda() +
  ggplot2::theme(legend.position = "none")

save_plot(p_value_dist, "06_value_distributions_by_layer.png", width = 8, height = 8)

# -----------------------------
# 7. PCA helper functions
# -----------------------------

prepare_matrix_for_pca <- function(mat, layer_name, top_n = 2000) {
  message("\nPreparing PCA matrix for ", layer_name)
  
  feature_var <- matrixStats::colVars(mat, na.rm = TRUE)
  names(feature_var) <- colnames(mat)
  
  keep <- !is.na(feature_var) & feature_var > 0
  
  mat <- mat[, keep, drop = FALSE]
  feature_var <- feature_var[keep]
  
  n_select <- min(top_n, ncol(mat))
  
  selected_features <- names(sort(feature_var, decreasing = TRUE))[seq_len(n_select)]
  
  x <- mat[, selected_features, drop = FALSE]
  
  feature_medians <- matrixStats::colMedians(x, na.rm = TRUE)
  names(feature_medians) <- colnames(x)
  
  missing_index <- which(is.na(x), arr.ind = TRUE)
  
  if (nrow(missing_index) > 0) {
    x[missing_index] <- feature_medians[missing_index[, 2]]
  }
  
  x_scaled <- scale(x, center = TRUE, scale = TRUE)
  
  valid_scaled <- apply(x_scaled, 2, function(z) all(is.finite(z)))
  
  x_scaled <- x_scaled[, valid_scaled, drop = FALSE]
  
  message(layer_name, " PCA features used: ", ncol(x_scaled))
  
  list(
    x_scaled = x_scaled,
    selected_features = colnames(x_scaled)
  )
}

run_pca <- function(x_scaled, layer_name) {
  pca <- prcomp(
    x_scaled,
    center = FALSE,
    scale. = FALSE
  )
  
  percent_var <- round((pca$sdev^2 / sum(pca$sdev^2)) * 100, 2)
  
  scores <- as.data.frame(pca$x[, 1:5, drop = FALSE]) %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::mutate(layer = layer_name)
  
  list(
    pca = pca,
    scores = scores,
    percent_var = percent_var
  )
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

# -----------------------------
# 8. PCA by subtype
# -----------------------------

plot_pca_by_subtype <- function(pca_scores, percent_var, layer_name) {
  p <- pca_scores %>%
    dplyr::left_join(labels, by = "sample_id") %>%
    ggplot2::ggplot(
      ggplot2::aes(x = PC1, y = PC2, color = subtype)
    ) +
    ggplot2::geom_point(size = 2.8, alpha = 0.82) +
    ggplot2::stat_ellipse(
      ggplot2::aes(group = subtype),
      linewidth = 0.9,
      alpha = 0.7,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = subtype_colors) +
    ggplot2::labs(
      title = paste0(layer_name, " PCA colored by subtype"),
      subtitle = "PCA uses temporary EDA-only imputation, top-variable filtering and scaling",
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)"),
      color = "Subtype"
    ) +
    theme_eda()
  
  save_plot(
    p,
    paste0("07_pca_subtype_", safe_filename(layer_name), ".png"),
    width = 8,
    height = 6
  )
  
  p
}

# -----------------------------
# 9. PCA by metadata/source site
# -----------------------------

make_top_category <- function(x, top_n = 12, other_label = "Other / small sites") {
  x <- as.character(x)
  
  counts <- sort(table(x), decreasing = TRUE)
  top_levels <- names(counts)[seq_len(min(top_n, length(counts)))]
  
  out <- ifelse(x %in% top_levels, x, other_label)
  
  factor(out, levels = c(top_levels, other_label))
}

plot_pca_by_metadata <- function(pca_scores,
                                 percent_var,
                                 layer_name,
                                 metadata_col,
                                 top_n = 12) {
  plot_df <- pca_scores %>%
    dplyr::left_join(metadata, by = "sample_id") %>%
    dplyr::mutate(
      metadata_group = make_top_category(.data[[metadata_col]], top_n = top_n)
    )
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = metadata_group,
      shape = subtype
    )
  ) +
    ggplot2::geom_point(size = 2.6, alpha = 0.85) +
    ggplot2::scale_color_viridis_d(option = "turbo", end = 0.92) +
    ggplot2::labs(
      title = paste0(layer_name, " PCA colored by ", metadata_col),
      subtitle = paste0(
        "Top ", top_n,
        " source groups shown; smaller groups collapsed to Other"
      ),
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)"),
      color = metadata_col,
      shape = "Subtype",
      caption = "Used only for source-site/batch inspection. No correction is applied here."
    ) +
    theme_eda(base_size = 12) +
    ggplot2::theme(
      legend.key.height = ggplot2::unit(0.45, "cm")
    )
  
  save_plot(
    p,
    paste0(
      "08_pca_metadata_",
      safe_filename(layer_name),
      "_",
      safe_filename(metadata_col),
      ".png"
    ),
    width = 10,
    height = 7
  )
  
  p
}

# -----------------------------
# 10. Run PCA
# -----------------------------

pca_feature_plan <- tibble::tibble(
  layer = c("mRNA", "DNA methylation", "RPPA"),
  top_features_for_pca = c(2000, 5000, ncol(rppa))
)

readr::write_tsv(
  pca_feature_plan,
  file.path(tables_dir, "eda_pca_feature_plan.tsv")
)

pca_mrna_input <- prepare_matrix_for_pca(
  mrna,
  layer_name = "mRNA",
  top_n = 2000
)

pca_methylation_input <- prepare_matrix_for_pca(
  methylation,
  layer_name = "DNA methylation",
  top_n = 5000
)

pca_rppa_input <- prepare_matrix_for_pca(
  rppa,
  layer_name = "RPPA",
  top_n = ncol(rppa)
)

pca_mrna <- run_pca(pca_mrna_input$x_scaled, "mRNA")
pca_methylation <- run_pca(pca_methylation_input$x_scaled, "DNA methylation")
pca_rppa <- run_pca(pca_rppa_input$x_scaled, "RPPA")

pca_scores_all <- dplyr::bind_rows(
  pca_mrna$scores,
  pca_methylation$scores,
  pca_rppa$scores
)

pca_variance_summary <- dplyr::bind_rows(
  tibble::tibble(
    layer = "mRNA",
    PC = paste0("PC", seq_along(pca_mrna$percent_var)),
    percent_variance = pca_mrna$percent_var
  ),
  tibble::tibble(
    layer = "DNA methylation",
    PC = paste0("PC", seq_along(pca_methylation$percent_var)),
    percent_variance = pca_methylation$percent_var
  ),
  tibble::tibble(
    layer = "RPPA",
    PC = paste0("PC", seq_along(pca_rppa$percent_var)),
    percent_variance = pca_rppa$percent_var
  )
)

readr::write_tsv(
  pca_scores_all,
  file.path(tables_dir, "eda_pca_scores.tsv")
)

readr::write_tsv(
  pca_variance_summary,
  file.path(tables_dir, "eda_pca_variance_summary.tsv")
)

plot_pca_by_subtype(pca_mrna$scores, pca_mrna$percent_var, "mRNA")
plot_pca_by_subtype(pca_methylation$scores, pca_methylation$percent_var, "DNA methylation")
plot_pca_by_subtype(pca_rppa$scores, pca_rppa$percent_var, "RPPA")

metadata_cols_for_eda <- c(
  "TCGA_TSS",
  "TISSUE_SOURCE_SITE"
)

metadata_cols_for_eda <- metadata_cols_for_eda[
  metadata_cols_for_eda %in% colnames(metadata)
]

for (col in metadata_cols_for_eda) {
  plot_pca_by_metadata(
    pca_mrna$scores,
    pca_mrna$percent_var,
    "mRNA",
    metadata_col = col,
    top_n = 12
  )
  
  plot_pca_by_metadata(
    pca_methylation$scores,
    pca_methylation$percent_var,
    "DNA methylation",
    metadata_col = col,
    top_n = 12
  )
  
  plot_pca_by_metadata(
    pca_rppa$scores,
    pca_rppa$percent_var,
    "RPPA",
    metadata_col = col,
    top_n = 12
  )
}

# -----------------------------
# 11. Source-site/subtype confounding analysis
# -----------------------------

check_confounding <- function(data, batch_col) {
  out <- data %>%
    dplyr::mutate(batch_value = as.character(.data[[batch_col]])) %>%
    dplyr::count(batch_value, subtype, name = "n") %>%
    tidyr::pivot_wider(
      names_from = subtype,
      values_from = n,
      values_fill = 0
    )
  
  if (!"LUAD" %in% colnames(out)) {
    out$LUAD <- 0
  }
  
  if (!"LUSC" %in% colnames(out)) {
    out$LUSC <- 0
  }
  
  out %>%
    dplyr::mutate(
      total = LUAD + LUSC,
      dominant_subtype = dplyr::case_when(
        LUAD > LUSC ~ "LUAD",
        LUSC > LUAD ~ "LUSC",
        TRUE ~ "Mixed equal"
      ),
      dominant_fraction = pmax(LUAD, LUSC) / total,
      only_one_subtype = LUAD == 0 | LUSC == 0,
      batch_variable = batch_col
    ) %>%
    dplyr::select(
      batch_variable,
      batch_value,
      LUAD,
      LUSC,
      total,
      dominant_subtype,
      dominant_fraction,
      only_one_subtype
    ) %>%
    dplyr::arrange(dplyr::desc(total))
}

confounding_tables <- lapply(metadata_cols_for_eda, function(col) {
  check_confounding(metadata, col)
}) %>%
  dplyr::bind_rows()

confounding_overview <- confounding_tables %>%
  dplyr::group_by(batch_variable) %>%
  dplyr::summarise(
    n_groups = dplyr::n(),
    groups_with_only_one_subtype = sum(only_one_subtype),
    percent_only_one_subtype = round(groups_with_only_one_subtype / n_groups * 100, 2),
    median_dominant_fraction = round(median(dominant_fraction), 3),
    min_dominant_fraction = round(min(dominant_fraction), 3),
    .groups = "drop"
  )

readr::write_tsv(
  confounding_tables,
  file.path(tables_dir, "eda_source_site_confounding_summary.tsv")
)

readr::write_tsv(
  confounding_overview,
  file.path(tables_dir, "eda_source_site_confounding_overview.tsv")
)

message("\nSource-site confounding overview:")
print(confounding_overview)

# -----------------------------
# 12. Visual confounding plots
# -----------------------------

# 12.1 TCGA_TSS stacked bar
if ("TCGA_TSS" %in% metadata_cols_for_eda) {
  tss_stacked <- metadata %>%
    dplyr::count(TCGA_TSS, subtype, name = "n") %>%
    dplyr::group_by(TCGA_TSS) %>%
    dplyr::mutate(total = sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      TCGA_TSS = stats::reorder(TCGA_TSS, total)
    )
  
  p_tss_stacked <- ggplot2::ggplot(
    tss_stacked,
    ggplot2::aes(x = TCGA_TSS, y = n, fill = subtype)
  ) +
    ggplot2::geom_col(width = 0.85, alpha = 0.94) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = subtype_colors) +
    ggplot2::labs(
      title = "TCGA tissue-source-site is confounded with subtype",
      subtitle = "Each TCGA_TSS group contains only LUAD or only LUSC samples",
      x = "TCGA tissue-source-site code",
      y = "Number of samples",
      fill = "Subtype",
      caption = "Because source site and subtype are not separable, site-based correction could remove biological signal."
    ) +
    theme_eda(base_size = 11)
  
  save_plot(
    p_tss_stacked,
    "09_confounding_TCGA_TSS_stacked_by_subtype.png",
    width = 9,
    height = 12
  )
}

# 12.2 Confounding overview plot
p_confounding_overview <- ggplot2::ggplot(
  confounding_overview,
  ggplot2::aes(
    x = batch_variable,
    y = percent_only_one_subtype,
    fill = batch_variable
  )
) +
  ggplot2::geom_col(width = 0.65, alpha = 0.92) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0(
        groups_with_only_one_subtype,
        "/",
        n_groups,
        " groups"
      )
    ),
    vjust = -0.35,
    fontface = "bold",
    size = 4.2
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 110),
    labels = function(x) paste0(x, "%")
  ) +
  ggplot2::scale_fill_viridis_d(option = "plasma", end = 0.85) +
  ggplot2::labs(
    title = "Source-site groups are not balanced across subtypes",
    subtitle = "High values indicate confounding between source/site and LUAD/LUSC label",
    x = "Metadata/source variable",
    y = "Groups containing only one subtype",
    fill = "Variable",
    caption = "This supports inspection of source effects, but not direct correction by source site."
  ) +
  theme_eda() +
  ggplot2::theme(legend.position = "none")

save_plot(
  p_confounding_overview,
  "10_confounding_overview_single_subtype_groups.png",
  width = 9,
  height = 5.5
)

# 12.3 Dominant subtype fraction
dominant_fraction_df <- confounding_tables %>%
  dplyr::mutate(
    dominant_fraction_percent = dominant_fraction * 100
  )

p_dominant_fraction <- ggplot2::ggplot(
  dominant_fraction_df,
  ggplot2::aes(
    x = dominant_fraction_percent,
    fill = dominant_subtype
  )
) +
  ggplot2::geom_histogram(
    binwidth = 5,
    boundary = 0,
    alpha = 0.88,
    color = "white"
  ) +
  ggplot2::facet_wrap(~batch_variable, ncol = 1) +
  ggplot2::scale_fill_manual(
    values = c(
      "LUAD" = subtype_colors["LUAD"],
      "LUSC" = subtype_colors["LUSC"],
      "Mixed equal" = "grey50"
    )
  ) +
  ggplot2::labs(
    title = "Dominant subtype fraction within each source group",
    subtitle = "Values near 100% mean a source group contains almost only one subtype",
    x = "Dominant subtype fraction within source group (%)",
    y = "Number of source groups",
    fill = "Dominant subtype",
    caption = "When source groups are subtype-specific, correction cannot distinguish source effects from biology."
  ) +
  theme_eda()

save_plot(
  p_dominant_fraction,
  "11_confounding_dominant_subtype_fraction.png",
  width = 9,
  height = 7
)

# -----------------------------
# 13. Optional mixed-source-site subset check
# -----------------------------
# This checks whether broad TISSUE_SOURCE_SITE names contain both LUAD and LUSC.
# It does NOT change the main dataset.
# It only creates a summary for possible later sensitivity analysis.

if ("TISSUE_SOURCE_SITE" %in% colnames(metadata)) {
  
  site_name_balance <- metadata %>%
    tibble::as_tibble() %>%
    dplyr::select(sample_id, TISSUE_SOURCE_SITE) %>%
    dplyr::left_join(
      labels %>%
        dplyr::select(sample_id, subtype),
      by = "sample_id"
    ) %>%
    dplyr::mutate(
      TISSUE_SOURCE_SITE = as.character(TISSUE_SOURCE_SITE),
      subtype = factor(as.character(subtype), levels = c("LUAD", "LUSC"))
    ) %>%
    dplyr::count(TISSUE_SOURCE_SITE, subtype, name = "n") %>%
    tidyr::pivot_wider(
      names_from = subtype,
      values_from = n,
      values_fill = 0
    )
  
  if (!"LUAD" %in% colnames(site_name_balance)) {
    site_name_balance$LUAD <- 0
  }
  
  if (!"LUSC" %in% colnames(site_name_balance)) {
    site_name_balance$LUSC <- 0
  }
  
  site_name_balance <- site_name_balance %>%
    dplyr::mutate(
      total = LUAD + LUSC,
      only_one_subtype = LUAD == 0 | LUSC == 0
    ) %>%
    dplyr::arrange(dplyr::desc(total))
  
  readr::write_tsv(
    site_name_balance,
    file.path(tables_dir, "eda_tissue_source_site_balance.tsv")
  )
  
  mixed_site_balance <- site_name_balance %>%
    dplyr::filter(LUAD > 0, LUSC > 0)
  
  mixed_site_summary <- tibble::tibble(
    variable = "TISSUE_SOURCE_SITE",
    mixed_site_groups = nrow(mixed_site_balance),
    mixed_site_samples = sum(mixed_site_balance$total),
    LUAD = sum(mixed_site_balance$LUAD),
    LUSC = sum(mixed_site_balance$LUSC)
  )
  
  readr::write_tsv(
    mixed_site_summary,
    file.path(tables_dir, "eda_mixed_source_site_subset_summary.tsv")
  )
  
  message("\nMixed-source-site subset summary:")
  print(mixed_site_summary)
  
  if (mixed_site_summary$mixed_site_samples > 0) {
    
    mixed_site_class_balance <- tibble::tibble(
      subtype = factor(c("LUAD", "LUSC"), levels = c("LUAD", "LUSC")),
      n = c(mixed_site_summary$LUAD, mixed_site_summary$LUSC)
    ) %>%
      dplyr::mutate(
        percent = round(n / sum(n) * 100, 2)
      )
    
    readr::write_tsv(
      mixed_site_class_balance,
      file.path(tables_dir, "eda_mixed_source_site_class_balance.tsv")
    )
    
    p_mixed_site <- ggplot2::ggplot(
      mixed_site_class_balance,
      ggplot2::aes(x = subtype, y = n, fill = subtype)
    ) +
      ggplot2::geom_col(width = 0.65, alpha = 0.92) +
      ggplot2::geom_text(
        ggplot2::aes(label = paste0(n, " samples\n", percent, "%")),
        vjust = -0.3,
        fontface = "bold",
        size = 4
      ) +
      ggplot2::scale_fill_manual(values = subtype_colors) +
      ggplot2::ylim(0, max(mixed_site_class_balance$n) * 1.2 + 1) +
      ggplot2::labs(
        title = "Mixed source-site sensitivity subset",
        subtitle = "Samples from broad source-site names containing both LUAD and LUSC",
        x = "Subtype",
        y = "Number of samples",
        caption = "This subset can be used later as a robustness check, not as the main dataset."
      ) +
      theme_eda() +
      ggplot2::theme(legend.position = "none")
    
    save_plot(
      p_mixed_site,
      "12_mixed_source_site_subset_class_balance.png",
      width = 7,
      height = 5
    )
    
  } else {
    
    message("\nNo mixed TISSUE_SOURCE_SITE groups found.")
    message("Skipping mixed-source-site subset plot.")
  }
}

# -----------------------------
# 14. Final EDA interpretation table
# -----------------------------

eda_conclusion <- tibble::tibble(
  section = c(
    "Class balance",
    "Feature counts",
    "Missingness",
    "PCA by subtype",
    "PCA by source/site",
    "Batch/source correction decision"
  ),
  conclusion = c(
    paste0(
      "The matched dataset contains ",
      class_balance$n[class_balance$subtype == "LUAD"],
      " LUAD and ",
      class_balance$n[class_balance$subtype == "LUSC"],
      " LUSC samples."
    ),
    paste0(
      "After basic QC, the dataset contains ",
      scales::comma(ncol(mrna)),
      " mRNA features, ",
      scales::comma(ncol(methylation)),
      " DNA methylation features and ",
      scales::comma(ncol(rppa)),
      " RPPA features."
    ),
    "After basic QC, mRNA and RPPA have no missing values, while DNA methylation has low remaining missingness. Remaining methylation missing values are not globally imputed and will be handled inside cross-validation folds.",
    "PCA is used to inspect whether major variation separates LUAD and LUSC within each omics layer.",
    "PCA colored by source/site variables is used to inspect whether sample clustering may also reflect tissue source or collection site.",
    "Source/site variables are inspected but not corrected because they are strongly confounded with subtype. Correcting them could remove biological LUAD/LUSC signal."
  )
)

readr::write_tsv(
  eda_conclusion,
  file.path(tables_dir, "eda_interpretation_summary.tsv")
)

message("\nEDA interpretation summary:")
print(eda_conclusion)

# -----------------------------
# 15. Save compact EDA object
# -----------------------------

eda_results <- list(
  class_balance = class_balance,
  dimension_summary = dimension_summary,
  missing_summary_after_qc = missing_summary_after_qc,
  pca_variance_summary = pca_variance_summary,
  pca_scores = pca_scores_all,
  confounding_overview = confounding_overview,
  confounding_tables = confounding_tables,
  eda_conclusion = eda_conclusion,
  note = paste(
    "EDA PCA used temporary imputation, scaling and top-variable selection only for visualization.",
    "No batch correction was applied.",
    "Source/site variables were inspected and found to be confounded with subtype, so they should not be used for direct correction."
  )
)

saveRDS(
  eda_results,
  file.path(processed_dir, "03_eda_results.rds"),
  compress = FALSE
)

message("\nSaved EDA results object to:")
message(file.path(processed_dir, "03_eda_results.rds"))

message("\nEDA figures saved in:")
message(figures_eda_dir)

message("\n03_eda.R completed successfully.")