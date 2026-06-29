############################################################
# 05_early_fusion_models.R
# Early-fusion multi-omics baseline modelling
# NSCLC multi-omics project: LUAD vs LUSC
#
# Input:
#   data/processed/02_qc_model_input_unimputed.rds
#   data/processed/04_single_omics_results.rds   optional, for comparison
#
# Output:
#   data/processed/05_early_fusion_results.rds
#   results/tables/early_fusion_*.tsv
#   results/figures/early_fusion/*.png
#
# Model:
#   - Early-fusion Elastic Net
#
# Integration strategy:
#   - For each CV fold, preprocess each omics layer using training data only
#   - Concatenate the processed mRNA, DNA methylation and RPPA matrices
#   - Fit Elastic Net on the combined feature matrix
#
# Leakage control:
#   Inside each outer CV fold, using training data only:
#     - median imputation
#     - variance-based top-feature filtering per omics layer
#     - scaling per omics layer
#     - model fitting and lambda tuning
#
# The test fold is transformed only using training-fold parameters.
############################################################

# -----------------------------
# 1. Setup
# -----------------------------

source("R/00_setup.R")

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "matrixStats",
  "glmnet",
  "pROC",
  "readr",
  "tibble",
  "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(matrixStats)
  library(glmnet)
  library(pROC)
  library(readr)
  library(tibble)
  library(scales)
})

qc_rds_path <- file.path(processed_dir, "02_qc_model_input_unimputed.rds")

if (!file.exists(qc_rds_path)) {
  stop(
    "02_qc_model_input_unimputed.rds not found.\n",
    "Run source('R/02_preprocessing.R') first."
  )
}

qc_data <- readRDS(qc_rds_path)

labels <- qc_data$labels
mrna <- qc_data$mrna
methylation <- qc_data$methylation
rppa <- qc_data$rppa

labels$sample_id <- as.character(labels$sample_id)
labels$subtype <- factor(labels$subtype, levels = c("LUAD", "LUSC"))

stopifnot(identical(labels$sample_id, rownames(mrna)))
stopifnot(identical(labels$sample_id, rownames(methylation)))
stopifnot(identical(labels$sample_id, rownames(rppa)))

figures_early_dir <- file.path(figures_dir, "early_fusion")
dir.create(figures_early_dir, recursive = TRUE, showWarnings = FALSE)

message("Loaded QC modelling object for early fusion.")
message("Samples: ", nrow(labels))
message("mRNA features: ", ncol(mrna))
message("DNA methylation features: ", ncol(methylation))
message("RPPA features: ", ncol(rppa))

# -----------------------------
# 2. Modelling plan
# -----------------------------

set.seed(123)

outer_k <- 5
inner_k_glmnet <- 5
glmnet_alpha <- 0.5
glmnet_lambda_choice <- "lambda.min"

# Fold-wise top-variable filtering per layer.
# These values keep the integrated matrix computationally manageable.
# They are not global biological feature selection.
early_fusion_feature_plan <- tibble::tibble(
  layer = c("mRNA", "DNA methylation", "RPPA"),
  top_n = c(3000, 5000, ncol(rppa))
)

readr::write_tsv(
  early_fusion_feature_plan,
  file.path(tables_dir, "early_fusion_feature_plan.tsv")
)

omics_list <- list(
  "mRNA" = mrna,
  "DNA methylation" = methylation,
  "RPPA" = rppa
)

# -----------------------------
# 3. Plot theme and colors
# -----------------------------

model_colors <- c(
  "Elastic Net" = "#1F77B4",
  "Random Forest" = "#D62728",
  "Early Fusion Elastic Net" = "#2CA02C"
)

layer_colors <- c(
  "mRNA" = "#1F77B4",
  "DNA methylation" = "#2CA02C",
  "RPPA" = "#9467BD"
)

subtype_colors <- c(
  "LUAD" = "#1F77B4",
  "LUSC" = "#D62728"
)

theme_model <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 3),
      plot.subtitle = ggplot2::element_text(size = base_size, color = "grey30"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "grey20"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(color = "grey35", size = base_size - 2)
    )
}

save_plot <- function(plot, filename, width = 8, height = 6, dpi = 320) {
  ggplot2::ggsave(
    filename = file.path(figures_early_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

# -----------------------------
# 4. Helper functions
# -----------------------------

safe_divide <- function(a, b) {
  ifelse(b == 0, NA_real_, a / b)
}

create_stratified_folds <- function(y, k = 5, seed = 123) {
  set.seed(seed)

  y <- factor(y, levels = c("LUAD", "LUSC"))
  fold_id <- integer(length(y))

  for (cls in levels(y)) {
    idx <- which(y == cls)
    idx <- sample(idx)
    fold_id[idx] <- rep(seq_len(k), length.out = length(idx))
  }

  split(seq_along(y), fold_id)
}

impute_with_medians <- function(x, medians) {
  missing_index <- which(is.na(x), arr.ind = TRUE)

  if (nrow(missing_index) > 0) {
    x[missing_index] <- medians[missing_index[, 2]]
  }

  x
}

make_unique_prefixed_names <- function(layer_name, features) {
  prefix <- dplyr::case_when(
    layer_name == "mRNA" ~ "mRNA",
    layer_name == "DNA methylation" ~ "METH",
    layer_name == "RPPA" ~ "RPPA",
    TRUE ~ safe_filename(layer_name)
  )

  make.unique(paste(prefix, features, sep = "__"))
}

fit_layer_preprocess_train <- function(x_train,
                                       layer_name,
                                       top_n) {
  feature_var <- matrixStats::colVars(x_train, na.rm = TRUE)
  names(feature_var) <- colnames(x_train)

  keep <- !is.na(feature_var) & feature_var > 0

  x_train <- x_train[, keep, drop = FALSE]
  feature_var <- feature_var[keep]

  n_select <- min(top_n, ncol(x_train))

  selected_features <- names(
    sort(feature_var, decreasing = TRUE)
  )[seq_len(n_select)]

  x_train <- x_train[, selected_features, drop = FALSE]

  feature_medians <- matrixStats::colMedians(x_train, na.rm = TRUE)
  names(feature_medians) <- colnames(x_train)
  feature_medians[!is.finite(feature_medians)] <- 0

  x_train <- impute_with_medians(x_train, feature_medians)

  feature_sd_after_impute <- matrixStats::colSds(x_train, na.rm = TRUE)
  valid_after_impute <- is.finite(feature_sd_after_impute) &
    feature_sd_after_impute > 0

  x_train <- x_train[, valid_after_impute, drop = FALSE]
  feature_medians <- feature_medians[valid_after_impute]
  selected_features <- selected_features[valid_after_impute]

  feature_center <- colMeans(x_train)
  feature_scale <- matrixStats::colSds(x_train)
  feature_scale[!is.finite(feature_scale) | feature_scale == 0] <- 1

  x_train <- sweep(x_train, 2, feature_center, "-")
  x_train <- sweep(x_train, 2, feature_scale, "/")

  prefixed_features <- make_unique_prefixed_names(layer_name, selected_features)

  colnames(x_train) <- prefixed_features

  feature_map <- tibble::tibble(
    fused_feature = prefixed_features,
    layer = layer_name,
    original_feature = selected_features
  )

  list(
    x_train = x_train,
    layer_name = layer_name,
    selected_features = selected_features,
    medians = feature_medians,
    center = feature_center,
    scale = feature_scale,
    feature_map = feature_map
  )
}

apply_layer_preprocess_test <- function(x_test, prep) {
  x_test <- x_test[, prep$selected_features, drop = FALSE]

  x_test <- impute_with_medians(x_test, prep$medians)

  x_test <- sweep(x_test, 2, prep$center, "-")
  x_test <- sweep(x_test, 2, prep$scale, "/")

  colnames(x_test) <- prep$feature_map$fused_feature

  x_test
}

fit_fusion_preprocess_train <- function(omics_list,
                                        train_idx,
                                        feature_plan) {
  prep_list <- list()
  train_layers <- list()
  feature_maps <- list()
  feature_counts <- list()

  for (layer_name in names(omics_list)) {
    top_n <- feature_plan %>%
      dplyr::filter(layer == layer_name) %>%
      dplyr::pull(top_n)

    x_train_layer <- omics_list[[layer_name]][train_idx, , drop = FALSE]

    prep <- fit_layer_preprocess_train(
      x_train = x_train_layer,
      layer_name = layer_name,
      top_n = top_n
    )

    prep_list[[layer_name]] <- prep
    train_layers[[layer_name]] <- prep$x_train
    feature_maps[[layer_name]] <- prep$feature_map

    feature_counts[[layer_name]] <- tibble::tibble(
      layer = layer_name,
      n_features_used = ncol(prep$x_train)
    )
  }

  x_train_fused <- do.call(cbind, train_layers)

  list(
    x_train = x_train_fused,
    prep_list = prep_list,
    feature_map = dplyr::bind_rows(feature_maps),
    feature_counts = dplyr::bind_rows(feature_counts)
  )
}

apply_fusion_preprocess_test <- function(omics_list,
                                         test_idx,
                                         fusion_prep) {
  test_layers <- list()

  for (layer_name in names(fusion_prep$prep_list)) {
    x_test_layer <- omics_list[[layer_name]][test_idx, , drop = FALSE]

    test_layers[[layer_name]] <- apply_layer_preprocess_test(
      x_test = x_test_layer,
      prep = fusion_prep$prep_list[[layer_name]]
    )
  }

  do.call(cbind, test_layers)
}

compute_binary_metrics <- function(actual, predicted, probability_lusc) {
  actual <- factor(actual, levels = c("LUAD", "LUSC"))
  predicted <- factor(predicted, levels = c("LUAD", "LUSC"))

  cm <- table(actual, predicted)

  tn <- as.numeric(cm["LUAD", "LUAD"])
  fp <- as.numeric(cm["LUAD", "LUSC"])
  fn <- as.numeric(cm["LUSC", "LUAD"])
  tp <- as.numeric(cm["LUSC", "LUSC"])

  precision_lusc <- safe_divide(tp, tp + fp)
  recall_lusc <- safe_divide(tp, tp + fn)

  f1_lusc <- safe_divide(
    2 * precision_lusc * recall_lusc,
    precision_lusc + recall_lusc
  )

  precision_luad <- safe_divide(tn, tn + fn)
  recall_luad <- safe_divide(tn, tn + fp)

  f1_luad <- safe_divide(
    2 * precision_luad * recall_luad,
    precision_luad + recall_luad
  )

  sensitivity_lusc <- recall_lusc
  specificity_luad <- recall_luad

  accuracy <- safe_divide(tp + tn, tp + tn + fp + fn)

  balanced_accuracy <- mean(
    c(sensitivity_lusc, specificity_luad),
    na.rm = TRUE
  )

  macro_f1 <- mean(
    c(f1_luad, f1_lusc),
    na.rm = TRUE
  )

  auc_value <- as.numeric(
    pROC::auc(
      response = actual,
      predictor = probability_lusc,
      levels = c("LUAD", "LUSC"),
      direction = "<",
      quiet = TRUE
    )
  )

  tibble::tibble(
    ROC_AUC = auc_value,
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    F1_LUAD = f1_luad,
    F1_LUSC = f1_lusc,
    macro_F1 = macro_f1,
    sensitivity_LUSC = sensitivity_lusc,
    specificity_LUAD = specificity_luad
  )
}

make_roc_data <- function(prediction_df) {
  roc_obj <- pROC::roc(
    response = factor(prediction_df$actual, levels = c("LUAD", "LUSC")),
    predictor = prediction_df$probability_lusc,
    levels = c("LUAD", "LUSC"),
    direction = "<",
    quiet = TRUE
  )

  tibble::tibble(
    model_label = "Early fusion + Elastic Net",
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities
  ) %>%
    dplyr::arrange(FPR, TPR)
}

# -----------------------------
# 5. Run early-fusion model
# -----------------------------

folds <- create_stratified_folds(
  y = labels$subtype,
  k = outer_k,
  seed = 123
)

all_predictions <- list()
all_fold_metrics <- list()
all_feature_maps <- list()
all_feature_counts <- list()
all_coefficients <- list()

for (fold_id in seq_along(folds)) {
  message("\n========================================")
  message("Early fusion | Fold: ", fold_id, "/", outer_k)
  message("========================================")

  test_idx <- folds[[fold_id]]
  train_idx <- setdiff(seq_len(nrow(labels)), test_idx)

  y_train <- labels$subtype[train_idx]
  y_test <- labels$subtype[test_idx]
  y_train_binary <- ifelse(y_train == "LUSC", 1, 0)

  fusion_prep <- fit_fusion_preprocess_train(
    omics_list = omics_list,
    train_idx = train_idx,
    feature_plan = early_fusion_feature_plan
  )

  x_train <- fusion_prep$x_train
  x_test <- apply_fusion_preprocess_test(
    omics_list = omics_list,
    test_idx = test_idx,
    fusion_prep = fusion_prep
  )

  message("Fused training matrix: ", nrow(x_train), " samples x ", ncol(x_train), " features")

  cv_fit <- glmnet::cv.glmnet(
    x = x_train,
    y = y_train_binary,
    family = "binomial",
    alpha = glmnet_alpha,
    nfolds = inner_k_glmnet,
    type.measure = "auc",
    standardize = FALSE
  )

  probability_lusc <- as.numeric(
    predict(
      cv_fit,
      newx = x_test,
      s = glmnet_lambda_choice,
      type = "response"
    )
  )

  predicted <- ifelse(probability_lusc >= 0.5, "LUSC", "LUAD")

  prediction_df <- tibble::tibble(
    model_label = "Early fusion + Elastic Net",
    model = "Elastic Net",
    integration = "Early fusion",
    fold = fold_id,
    sample_id = labels$sample_id[test_idx],
    actual = as.character(y_test),
    predicted = predicted,
    probability_lusc = probability_lusc
  )

  fold_metrics <- compute_binary_metrics(
    actual = prediction_df$actual,
    predicted = prediction_df$predicted,
    probability_lusc = prediction_df$probability_lusc
  ) %>%
    dplyr::mutate(
      model_label = "Early fusion + Elastic Net",
      model = "Elastic Net",
      integration = "Early fusion",
      fold = fold_id,
      n_features_used = ncol(x_train),
      lambda = cv_fit[[glmnet_lambda_choice]],
      .before = 1
    )

  coef_mat <- as.matrix(
    coef(cv_fit, s = glmnet_lambda_choice)
  )

  coef_df <- tibble::tibble(
    fused_feature = rownames(coef_mat),
    coefficient = as.numeric(coef_mat[, 1])
  ) %>%
    dplyr::filter(fused_feature != "(Intercept)", coefficient != 0) %>%
    dplyr::left_join(fusion_prep$feature_map, by = "fused_feature") %>%
    dplyr::mutate(
      abs_coefficient = abs(coefficient),
      fold = fold_id,
      model_label = "Early fusion + Elastic Net",
      model = "Elastic Net",
      integration = "Early fusion"
    ) %>%
    dplyr::arrange(dplyr::desc(abs_coefficient))

  feature_counts <- fusion_prep$feature_counts %>%
    dplyr::mutate(
      fold = fold_id,
      model_label = "Early fusion + Elastic Net"
    )

  all_predictions[[fold_id]] <- prediction_df
  all_fold_metrics[[fold_id]] <- fold_metrics
  all_feature_maps[[fold_id]] <- fusion_prep$feature_map %>%
    dplyr::mutate(fold = fold_id)
  all_feature_counts[[fold_id]] <- feature_counts
  all_coefficients[[fold_id]] <- coef_df
}

predictions_all <- dplyr::bind_rows(all_predictions) %>%
  dplyr::mutate(
    actual = factor(actual, levels = c("LUAD", "LUSC")),
    predicted = factor(predicted, levels = c("LUAD", "LUSC"))
  )

fold_metrics_all <- dplyr::bind_rows(all_fold_metrics)
feature_maps_all <- dplyr::bind_rows(all_feature_maps)
feature_counts_all <- dplyr::bind_rows(all_feature_counts)
coefficients_all <- dplyr::bind_rows(all_coefficients)

# -----------------------------
# 6. Metrics and summaries
# -----------------------------

overall_metrics <- compute_binary_metrics(
  actual = predictions_all$actual,
  predicted = predictions_all$predicted,
  probability_lusc = predictions_all$probability_lusc
) %>%
  dplyr::mutate(
    model_label = "Early fusion + Elastic Net",
    model = "Elastic Net",
    integration = "Early fusion",
    n_samples = nrow(predictions_all),
    .before = 1
  )

metric_summary <- fold_metrics_all %>%
  dplyr::summarise(
    model_label = "Early fusion + Elastic Net",
    mean_ROC_AUC = round(mean(ROC_AUC, na.rm = TRUE), 4),
    sd_ROC_AUC = round(sd(ROC_AUC, na.rm = TRUE), 4),
    mean_balanced_accuracy = round(mean(balanced_accuracy, na.rm = TRUE), 4),
    sd_balanced_accuracy = round(sd(balanced_accuracy, na.rm = TRUE), 4),
    mean_macro_F1 = round(mean(macro_F1, na.rm = TRUE), 4),
    sd_macro_F1 = round(sd(macro_F1, na.rm = TRUE), 4),
    mean_F1_LUAD = round(mean(F1_LUAD, na.rm = TRUE), 4),
    mean_F1_LUSC = round(mean(F1_LUSC, na.rm = TRUE), 4)
  )

classes <- factor(c("LUAD", "LUSC"), levels = c("LUAD", "LUSC"))

confusion_counts_raw <- predictions_all %>%
  dplyr::count(actual, predicted, name = "n")

confusion_counts <- tidyr::expand_grid(
  actual = classes,
  predicted = classes
) %>%
  dplyr::left_join(
    confusion_counts_raw,
    by = c("actual", "predicted")
  ) %>%
  dplyr::mutate(
    n = tidyr::replace_na(n, 0L),
    model_label = "Early fusion + Elastic Net"
  )

coefficient_summary <- coefficients_all %>%
  dplyr::group_by(layer, original_feature, fused_feature) %>%
  dplyr::summarise(
    selected_in_folds = dplyr::n_distinct(fold),
    mean_abs_coefficient = mean(abs_coefficient, na.rm = TRUE),
    mean_coefficient = mean(coefficient, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(selected_in_folds),
    dplyr::desc(mean_abs_coefficient)
  )

selected_feature_count_by_layer <- coefficient_summary %>%
  dplyr::group_by(layer) %>%
  dplyr::summarise(
    selected_features = dplyr::n(),
    mean_selected_in_folds = round(mean(selected_in_folds), 2),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(selected_features))

# -----------------------------
# 7. Save tables
# -----------------------------

readr::write_tsv(
  predictions_all,
  file.path(tables_dir, "early_fusion_out_of_fold_predictions.tsv")
)

readr::write_tsv(
  fold_metrics_all,
  file.path(tables_dir, "early_fusion_fold_metrics.tsv")
)

readr::write_tsv(
  overall_metrics,
  file.path(tables_dir, "early_fusion_overall_metrics.tsv")
)

readr::write_tsv(
  metric_summary,
  file.path(tables_dir, "early_fusion_metric_summary.tsv")
)

readr::write_tsv(
  confusion_counts,
  file.path(tables_dir, "early_fusion_confusion_matrix.tsv")
)

readr::write_tsv(
  feature_maps_all,
  file.path(tables_dir, "early_fusion_feature_map_by_fold.tsv")
)

readr::write_tsv(
  feature_counts_all,
  file.path(tables_dir, "early_fusion_feature_counts_by_fold.tsv")
)

readr::write_tsv(
  coefficients_all,
  file.path(tables_dir, "early_fusion_coefficients_by_fold.tsv")
)

readr::write_tsv(
  coefficient_summary,
  file.path(tables_dir, "early_fusion_coefficient_summary.tsv")
)

readr::write_tsv(
  selected_feature_count_by_layer,
  file.path(tables_dir, "early_fusion_selected_feature_count_by_layer.tsv")
)

message("\nEarly-fusion overall metrics:")
print(overall_metrics)

# -----------------------------
# 8. Visualizations
# -----------------------------

roc_data <- make_roc_data(predictions_all)

readr::write_tsv(
  roc_data,
  file.path(tables_dir, "early_fusion_roc_curve_data.tsv")
)

p_roc <- ggplot2::ggplot(
  roc_data,
  ggplot2::aes(x = FPR, y = TPR)
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey55"
  ) +
  ggplot2::geom_line(color = "#2CA02C", linewidth = 1.2, alpha = 0.95) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Early-fusion ROC curve",
    subtitle = "Out-of-fold predictions from 5-fold stratified cross-validation",
    x = "False positive rate",
    y = "True positive rate",
    caption = "mRNA, DNA methylation and RPPA were concatenated after fold-wise preprocessing."
  ) +
  theme_model()

save_plot(
  p_roc,
  "01_early_fusion_roc_curve.png",
  width = 7,
  height = 6
)

metric_plot_df <- overall_metrics %>%
  dplyr::select(
    ROC_AUC,
    balanced_accuracy,
    macro_F1,
    F1_LUAD,
    F1_LUSC
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(
      metric,
      ROC_AUC = "ROC-AUC",
      balanced_accuracy = "Balanced accuracy",
      macro_F1 = "Macro-F1",
      F1_LUAD = "F1 LUAD",
      F1_LUSC = "F1 LUSC"
    )
  )

p_metrics <- ggplot2::ggplot(
  metric_plot_df,
  ggplot2::aes(x = reorder(metric, value), y = value)
) +
  ggplot2::geom_col(fill = "#2CA02C", width = 0.7, alpha = 0.92) +
  ggplot2::coord_flip() +
  ggplot2::geom_text(
    ggplot2::aes(label = round(value, 3)),
    hjust = -0.1,
    fontface = "bold",
    size = 4
  ) +
  ggplot2::scale_y_continuous(limits = c(0, 1.08)) +
  ggplot2::labs(
    title = "Early-fusion Elastic Net performance",
    subtitle = "Performance based on out-of-fold predictions",
    x = "Metric",
    y = "Metric value"
  ) +
  theme_model()

save_plot(
  p_metrics,
  "02_early_fusion_metric_summary.png",
  width = 8,
  height = 5
)

p_confusion <- ggplot2::ggplot(
  confusion_counts,
  ggplot2::aes(x = predicted, y = actual, fill = n)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.8) +
  ggplot2::geom_text(
    ggplot2::aes(label = n),
    size = 6,
    fontface = "bold",
    color = "white"
  ) +
  ggplot2::scale_fill_viridis_c(option = "magma", direction = -1) +
  ggplot2::labs(
    title = "Early-fusion confusion matrix",
    subtitle = "Counts are based on out-of-fold predictions",
    x = "Predicted subtype",
    y = "Actual subtype",
    fill = "Count"
  ) +
  theme_model()

save_plot(
  p_confusion,
  "03_early_fusion_confusion_matrix.png",
  width = 6.5,
  height = 5.5
)

p_layer_features <- ggplot2::ggplot(
  selected_feature_count_by_layer,
  ggplot2::aes(x = reorder(layer, selected_features), y = selected_features, fill = layer)
) +
  ggplot2::geom_col(width = 0.7, alpha = 0.92) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = layer_colors) +
  ggplot2::geom_text(
    ggplot2::aes(label = selected_features),
    hjust = -0.1,
    fontface = "bold",
    size = 4
  ) +
  ggplot2::labs(
    title = "Selected early-fusion features by omics layer",
    subtitle = "Non-zero Elastic Net coefficients aggregated across folds",
    x = "Omics layer",
    y = "Number of selected features",
    fill = "Layer"
  ) +
  theme_model() +
  ggplot2::theme(legend.position = "none")

save_plot(
  p_layer_features,
  "04_early_fusion_selected_features_by_layer.png",
  width = 8,
  height = 5
)

top_features <- coefficient_summary %>%
  dplyr::slice_max(
    order_by = mean_abs_coefficient,
    n = 30,
    with_ties = FALSE
  ) %>%
  dplyr::mutate(
    feature_label = paste0(layer, ": ", original_feature),
    feature_label = reorder(feature_label, mean_abs_coefficient)
  )

p_top_features <- ggplot2::ggplot(
  top_features,
  ggplot2::aes(x = feature_label, y = mean_abs_coefficient, fill = layer)
) +
  ggplot2::geom_col(width = 0.75, alpha = 0.92) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = layer_colors) +
  ggplot2::labs(
    title = "Top early-fusion Elastic Net features",
    subtitle = "Ranked by mean absolute coefficient across CV folds",
    x = "Feature",
    y = "Mean absolute coefficient",
    fill = "Layer"
  ) +
  theme_model(base_size = 11)

save_plot(
  p_top_features,
  "05_early_fusion_top_features.png",
  width = 10,
  height = 8
)

# -----------------------------
# 9. Optional comparison with single-omics baselines
# -----------------------------

single_results_path <- file.path(processed_dir, "04_single_omics_results.rds")

comparison_metrics <- NULL

if (file.exists(single_results_path)) {
  single_results <- readRDS(single_results_path)

  single_overall <- single_results$overall_metrics %>%
    dplyr::mutate(
      analysis = "Single-omics",
      model_label = paste(layer, model, sep = " + ")
    ) %>%
    dplyr::select(
      analysis,
      model_label,
      ROC_AUC,
      balanced_accuracy,
      macro_F1,
      F1_LUAD,
      F1_LUSC
    )

  early_overall <- overall_metrics %>%
    dplyr::mutate(
      analysis = "Multi-omics early fusion"
    ) %>%
    dplyr::select(
      analysis,
      model_label,
      ROC_AUC,
      balanced_accuracy,
      macro_F1,
      F1_LUAD,
      F1_LUSC
    )

  comparison_metrics <- dplyr::bind_rows(
    single_overall,
    early_overall
  ) %>%
    dplyr::arrange(dplyr::desc(ROC_AUC))

  readr::write_tsv(
    comparison_metrics,
    file.path(tables_dir, "early_fusion_vs_single_omics_metrics.tsv")
  )

  p_compare_auc <- ggplot2::ggplot(
    comparison_metrics,
    ggplot2::aes(x = reorder(model_label, ROC_AUC), y = ROC_AUC, fill = analysis)
  ) +
    ggplot2::geom_col(width = 0.7, alpha = 0.92) +
    ggplot2::coord_flip() +
    ggplot2::geom_text(
      ggplot2::aes(label = round(ROC_AUC, 3)),
      hjust = -0.1,
      fontface = "bold",
      size = 3.6
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1.08)) +
    ggplot2::scale_fill_manual(
      values = c(
        "Single-omics" = "#7F7F7F",
        "Multi-omics early fusion" = "#2CA02C"
      )
    ) +
    ggplot2::labs(
      title = "Early fusion compared with single-omics baselines",
      subtitle = "ROC-AUC from out-of-fold predictions",
      x = "Model",
      y = "ROC-AUC",
      fill = "Analysis"
    ) +
    theme_model()

  save_plot(
    p_compare_auc,
    "06_early_fusion_vs_single_omics_auc.png",
    width = 10,
    height = 7
  )
}

# -----------------------------
# 10. Interpretation and save object
# -----------------------------

best_text <- paste0(
  "Early-fusion Elastic Net achieved ROC-AUC = ",
  round(overall_metrics$ROC_AUC, 3),
  ", balanced accuracy = ",
  round(overall_metrics$balanced_accuracy, 3),
  ", and macro-F1 = ",
  round(overall_metrics$macro_F1, 3),
  "."
)

comparison_text <- if (!is.null(comparison_metrics)) {
  best_available <- comparison_metrics %>%
    dplyr::arrange(dplyr::desc(ROC_AUC)) %>%
    dplyr::slice(1)

  paste0(
    "When compared with available single-omics baselines, the top ROC-AUC model was ",
    best_available$model_label,
    " with ROC-AUC = ",
    round(best_available$ROC_AUC, 3),
    "."
  )
} else {
  "Single-omics results were not found, so cross-model comparison was skipped."
}

early_fusion_interpretation <- tibble::tibble(
  section = c(
    "Early-fusion strategy",
    "Leakage control",
    "Performance summary",
    "Comparison purpose",
    "Next analysis step"
  ),
  conclusion = c(
    "mRNA, DNA methylation and RPPA were preprocessed separately inside each CV fold and then concatenated into one integrated feature matrix.",
    "Median imputation, scaling and top-variable filtering were fitted only on the training fold, then applied to the held-out fold.",
    best_text,
    comparison_text,
    "The next step is DIABLO, which performs supervised multi-block integration and feature selection instead of simple concatenation."
  )
)

readr::write_tsv(
  early_fusion_interpretation,
  file.path(tables_dir, "early_fusion_interpretation_summary.tsv")
)

message("\nEarly-fusion interpretation summary:")
print(early_fusion_interpretation)

early_fusion_results <- list(
  feature_plan = early_fusion_feature_plan,
  predictions = predictions_all,
  fold_metrics = fold_metrics_all,
  overall_metrics = overall_metrics,
  metric_summary = metric_summary,
  confusion_counts = confusion_counts,
  roc_data = roc_data,
  feature_maps_by_fold = feature_maps_all,
  feature_counts_by_fold = feature_counts_all,
  coefficients_by_fold = coefficients_all,
  coefficient_summary = coefficient_summary,
  selected_feature_count_by_layer = selected_feature_count_by_layer,
  comparison_metrics = comparison_metrics,
  interpretation = early_fusion_interpretation,
  note = paste(
    "Early fusion concatenated mRNA, DNA methylation and RPPA after fold-wise preprocessing.",
    "No clinical/source-site metadata variables were used as predictors.",
    "All imputation, scaling and feature filtering were fitted inside cross-validation training folds only.",
    "This is a multi-omics integration baseline before DIABLO."
  )
)

saveRDS(
  early_fusion_results,
  file.path(processed_dir, "05_early_fusion_results.rds"),
  compress = FALSE
)

message("\nSaved early-fusion results to:")
message(file.path(processed_dir, "05_early_fusion_results.rds"))

message("\nEarly-fusion figures saved in:")
message(figures_early_dir)

message("\n05_early_fusion_models.R completed successfully.")
