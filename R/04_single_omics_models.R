############################################################
# 04_single_omics_models.R
# Single-omics baseline modelling
# NSCLC multi-omics project: LUAD vs LUSC
#
# Input:
#   data/processed/02_qc_model_input_unimputed.rds
#
# Output:
#   data/processed/04_single_omics_results.rds
#   results/tables/single_omics_*.tsv
#   results/figures/single_omics/*.png
#
# Models:
#   - Elastic Net
#   - Random Forest
#
# Omics layers:
#   - mRNA
#   - DNA methylation
#   - RPPA
#
# Important leakage control:
#   Inside each outer CV fold, using training data only:
#   - median imputation
#   - variance-based top-feature filtering
#   - scaling for Elastic Net
#   - model fitting
#
#   The test fold is transformed only using training-fold parameters.
############################################################


# -----------------------------
# 1. Setup
# -----------------------------

source("R/00_setup.R")

required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "matrixStats",
  "glmnet", "ranger", "pROC", "readr", "tibble", "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(matrixStats)
  library(glmnet)
  library(ranger)
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

qc_data     <- readRDS(qc_rds_path)
labels      <- qc_data$labels
mrna        <- qc_data$mrna
methylation <- qc_data$methylation
rppa        <- qc_data$rppa

labels$sample_id <- as.character(labels$sample_id)
labels$subtype   <- factor(labels$subtype, levels = c("LUAD", "LUSC"))

stopifnot(identical(labels$sample_id, rownames(mrna)))
stopifnot(identical(labels$sample_id, rownames(methylation)))
stopifnot(identical(labels$sample_id, rownames(rppa)))

figures_single_dir <- file.path(figures_dir, "single_omics")
dir.create(figures_single_dir, recursive = TRUE, showWarnings = FALSE)

message("Loaded QC modelling object.")
message("Samples: ",              nrow(labels))
message("mRNA features: ",        ncol(mrna))
message("DNA methylation features: ", ncol(methylation))
message("RPPA features: ",        ncol(rppa))


# -----------------------------
# 2. Modelling plan
# -----------------------------

set.seed(123)

outer_k             <- 5
inner_k_glmnet      <- 5
glmnet_alpha        <- 0.5
glmnet_lambda_choice <- "lambda.min"
rf_num_trees        <- 500

# Top-variable filtering is done INSIDE each CV training fold.
# These are computational limits, not biological conclusions.
model_plan <- tibble::tibble(
  layer                = c("mRNA", "DNA methylation", "RPPA"),
  elastic_net_top_n    = c(5000, 10000, ncol(rppa)),
  random_forest_top_n  = c(2000, 3000,  ncol(rppa))
)

readr::write_tsv(
  model_plan,
  file.path(tables_dir, "single_omics_model_plan.tsv")
)

omics_list <- list(
  "mRNA"            = mrna,
  "DNA methylation" = methylation,
  "RPPA"            = rppa
)


# -----------------------------
# 3. Plot theme and colours
# -----------------------------

model_colors <- c(
  "Elastic Net"   = "#1F77B4",
  "Random Forest" = "#D62728"
)

subtype_colors <- c(
  "LUAD" = "#1F77B4",
  "LUSC" = "#D62728"
)

theme_model <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = base_size + 3),
      plot.subtitle = ggplot2::element_text(size = base_size, color = "grey30"),
      axis.title    = ggplot2::element_text(face = "bold"),
      axis.text     = ggplot2::element_text(color = "grey20"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title  = ggplot2::element_text(face = "bold"),
      plot.caption  = ggplot2::element_text(color = "grey35", size = base_size - 2)
    )
}

save_plot <- function(plot, filename, width = 8, height = 6, dpi = 320) {
  ggplot2::ggsave(
    filename = file.path(figures_single_dir, filename),
    plot     = plot,
    width    = width,
    height   = height,
    dpi      = dpi,
    bg       = "white"
  )
}

safe_filename <- function(x) gsub("[^A-Za-z0-9]+", "_", x)


# -----------------------------
# 4. Helper functions
# -----------------------------

safe_divide <- function(a, b) ifelse(b == 0, NA_real_, a / b)

create_stratified_folds <- function(y, k = 5, seed = 123) {
  set.seed(seed)
  y       <- factor(y, levels = c("LUAD", "LUSC"))
  fold_id <- integer(length(y))
  for (cls in levels(y)) {
    idx          <- which(y == cls)
    idx          <- sample(idx)
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

fit_preprocess_train <- function(x_train, top_n, scale_features = TRUE) {
  feature_var        <- matrixStats::colVars(x_train, na.rm = TRUE)
  names(feature_var) <- colnames(x_train)
  
  keep        <- !is.na(feature_var) & feature_var > 0
  x_train     <- x_train[, keep, drop = FALSE]
  feature_var <- feature_var[keep]
  
  n_select          <- min(top_n, ncol(x_train))
  selected_features <- names(sort(feature_var, decreasing = TRUE))[seq_len(n_select)]
  x_train           <- x_train[, selected_features, drop = FALSE]
  
  feature_medians                             <- matrixStats::colMedians(x_train, na.rm = TRUE)
  names(feature_medians)                      <- colnames(x_train)
  feature_medians[!is.finite(feature_medians)] <- 0
  x_train <- impute_with_medians(x_train, feature_medians)
  
  feature_sd_after_impute <- matrixStats::colSds(x_train, na.rm = TRUE)
  valid_after_impute      <- is.finite(feature_sd_after_impute) & feature_sd_after_impute > 0
  x_train           <- x_train[, valid_after_impute, drop = FALSE]
  feature_medians   <- feature_medians[valid_after_impute]
  selected_features <- selected_features[valid_after_impute]
  
  if (scale_features) {
    feature_center <- colMeans(x_train)
    feature_scale  <- matrixStats::colSds(x_train)
    feature_scale[!is.finite(feature_scale) | feature_scale == 0] <- 1
    x_train <- sweep(x_train, 2, feature_center, "-")
    x_train <- sweep(x_train, 2, feature_scale,  "/")
  } else {
    feature_center <- rep(0, ncol(x_train))
    feature_scale  <- rep(1, ncol(x_train))
    names(feature_center) <- colnames(x_train)
    names(feature_scale)  <- colnames(x_train)
  }
  
  list(
    x_train           = x_train,
    selected_features = selected_features,
    medians           = feature_medians,
    center            = feature_center,
    scale             = feature_scale,
    scale_features    = scale_features
  )
}

apply_preprocess_test <- function(x_test, prep) {
  x_test <- x_test[, prep$selected_features, drop = FALSE]
  x_test <- impute_with_medians(x_test, prep$medians)
  if (isTRUE(prep$scale_features)) {
    x_test <- sweep(x_test, 2, prep$center, "-")
    x_test <- sweep(x_test, 2, prep$scale,  "/")
  }
  x_test
}

compute_binary_metrics <- function(actual, predicted, probability_lusc) {
  actual    <- factor(actual,    levels = c("LUAD", "LUSC"))
  predicted <- factor(predicted, levels = c("LUAD", "LUSC"))
  
  cm <- table(actual, predicted)
  tn <- as.numeric(cm["LUAD", "LUAD"])
  fp <- as.numeric(cm["LUAD", "LUSC"])
  fn <- as.numeric(cm["LUSC", "LUAD"])
  tp <- as.numeric(cm["LUSC", "LUSC"])
  
  precision_lusc <- safe_divide(tp, tp + fp)
  recall_lusc    <- safe_divide(tp, tp + fn)
  f1_lusc        <- safe_divide(
    2 * precision_lusc * recall_lusc,
    precision_lusc + recall_lusc
  )
  
  precision_luad <- safe_divide(tn, tn + fn)
  recall_luad    <- safe_divide(tn, tn + fp)
  f1_luad        <- safe_divide(
    2 * precision_luad * recall_luad,
    precision_luad + recall_luad
  )
  
  sensitivity_lusc <- recall_lusc
  specificity_luad <- recall_luad
  accuracy         <- safe_divide(tp + tn, tp + tn + fp + fn)
  balanced_accuracy <- mean(c(sensitivity_lusc, specificity_luad), na.rm = TRUE)
  macro_f1         <- mean(c(f1_luad, f1_lusc), na.rm = TRUE)
  
  auc_value <- as.numeric(
    pROC::auc(
      response  = actual,
      predictor = probability_lusc,
      levels    = c("LUAD", "LUSC"),
      direction = "<",
      quiet     = TRUE
    )
  )
  
  tibble::tibble(
    ROC_AUC           = auc_value,
    accuracy          = accuracy,
    balanced_accuracy = balanced_accuracy,
    F1_LUAD           = f1_luad,
    F1_LUSC           = f1_lusc,
    macro_F1          = macro_f1,
    sensitivity_LUSC  = sensitivity_lusc,
    specificity_LUAD  = specificity_luad
  )
}

make_roc_data <- function(prediction_df) {
  prediction_df %>%
    dplyr::group_by(layer, model) %>%
    dplyr::group_modify(function(.x, .y) {
      roc_obj <- pROC::roc(
        response  = factor(.x$actual, levels = c("LUAD", "LUSC")),
        predictor = .x$probability_lusc,
        levels    = c("LUAD", "LUSC"),
        direction = "<",
        quiet     = TRUE
      )
      tibble::tibble(
        FPR = 1 - roc_obj$specificities,
        TPR = roc_obj$sensitivities
      ) %>%
        dplyr::arrange(FPR, TPR)
    }) %>%
    dplyr::ungroup()
}


# -----------------------------
# 5. Model fitting functions
# -----------------------------

fit_elastic_net_fold <- function(
    x, y, fold_id, train_idx, test_idx,
    layer_name, top_n, seed = 123
) {
  set.seed(seed + fold_id)
  
  x_train_raw <- x[train_idx, , drop = FALSE]
  x_test_raw  <- x[test_idx,  , drop = FALSE]
  y_train     <- y[train_idx]
  y_test      <- y[test_idx]
  
  prep    <- fit_preprocess_train(x_train = x_train_raw, top_n = top_n, scale_features = TRUE)
  x_train <- prep$x_train
  x_test  <- apply_preprocess_test(x_test_raw, prep)
  
  y_train_binary <- ifelse(y_train == "LUSC", 1, 0)
  
  cv_fit <- glmnet::cv.glmnet(
    x            = x_train,
    y            = y_train_binary,
    family       = "binomial",
    alpha        = glmnet_alpha,
    nfolds       = inner_k_glmnet,
    type.measure = "auc",
    standardize  = FALSE
  )
  
  probability_lusc <- as.numeric(
    predict(cv_fit, newx = x_test, s = glmnet_lambda_choice, type = "response")
  )
  predicted <- ifelse(probability_lusc >= 0.5, "LUSC", "LUAD")
  
  prediction_df <- tibble::tibble(
    layer            = layer_name,
    model            = "Elastic Net",
    fold             = fold_id,
    sample_id        = labels$sample_id[test_idx],
    actual           = as.character(y_test),
    predicted        = predicted,
    probability_lusc = probability_lusc
  )
  
  fold_metrics <- compute_binary_metrics(
    actual           = prediction_df$actual,
    predicted        = prediction_df$predicted,
    probability_lusc = prediction_df$probability_lusc
  ) %>%
    dplyr::mutate(
      layer           = layer_name,
      model           = "Elastic Net",
      fold            = fold_id,
      n_features_used = ncol(x_train),
      lambda          = cv_fit[[glmnet_lambda_choice]],
      .before = 1
    )
  
  coef_mat <- as.matrix(coef(cv_fit, s = glmnet_lambda_choice))
  coef_df  <- tibble::tibble(
    feature     = rownames(coef_mat),
    coefficient = as.numeric(coef_mat[, 1])
  ) %>%
    dplyr::filter(feature != "(Intercept)", coefficient != 0) %>%
    dplyr::mutate(
      abs_coefficient = abs(coefficient),
      layer = layer_name,
      model = "Elastic Net",
      fold  = fold_id
    ) %>%
    dplyr::arrange(dplyr::desc(abs_coefficient))
  
  list(
    predictions       = prediction_df,
    fold_metrics      = fold_metrics,
    feature_importance = coef_df
  )
}

fit_random_forest_fold <- function(
    x, y, fold_id, train_idx, test_idx,
    layer_name, top_n, seed = 123
) {
  set.seed(seed + 1000 + fold_id)
  
  x_train_raw <- x[train_idx, , drop = FALSE]
  x_test_raw  <- x[test_idx,  , drop = FALSE]
  y_train     <- y[train_idx]
  y_test      <- y[test_idx]
  
  prep    <- fit_preprocess_train(x_train = x_train_raw, top_n = top_n, scale_features = FALSE)
  x_train <- prep$x_train
  x_test  <- apply_preprocess_test(x_test_raw, prep)
  
  rf_feature_names <- make.names(colnames(x_train), unique = TRUE)
  feature_map      <- tibble::tibble(
    rf_feature = rf_feature_names,
    feature    = colnames(x_train)
  )
  colnames(x_train) <- rf_feature_names
  colnames(x_test)  <- rf_feature_names
  
  rf_train_df         <- as.data.frame(x_train, check.names = FALSE)
  rf_train_df$subtype <- factor(y_train, levels = c("LUAD", "LUSC"))
  rf_test_df          <- as.data.frame(x_test, check.names = FALSE)
  
  mtry_value <- max(1, floor(sqrt(ncol(x_train))))
  
  rf_fit <- ranger::ranger(
    formula   = subtype ~ .,
    data      = rf_train_df,
    probability = TRUE,
    num.trees = rf_num_trees,
    mtry      = mtry_value,
    importance = "impurity",
    seed      = seed + 1000 + fold_id
  )
  
  rf_pred          <- predict(rf_fit, data = rf_test_df)$predictions
  probability_lusc <- as.numeric(rf_pred[, "LUSC"])
  predicted        <- ifelse(probability_lusc >= 0.5, "LUSC", "LUAD")
  
  prediction_df <- tibble::tibble(
    layer            = layer_name,
    model            = "Random Forest",
    fold             = fold_id,
    sample_id        = labels$sample_id[test_idx],
    actual           = as.character(y_test),
    predicted        = predicted,
    probability_lusc = probability_lusc
  )
  
  fold_metrics <- compute_binary_metrics(
    actual           = prediction_df$actual,
    predicted        = prediction_df$predicted,
    probability_lusc = prediction_df$probability_lusc
  ) %>%
    dplyr::mutate(
      layer           = layer_name,
      model           = "Random Forest",
      fold            = fold_id,
      n_features_used = ncol(x_train),
      mtry            = mtry_value,
      num_trees       = rf_num_trees,
      .before = 1
    )
  
  importance_vec <- ranger::importance(rf_fit)
  importance_df  <- tibble::tibble(
    rf_feature = names(importance_vec),
    importance = as.numeric(importance_vec)
  ) %>%
    dplyr::left_join(feature_map, by = "rf_feature") %>%
    dplyr::select(feature, importance) %>%
    dplyr::arrange(dplyr::desc(importance)) %>%
    dplyr::slice_head(n = 100) %>%
    dplyr::mutate(layer = layer_name, model = "Random Forest", fold = fold_id)
  
  list(
    predictions        = prediction_df,
    fold_metrics       = fold_metrics,
    feature_importance = importance_df
  )
}


# -----------------------------
# 6. Run single-omics models
# -----------------------------

folds <- create_stratified_folds(y = labels$subtype, k = outer_k, seed = 123)

all_predictions    <- list()
all_fold_metrics   <- list()
all_elastic_features <- list()
all_rf_features    <- list()
counter            <- 1

for (layer_name in names(omics_list)) {
  message("\n========================================")
  message("Running single-omics models for: ", layer_name)
  message("========================================")
  
  x          <- omics_list[[layer_name]]
  y          <- labels$subtype
  layer_plan <- model_plan %>% dplyr::filter(layer == layer_name)
  elastic_top_n <- layer_plan$elastic_net_top_n[1]
  rf_top_n      <- layer_plan$random_forest_top_n[1]
  
  for (fold_id in seq_along(folds)) {
    message("\nLayer: ", layer_name, " | Fold: ", fold_id, "/", outer_k)
    
    test_idx  <- folds[[fold_id]]
    train_idx <- setdiff(seq_len(nrow(x)), test_idx)
    
    message("Elastic Net...")
    elastic_result <- fit_elastic_net_fold(
      x = x, y = y, fold_id = fold_id,
      train_idx = train_idx, test_idx = test_idx,
      layer_name = layer_name, top_n = elastic_top_n, seed = 123
    )
    
    message("Random Forest...")
    rf_result <- fit_random_forest_fold(
      x = x, y = y, fold_id = fold_id,
      train_idx = train_idx, test_idx = test_idx,
      layer_name = layer_name, top_n = rf_top_n, seed = 123
    )
    
    all_predictions[[counter]]      <- elastic_result$predictions
    all_fold_metrics[[counter]]     <- elastic_result$fold_metrics
    all_elastic_features[[counter]] <- elastic_result$feature_importance
    counter <- counter + 1
    
    all_predictions[[counter]]  <- rf_result$predictions
    all_fold_metrics[[counter]] <- rf_result$fold_metrics
    all_rf_features[[counter]]  <- rf_result$feature_importance
    counter <- counter + 1
  }
}

predictions_all <- dplyr::bind_rows(all_predictions) %>%
  dplyr::mutate(
    actual    = factor(actual,    levels = c("LUAD", "LUSC")),
    predicted = factor(predicted, levels = c("LUAD", "LUSC")),
    layer     = factor(layer,     levels = c("mRNA", "DNA methylation", "RPPA")),
    model     = factor(model,     levels = c("Elastic Net", "Random Forest"))
  )

fold_metrics_all <- dplyr::bind_rows(all_fold_metrics) %>%
  dplyr::mutate(
    layer = factor(layer, levels = c("mRNA", "DNA methylation", "RPPA")),
    model = factor(model, levels = c("Elastic Net", "Random Forest"))
  )

elastic_features_all <- dplyr::bind_rows(all_elastic_features)
rf_features_all      <- dplyr::bind_rows(all_rf_features)


# -----------------------------
# 7. Overall metrics and confusion matrices
# -----------------------------

overall_metrics <- predictions_all %>%
  dplyr::group_by(layer, model) %>%
  dplyr::group_modify(function(.x, .y) {
    compute_binary_metrics(
      actual           = .x$actual,
      predicted        = .x$predicted,
      probability_lusc = .x$probability_lusc
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(ROC_AUC))

metric_summary <- fold_metrics_all %>%
  dplyr::group_by(layer, model) %>%
  dplyr::summarise(
    mean_ROC_AUC          = round(mean(ROC_AUC,          na.rm = TRUE), 4),
    sd_ROC_AUC            = round(sd(ROC_AUC,            na.rm = TRUE), 4),
    mean_balanced_accuracy = round(mean(balanced_accuracy, na.rm = TRUE), 4),
    sd_balanced_accuracy  = round(sd(balanced_accuracy,  na.rm = TRUE), 4),
    mean_macro_F1         = round(mean(macro_F1,         na.rm = TRUE), 4),
    sd_macro_F1           = round(sd(macro_F1,           na.rm = TRUE), 4),
    mean_F1_LUAD          = round(mean(F1_LUAD,          na.rm = TRUE), 4),
    mean_F1_LUSC          = round(mean(F1_LUSC,          na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(mean_ROC_AUC))

classes <- factor(c("LUAD", "LUSC"), levels = c("LUAD", "LUSC"))

confusion_counts_raw <- predictions_all %>%
  dplyr::count(layer, model, actual, predicted, name = "n")

confusion_template <- tidyr::expand_grid(
  layer     = levels(predictions_all$layer),
  model     = levels(predictions_all$model),
  actual    = classes,
  predicted = classes
)

confusion_counts <- confusion_template %>%
  dplyr::left_join(
    confusion_counts_raw,
    by = c("layer", "model", "actual", "predicted")
  ) %>%
  dplyr::mutate(n = tidyr::replace_na(n, 0L))

readr::write_tsv(predictions_all,   file.path(tables_dir, "single_omics_out_of_fold_predictions.tsv"))
readr::write_tsv(fold_metrics_all,  file.path(tables_dir, "single_omics_fold_metrics.tsv"))
readr::write_tsv(overall_metrics,   file.path(tables_dir, "single_omics_overall_metrics.tsv"))
readr::write_tsv(metric_summary,    file.path(tables_dir, "single_omics_metric_summary.tsv"))
readr::write_tsv(confusion_counts,  file.path(tables_dir, "single_omics_confusion_matrices.tsv"))

message("\nOverall single-omics metrics:")
print(overall_metrics)


# -----------------------------
# 8. Feature importance summaries
# -----------------------------

elastic_feature_summary <- elastic_features_all %>%
  dplyr::group_by(layer, model, feature) %>%
  dplyr::summarise(
    selected_in_folds    = dplyr::n_distinct(fold),
    mean_abs_coefficient = mean(abs_coefficient, na.rm = TRUE),
    mean_coefficient     = mean(coefficient,     na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(layer, dplyr::desc(selected_in_folds), dplyr::desc(mean_abs_coefficient))

rf_feature_summary <- rf_features_all %>%
  dplyr::group_by(layer, model, feature) %>%
  dplyr::summarise(
    selected_in_folds = dplyr::n_distinct(fold),
    mean_importance   = mean(importance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(layer, dplyr::desc(mean_importance))

readr::write_tsv(elastic_features_all,    file.path(tables_dir, "single_omics_elastic_net_selected_features_by_fold.tsv"))
readr::write_tsv(elastic_feature_summary, file.path(tables_dir, "single_omics_elastic_net_feature_summary.tsv"))
readr::write_tsv(rf_features_all,         file.path(tables_dir, "single_omics_random_forest_top100_importance_by_fold.tsv"))
readr::write_tsv(rf_feature_summary,      file.path(tables_dir, "single_omics_random_forest_feature_summary.tsv"))


# -----------------------------
# 9. ROC curve visualisation
# -----------------------------

roc_data <- make_roc_data(predictions_all)
readr::write_tsv(roc_data, file.path(tables_dir, "single_omics_roc_curve_data.tsv"))

p_roc <- ggplot2::ggplot(
  roc_data,
  ggplot2::aes(x = FPR, y = TPR, color = model)
) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
  ggplot2::geom_line(linewidth = 1.15, alpha = 0.95) +
  ggplot2::facet_wrap(~layer, nrow = 1) +
  ggplot2::scale_color_manual(values = model_colors) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title    = "Single-omics ROC curves",
    subtitle = "Out-of-fold predictions from 5-fold stratified cross-validation",
    x        = "False positive rate",
    y        = "True positive rate",
    color    = "Model",
    caption  = "LUSC is treated as the positive class only for ROC probability direction."
  ) +
  theme_model()

save_plot(p_roc, "01_single_omics_roc_curves.png", width = 12, height = 5)


# -----------------------------
# 10. Metric comparison plots
# -----------------------------

metric_plot_df <- overall_metrics %>%
  dplyr::select(layer, model, ROC_AUC, balanced_accuracy, macro_F1, F1_LUAD, F1_LUSC) %>%
  tidyr::pivot_longer(
    cols      = c(ROC_AUC, balanced_accuracy, macro_F1, F1_LUAD, F1_LUSC),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(
      metric,
      ROC_AUC           = "ROC-AUC",
      balanced_accuracy = "Balanced accuracy",
      macro_F1          = "Macro-F1",
      F1_LUAD           = "F1 LUAD",
      F1_LUSC           = "F1 LUSC"
    )
  )

p_metrics <- ggplot2::ggplot(
  metric_plot_df,
  ggplot2::aes(x = layer, y = value, fill = model)
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.72),
    width = 0.62, alpha = 0.92
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = round(value, 3)),
    position = ggplot2::position_dodge(width = 0.72),
    vjust = -0.25, size = 3.2, fontface = "bold"
  ) +
  ggplot2::facet_wrap(~metric, ncol = 2) +
  ggplot2::scale_fill_manual(values = model_colors) +
  ggplot2::scale_y_continuous(limits = c(0, 1.08)) +
  ggplot2::labs(
    title    = "Single-omics model comparison",
    subtitle = "Performance from out-of-fold predictions",
    x        = "Omics layer",
    y        = "Metric value",
    fill     = "Model"
  ) +
  theme_model(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

save_plot(p_metrics, "02_single_omics_metric_comparison.png", width = 11, height = 8)

p_auc <- overall_metrics %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x    = reorder(paste(layer, model, sep = " + "), ROC_AUC),
      y    = ROC_AUC,
      fill = model
    )
  ) +
  ggplot2::geom_col(width = 0.7, alpha = 0.92) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = model_colors) +
  ggplot2::geom_text(
    ggplot2::aes(label = round(ROC_AUC, 3)),
    hjust = -0.1, fontface = "bold", size = 4
  ) +
  ggplot2::scale_y_continuous(limits = c(0, 1.08)) +
  ggplot2::labs(
    title    = "ROC-AUC ranking of single-omics models",
    subtitle = "Higher values indicate better LUAD/LUSC discrimination",
    x        = "Layer + model",
    y        = "ROC-AUC",
    fill     = "Model"
  ) +
  theme_model()

save_plot(p_auc, "03_single_omics_auc_ranking.png", width = 9, height = 6)


# -----------------------------
# 11. Confusion matrix visualisation
# -----------------------------

p_confusion <- ggplot2::ggplot(
  confusion_counts,
  ggplot2::aes(x = predicted, y = actual, fill = n)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.8) +
  ggplot2::geom_text(ggplot2::aes(label = n), size = 5, fontface = "bold", color = "white") +
  ggplot2::facet_grid(layer ~ model) +
  ggplot2::scale_fill_viridis_c(option = "magma", direction = -1) +
  ggplot2::labs(
    title    = "Single-omics confusion matrices",
    subtitle = "Counts are based on out-of-fold predictions",
    x        = "Predicted subtype",
    y        = "Actual subtype",
    fill     = "Count"
  ) +
  theme_model(base_size = 12)

save_plot(p_confusion, "04_single_omics_confusion_matrices.png", width = 9, height = 9)


# -----------------------------
# 12. Misclassification visualisation
# -----------------------------

misclassification_df <- predictions_all %>%
  dplyr::mutate(correct = actual == predicted) %>%
  dplyr::count(layer, model, correct, name = "n") %>%
  dplyr::mutate(result = ifelse(correct, "Correct", "Misclassified"))

p_misclass <- ggplot2::ggplot(
  misclassification_df,
  ggplot2::aes(x = layer, y = n, fill = result)
) +
  ggplot2::geom_col(width = 0.7, alpha = 0.92) +
  ggplot2::facet_wrap(~model) +
  ggplot2::scale_fill_manual(
    values = c("Correct" = "#2CA02C", "Misclassified" = "#D62728")
  ) +
  ggplot2::labs(
    title    = "Correct and misclassified samples by model",
    subtitle = "Based on out-of-fold predictions",
    x        = "Omics layer",
    y        = "Number of samples",
    fill     = "Prediction result"
  ) +
  theme_model() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

save_plot(p_misclass, "05_single_omics_misclassification_counts.png", width = 9, height = 5.5)


# -----------------------------
# 13. Feature importance plots
# -----------------------------

plot_elastic_features <- function(layer_name, top_n = 20) {
  df <- elastic_feature_summary %>%
    dplyr::filter(layer == layer_name) %>%
    dplyr::slice_max(order_by = mean_abs_coefficient, n = top_n, with_ties = FALSE) %>%
    dplyr::mutate(feature = reorder(feature, mean_abs_coefficient))
  
  if (nrow(df) == 0) {
    message("No Elastic Net features to plot for ", layer_name)
    return(NULL)
  }
  
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = feature, y = mean_abs_coefficient, fill = selected_in_folds)
  ) +
    ggplot2::geom_col(width = 0.75, alpha = 0.92) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_viridis_c(option = "plasma") +
    ggplot2::labs(
      title    = paste0("Top Elastic Net features: ", layer_name),
      subtitle = "Ranked by mean absolute coefficient across CV folds",
      x        = "Feature",
      y        = "Mean absolute coefficient",
      fill     = "Selected\nfolds"
    ) +
    theme_model(base_size = 11)
  
  save_plot(
    p,
    paste0("06_top_elastic_net_features_", safe_filename(layer_name), ".png"),
    width = 9, height = 7
  )
  p
}

plot_rf_features <- function(layer_name, top_n = 20) {
  df <- rf_feature_summary %>%
    dplyr::filter(layer == layer_name) %>%
    dplyr::slice_max(order_by = mean_importance, n = top_n, with_ties = FALSE) %>%
    dplyr::mutate(feature = reorder(feature, mean_importance))
  
  if (nrow(df) == 0) {
    message("No Random Forest features to plot for ", layer_name)
    return(NULL)
  }
  
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = feature, y = mean_importance, fill = selected_in_folds)
  ) +
    ggplot2::geom_col(width = 0.75, alpha = 0.92) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_viridis_c(option = "viridis") +
    ggplot2::labs(
      title    = paste0("Top Random Forest features: ", layer_name),
      subtitle = "Ranked by mean impurity importance across CV folds",
      x        = "Feature",
      y        = "Mean importance",
      fill     = "Selected\nfolds"
    ) +
    theme_model(base_size = 11)
  
  save_plot(
    p,
    paste0("07_top_random_forest_features_", safe_filename(layer_name), ".png"),
    width = 9, height = 7
  )
  p
}

for (layer_name in names(omics_list)) {
  plot_elastic_features(layer_name, top_n = 20)
  plot_rf_features(layer_name, top_n = 20)
}


# -----------------------------
# 14. Final interpretation table
# -----------------------------

best_model <- overall_metrics %>%
  dplyr::arrange(dplyr::desc(ROC_AUC)) %>%
  dplyr::slice(1)

single_omics_interpretation <- tibble::tibble(
  section = c(
    "Single-omics modelling",
    "Leakage control",
    "Best baseline model",
    "Why report macro-F1",
    "Next analysis step"
  ),
  conclusion = c(
    "Elastic Net and Random Forest were trained separately on mRNA, DNA methylation and RPPA to establish single-omics baseline performance.",
    "Median imputation, scaling and top-variable filtering were fitted inside each cross-validation training fold only, then applied to the held-out fold.",
    paste0(
      "The strongest single-omics baseline was ", as.character(best_model$layer),
      " + ", as.character(best_model$model),
      " with ROC-AUC = ", round(best_model$ROC_AUC, 3), "."
    ),
    "LUSC is used as the positive class for probability-based ROC calculations, but F1_LUAD, F1_LUSC and macro-F1 are reported to evaluate both subtypes fairly.",
    "The next step is multi-omics integration to test whether combining omics layers improves over the best single-omics baseline."
  )
)

readr::write_tsv(
  single_omics_interpretation,
  file.path(tables_dir, "single_omics_interpretation_summary.tsv")
)

message("\nSingle-omics interpretation summary:")
print(single_omics_interpretation)


# -----------------------------
# 15. Save result object
# -----------------------------

single_omics_results <- list(
  model_plan              = model_plan,
  predictions             = predictions_all,
  fold_metrics            = fold_metrics_all,
  overall_metrics         = overall_metrics,
  metric_summary          = metric_summary,
  confusion_counts        = confusion_counts,
  roc_data                = roc_data,
  elastic_features_by_fold = elastic_features_all,
  elastic_feature_summary = elastic_feature_summary,
  rf_features_by_fold     = rf_features_all,
  rf_feature_summary      = rf_feature_summary,
  interpretation          = single_omics_interpretation,
  note = paste(
    "Single-omics models were evaluated using 5-fold stratified outer cross-validation.",
    "All imputation, scaling and top-feature filtering were fitted inside training folds only.",
    "No clinical/source-site metadata variables were used as predictors.",
    "LUSC was used as the positive class for probability-based metrics,",
    "but metrics for both LUAD and LUSC were reported."
  )
)

saveRDS(
  single_omics_results,
  file.path(processed_dir, "04_single_omics_results.rds"),
  compress = FALSE
)

message("\nSaved single-omics modelling results to:")
message(file.path(processed_dir, "04_single_omics_results.rds"))
message("\nSingle-omics figures saved in:")
message(figures_single_dir)
message("\n04_single_omics_models.R completed successfully.")