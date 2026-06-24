############################################################
# 04_single_omics_models.R
# Single-omics baseline models:
# Elastic net logistic regression + Random forest
############################################################

# -----------------------------
# 1. Setup
# -----------------------------

source("R/00_setup.R")

model_packages <- c(
  "glmnet",
  "ranger",
  "pROC"
)

installed <- model_packages %in% rownames(installed.packages())

if (any(!installed)) {
  install.packages(model_packages[!installed])
}

suppressPackageStartupMessages({
  library(glmnet)
  library(ranger)
  library(pROC)
})

# -----------------------------
# 2. User settings
# -----------------------------

set.seed(123)

K_OUTER <- 5

RUN_ELASTIC_NET <- TRUE
RUN_RANDOM_FOREST <- TRUE

# For first test run, you can use 100.
# For final results, use 300 or 500.
RF_NUM_TREES <- 300

RF_NUM_THREADS <- max(1, parallel::detectCores() - 1)

message("Using ", RF_NUM_THREADS, " thread(s) for random forest.")

# -----------------------------
# 3. Load preprocessed data
# -----------------------------

processed_data <- readRDS(
  file.path(processed_dir, "02_preprocessed_data.rds")
)

labels <- processed_data$labels
mrna <- processed_data$mrna
methylation <- processed_data$methylation
rppa <- processed_data$rppa

stopifnot(identical(labels$sample_id, rownames(mrna)))
stopifnot(identical(labels$sample_id, rownames(methylation)))
stopifnot(identical(labels$sample_id, rownames(rppa)))

labels$subtype <- factor(labels$subtype, levels = c("LUAD", "LUSC"))

message("Loaded preprocessed data for single-omics modelling.")

# -----------------------------
# 4. Helper: make feature names safe
# -----------------------------
# Some RPPA or gene names may contain characters that can confuse model functions.
# We keep a map from safe feature names back to original feature names.

make_safe_matrix <- function(mat, layer_name) {
  original_names <- colnames(mat)
  safe_names <- make.names(original_names, unique = TRUE)
  
  colnames(mat) <- safe_names
  
  feature_map <- tibble::tibble(
    layer = layer_name,
    feature_safe = safe_names,
    feature_original = original_names
  )
  
  list(matrix = mat, feature_map = feature_map)
}

mrna_safe <- make_safe_matrix(mrna, "mRNA")
methylation_safe <- make_safe_matrix(methylation, "DNA methylation")
rppa_safe <- make_safe_matrix(rppa, "RPPA")

omics_list <- list(
  "mRNA" = mrna_safe$matrix,
  "DNA methylation" = methylation_safe$matrix,
  "RPPA" = rppa_safe$matrix
)

feature_maps <- dplyr::bind_rows(
  mrna_safe$feature_map,
  methylation_safe$feature_map,
  rppa_safe$feature_map
)

readr::write_tsv(
  feature_maps,
  file.path(tables_dir, "single_omics_feature_name_map.tsv")
)

# -----------------------------
# 5. Helper: stratified folds
# -----------------------------

make_stratified_folds <- function(y, k = 5, seed = 123) {
  set.seed(seed)
  
  y <- factor(y)
  fold_id <- rep(NA_integer_, length(y))
  
  for (class_name in levels(y)) {
    idx <- which(y == class_name)
    idx <- sample(idx)
    
    class_folds <- rep(seq_len(k), length.out = length(idx))
    fold_id[idx] <- class_folds
  }
  
  return(fold_id)
}

outer_folds <- make_stratified_folds(
  labels$subtype,
  k = K_OUTER,
  seed = 123
)

fold_balance <- tibble::tibble(
  sample_id = labels$sample_id,
  subtype = labels$subtype,
  outer_fold = outer_folds
) %>%
  dplyr::count(outer_fold, subtype)

message("\nOuter fold class balance:")
print(fold_balance)

readr::write_tsv(
  fold_balance,
  file.path(tables_dir, "single_omics_outer_fold_balance.tsv")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
)

# -----------------------------
# 6. Helper: metrics
# -----------------------------

safe_divide <- function(a, b) {
  ifelse(b == 0, NA_real_, a / b)
}

compute_binary_metrics <- function(actual, predicted, probability_lusc) {
  actual <- factor(actual, levels = c("LUAD", "LUSC"))
  predicted <- factor(predicted, levels = c("LUAD", "LUSC"))
  
  cm <- table(actual, predicted)
  
  tn <- cm["LUAD", "LUAD"]
  fp <- cm["LUAD", "LUSC"]
  fn <- cm["LUSC", "LUAD"]
  tp <- cm["LUSC", "LUSC"]
  
  sensitivity_lusc <- safe_divide(tp, tp + fn)
  specificity_luad <- safe_divide(tn, tn + fp)
  
  precision_lusc <- safe_divide(tp, tp + fp)
  recall_lusc <- sensitivity_lusc
  
  f1_lusc <- safe_divide(
    2 * precision_lusc * recall_lusc,
    precision_lusc + recall_lusc
  )
  
  accuracy <- safe_divide(tp + tn, tp + tn + fp + fn)
  balanced_accuracy <- mean(c(sensitivity_lusc, specificity_luad), na.rm = TRUE)
  
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
    F1_LUSC = f1_lusc,
    sensitivity_LUSC = sensitivity_lusc,
    specificity_LUAD = specificity_luad
  )
}

make_confusion_table <- function(predictions_df) {
  predictions_df %>%
    dplyr::group_by(layer, model) %>%
    dplyr::group_modify(~ {
      cm <- as.data.frame(
        table(
          actual = factor(.x$actual, levels = c("LUAD", "LUSC")),
          predicted = factor(.x$predicted, levels = c("LUAD", "LUSC"))
        )
      )
      
      tibble::as_tibble(cm)
    }) %>%
    dplyr::ungroup()
}

# -----------------------------
# 7. Elastic net model function
# -----------------------------

run_elastic_net_cv <- function(x, y, sample_ids, layer_name, outer_folds) {
  message("\nRunning elastic net for ", layer_name)
  
  all_predictions <- list()
  all_features <- list()
  
  for (fold in sort(unique(outer_folds))) {
    message("  Outer fold ", fold)
    
    test_idx <- which(outer_folds == fold)
    train_idx <- setdiff(seq_along(y), test_idx)
    
    x_train <- x[train_idx, , drop = FALSE]
    x_test <- x[test_idx, , drop = FALSE]
    
    y_train <- y[train_idx]
    y_test <- y[test_idx]
    
    y_train_binary <- ifelse(y_train == "LUSC", 1, 0)
    
    inner_foldid <- make_stratified_folds(
      y_train,
      k = 5,
      seed = 100 + fold
    )
    
    cv_fit <- glmnet::cv.glmnet(
      x = x_train,
      y = y_train_binary,
      family = "binomial",
      alpha = 0.5,
      type.measure = "auc",
      foldid = inner_foldid,
      standardize = TRUE
    )
    
    probability_lusc <- as.numeric(
      predict(
        cv_fit,
        newx = x_test,
        s = "lambda.min",
        type = "response"
      )
    )
    
    predicted <- ifelse(probability_lusc >= 0.5, "LUSC", "LUAD")
    
    pred_df <- tibble::tibble(
      sample_id = sample_ids[test_idx],
      actual = as.character(y_test),
      predicted = predicted,
      probability_LUSC = probability_lusc,
      outer_fold = fold,
      layer = layer_name,
      model = "Elastic net",
      lambda_min = cv_fit$lambda.min,
      lambda_1se = cv_fit$lambda.1se
    )
    
    coefs <- as.matrix(
      coef(cv_fit, s = "lambda.min")
    )
    
    coef_df <- tibble::tibble(
      feature_safe = rownames(coefs),
      coefficient = as.numeric(coefs[, 1])
    ) %>%
      dplyr::filter(feature_safe != "(Intercept)") %>%
      dplyr::filter(coefficient != 0) %>%
      dplyr::mutate(
        layer = layer_name,
        model = "Elastic net",
        outer_fold = fold,
        abs_coefficient = abs(coefficient)
      ) %>%
      dplyr::arrange(dplyr::desc(abs_coefficient))
    
    all_predictions[[as.character(fold)]] <- pred_df
    all_features[[as.character(fold)]] <- coef_df
  }
  
  predictions <- dplyr::bind_rows(all_predictions)
  selected_features <- dplyr::bind_rows(all_features)
  
  return(
    list(
      predictions = predictions,
      selected_features = selected_features
    )
  )
}

# -----------------------------
# 8. Random forest model function
# -----------------------------

run_random_forest_cv <- function(x, y, sample_ids, layer_name, outer_folds) {
  message("\nRunning random forest for ", layer_name)
  
  all_predictions <- list()
  all_importance <- list()
  
  for (fold in sort(unique(outer_folds))) {
    message("  Outer fold ", fold)
    
    test_idx <- which(outer_folds == fold)
    train_idx <- setdiff(seq_along(y), test_idx)
    
    x_train <- x[train_idx, , drop = FALSE]
    x_test <- x[test_idx, , drop = FALSE]
    
    y_train <- factor(y[train_idx], levels = c("LUAD", "LUSC"))
    y_test <- factor(y[test_idx], levels = c("LUAD", "LUSC"))
    
    x_train_df <- as.data.frame(x_train, check.names = FALSE)
    x_test_df <- as.data.frame(x_test, check.names = FALSE)
    
    mtry_value <- max(1, floor(sqrt(ncol(x_train))))
    
    rf_fit <- ranger::ranger(
      x = x_train_df,
      y = y_train,
      probability = TRUE,
      num.trees = RF_NUM_TREES,
      mtry = mtry_value,
      importance = "permutation",
      num.threads = RF_NUM_THREADS,
      seed = 200 + fold
    )
    
    rf_pred <- predict(
      rf_fit,
      data = x_test_df
    )
    
    probability_lusc <- rf_pred$predictions[, "LUSC"]
    predicted <- ifelse(probability_lusc >= 0.5, "LUSC", "LUAD")
    
    pred_df <- tibble::tibble(
      sample_id = sample_ids[test_idx],
      actual = as.character(y_test),
      predicted = predicted,
      probability_LUSC = probability_lusc,
      outer_fold = fold,
      layer = layer_name,
      model = "Random forest",
      mtry = mtry_value,
      num_trees = RF_NUM_TREES
    )
    
    importance_values <- ranger::importance(rf_fit)
    
    importance_df <- tibble::tibble(
      feature_safe = names(importance_values),
      importance = as.numeric(importance_values)
    ) %>%
      dplyr::mutate(
        layer = layer_name,
        model = "Random forest",
        outer_fold = fold
      ) %>%
      dplyr::arrange(dplyr::desc(importance)) %>%
      dplyr::slice_head(n = 100)
    
    all_predictions[[as.character(fold)]] <- pred_df
    all_importance[[as.character(fold)]] <- importance_df
  }
  
  predictions <- dplyr::bind_rows(all_predictions)
  importance <- dplyr::bind_rows(all_importance)
  
  return(
    list(
      predictions = predictions,
      importance = importance
    )
  )
}

# -----------------------------
# 9. Run models
# -----------------------------

all_predictions <- list()
elastic_net_features <- list()
random_forest_importance <- list()

y <- labels$subtype
sample_ids <- labels$sample_id

for (layer_name in names(omics_list)) {
  x <- omics_list[[layer_name]]
  
  message("\n=======================================")
  message("Layer: ", layer_name)
  message("Input dimensions: ", nrow(x), " samples x ", ncol(x), " features")
  message("=======================================")
  
  if (RUN_ELASTIC_NET) {
    en_result <- run_elastic_net_cv(
      x = x,
      y = y,
      sample_ids = sample_ids,
      layer_name = layer_name,
      outer_folds = outer_folds
    )
    
    all_predictions[[paste(layer_name, "Elastic net")]] <- en_result$predictions
    elastic_net_features[[layer_name]] <- en_result$selected_features
  }
  
  if (RUN_RANDOM_FOREST) {
    rf_result <- run_random_forest_cv(
      x = x,
      y = y,
      sample_ids = sample_ids,
      layer_name = layer_name,
      outer_folds = outer_folds
    )
    
    all_predictions[[paste(layer_name, "Random forest")]] <- rf_result$predictions
    random_forest_importance[[layer_name]] <- rf_result$importance
  }
}

single_omics_predictions <- dplyr::bind_rows(all_predictions)

single_omics_predictions <- single_omics_predictions %>%
  dplyr::mutate(
    actual = factor(actual, levels = c("LUAD", "LUSC")),
    predicted = factor(predicted, levels = c("LUAD", "LUSC"))
  )

# -----------------------------
# 10. Metrics
# -----------------------------

fold_metrics <- single_omics_predictions %>%
  dplyr::group_by(layer, model, outer_fold) %>%
  dplyr::group_modify(~ {
    compute_binary_metrics(
      actual = .x$actual,
      predicted = .x$predicted,
      probability_lusc = .x$probability_LUSC
    )
  }) %>%
  dplyr::ungroup()

overall_metrics <- single_omics_predictions %>%
  dplyr::group_by(layer, model) %>%
  dplyr::group_modify(~ {
    compute_binary_metrics(
      actual = .x$actual,
      predicted = .x$predicted,
      probability_lusc = .x$probability_LUSC
    )
  }) %>%
  dplyr::ungroup()

metric_summary <- fold_metrics %>%
  dplyr::group_by(layer, model) %>%
  dplyr::summarise(
    mean_ROC_AUC = mean(ROC_AUC, na.rm = TRUE),
    sd_ROC_AUC = sd(ROC_AUC, na.rm = TRUE),
    mean_balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    sd_balanced_accuracy = sd(balanced_accuracy, na.rm = TRUE),
    mean_F1_LUSC = mean(F1_LUSC, na.rm = TRUE),
    sd_F1_LUSC = sd(F1_LUSC, na.rm = TRUE),
    .groups = "drop"
  )

message("\nSingle-omics fold-level metric summary:")
print(metric_summary)

message("\nSingle-omics overall metrics:")
print(overall_metrics)

# -----------------------------
# 11. Confusion matrices
# -----------------------------

confusion_matrices <- make_confusion_table(single_omics_predictions)

message("\nConfusion matrices:")
print(confusion_matrices)

# -----------------------------
# 12. Add original feature names to importance tables
# -----------------------------

elastic_net_features_df <- dplyr::bind_rows(elastic_net_features)

if (nrow(elastic_net_features_df) > 0) {
  elastic_net_features_df <- elastic_net_features_df %>%
    dplyr::left_join(
      feature_maps,
      by = c("layer", "feature_safe")
    )
}

random_forest_importance_df <- dplyr::bind_rows(random_forest_importance)

if (nrow(random_forest_importance_df) > 0) {
  random_forest_importance_df <- random_forest_importance_df %>%
    dplyr::left_join(
      feature_maps,
      by = c("layer", "feature_safe")
    )
}

# -----------------------------
# 13. Save tables
# -----------------------------

readr::write_tsv(
  single_omics_predictions,
  file.path(tables_dir, "single_omics_predictions.tsv")
)

readr::write_tsv(
  fold_metrics,
  file.path(tables_dir, "single_omics_fold_metrics.tsv")
)

readr::write_tsv(
  overall_metrics,
  file.path(tables_dir, "single_omics_overall_metrics.tsv")
)

readr::write_tsv(
  metric_summary,
  file.path(tables_dir, "single_omics_metric_summary.tsv")
)

readr::write_tsv(
  confusion_matrices,
  file.path(tables_dir, "single_omics_confusion_matrices.tsv")
)

readr::write_tsv(
  elastic_net_features_df,
  file.path(tables_dir, "single_omics_elastic_net_selected_features.tsv")
)

readr::write_tsv(
  random_forest_importance_df,
  file.path(tables_dir, "single_omics_random_forest_top_importance.tsv")
)

# -----------------------------
# 14. Performance plots
# -----------------------------

p_auc <- ggplot2::ggplot(
  metric_summary,
  ggplot2::aes(x = layer, y = mean_ROC_AUC, fill = model)
) +
  ggplot2::geom_col(position = "dodge", alpha = 0.85) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = mean_ROC_AUC - sd_ROC_AUC,
      ymax = mean_ROC_AUC + sd_ROC_AUC
    ),
    position = ggplot2::position_dodge(width = 0.9),
    width = 0.2
  ) +
  ggplot2::coord_cartesian(ylim = c(0.5, 1.0)) +
  ggplot2::labs(
    title = "Single-omics model performance",
    subtitle = "Mean ROC-AUC across outer folds",
    x = "Omics layer",
    y = "Mean ROC-AUC",
    fill = "Model"
  ) +
  ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(
  file.path(figures_dir, "single_omics_auc_comparison.png"),
  p_auc,
  width = 8,
  height = 5,
  dpi = 300
)

p_balacc <- ggplot2::ggplot(
  metric_summary,
  ggplot2::aes(x = layer, y = mean_balanced_accuracy, fill = model)
) +
  ggplot2::geom_col(position = "dodge", alpha = 0.85) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = mean_balanced_accuracy - sd_balanced_accuracy,
      ymax = mean_balanced_accuracy + sd_balanced_accuracy
    ),
    position = ggplot2::position_dodge(width = 0.9),
    width = 0.2
  ) +
  ggplot2::coord_cartesian(ylim = c(0.5, 1.0)) +
  ggplot2::labs(
    title = "Single-omics model performance",
    subtitle = "Mean balanced accuracy across outer folds",
    x = "Omics layer",
    y = "Mean balanced accuracy",
    fill = "Model"
  ) +
  ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(
  file.path(figures_dir, "single_omics_balanced_accuracy_comparison.png"),
  p_balacc,
  width = 8,
  height = 5,
  dpi = 300
)

print(p_auc)
print(p_balacc)

# -----------------------------
# 15. Save model result object
# -----------------------------

single_omics_results <- list(
  predictions = single_omics_predictions,
  fold_metrics = fold_metrics,
  overall_metrics = overall_metrics,
  metric_summary = metric_summary,
  confusion_matrices = confusion_matrices,
  elastic_net_features = elastic_net_features_df,
  random_forest_importance = random_forest_importance_df,
  outer_folds = outer_folds,
  fold_balance = fold_balance
)

saveRDS(
  single_omics_results,
  file.path(processed_dir, "04_single_omics_results.rds"),
  compress = FALSE
)

message("\nSaved single-omics results to:")
message(file.path(processed_dir, "04_single_omics_results.rds"))

message("\n04_single_omics_models.R completed successfully.")