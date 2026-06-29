############################################################
# 06_diablo_models.R
# DIABLO / multiblock sPLS-DA multi-omics integration
# NSCLC multi-omics project: LUAD vs LUSC
#
# Input:
#   data/processed/02_qc_model_input_unimputed.rds
#   data/processed/04_single_omics_results.rds      optional, for comparison
#   data/processed/05_early_fusion_results.rds      optional, for comparison
#
# Output:
#   data/processed/06_diablo_results.rds
#   results/tables/diablo_*.tsv
#   results/figures/diablo/*.png
#
# Model:
#   - DIABLO: supervised multiblock sparse PLS-DA via mixOmics::block.splsda
#
# Integration strategy:
#   - Keep mRNA, DNA methylation, and RPPA as separate omics blocks
#   - Learn shared supervised latent components across blocks
#   - Select sparse feature sets from each omics layer
#
# Leakage control:
#   Inside each outer CV fold, using training data only:
#     - variance-based prefiltering per block
#     - median imputation per block
#     - scaling per block
#     - optional inner-CV DIABLO keepX tuning
#     - model fitting
#
# The test fold is transformed only using training-fold parameters.
# Full-data final model is created only for visualization/feature interpretation,
# not for reporting predictive performance.
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
  "mixOmics",
  "pROC",
  "readr",
  "tibble",
  "scales",
  "purrr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall missing CRAN packages with install.packages().",
    "\nInstall mixOmics with:",
    "\n  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')",
    "\n  BiocManager::install('mixOmics')"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(matrixStats)
  library(mixOmics)
  library(pROC)
  library(readr)
  library(tibble)
  library(scales)
  library(purrr)
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

figures_diablo_dir <- file.path(figures_dir, "diablo")
dir.create(figures_diablo_dir, recursive = TRUE, showWarnings = FALSE)

message("Loaded QC modelling object for DIABLO.")
message("Samples: ", nrow(labels))
message("mRNA features: ", ncol(mrna))
message("DNA methylation features: ", ncol(methylation))
message("RPPA features: ", ncol(rppa))

# -----------------------------
# 2. DIABLO modelling plan
# -----------------------------

set.seed(123)

outer_k <- 5
inner_k_diablo <- 5
ncomp_diablo <- 2
prediction_distance <- "centroids.dist"
design_strength <- 0.1

# Keep this TRUE for the final project run.
# Set FALSE only if runtime becomes too high; then the fixed keepX below is used.
diablo_run_tuning <- TRUE

# Fold-wise top-variable filtering before DIABLO.
# This is done inside each CV training fold only and is used for computational feasibility.
# It is NOT global supervised feature selection.
diablo_prefilter_plan <- tibble::tibble(
  layer = c("mRNA", "DNA methylation", "RPPA"),
  top_n = c(500, 1000, ncol(rppa))
)

# Small keepX grid for inner-CV DIABLO tuning.
# These values control how many variables DIABLO keeps per component per omics block.
diablo_keepX_grid <- list(
  "mRNA" = c(10, 25, 50, 100),
  "DNA methylation" = c(10, 25, 50, 100),
  "RPPA" = c(5, 10, 25, 50)
)

# Fallback if tuning is disabled or fails in a fold.
diablo_fixed_keepX <- list(
  "mRNA" = rep(50, ncomp_diablo),
  "DNA methylation" = rep(50, ncomp_diablo),
  "RPPA" = rep(25, ncomp_diablo)
)

readr::write_tsv(
  diablo_prefilter_plan,
  file.path(tables_dir, "diablo_prefilter_plan.tsv")
)

readr::write_tsv(
  tibble::tibble(
    layer = rep(names(diablo_keepX_grid), lengths(diablo_keepX_grid)),
    keepX_candidate = unlist(diablo_keepX_grid, use.names = FALSE)
  ),
  file.path(tables_dir, "diablo_keepX_grid.tsv")
)

omics_list <- list(
  "mRNA" = mrna,
  "DNA methylation" = methylation,
  "RPPA" = rppa
)

make_design_matrix <- function(block_names, strength = 0.1) {
  design <- matrix(
    strength,
    nrow = length(block_names),
    ncol = length(block_names),
    dimnames = list(block_names, block_names)
  )
  diag(design) <- 0
  design
}

diablo_design <- make_design_matrix(names(omics_list), design_strength)

readr::write_tsv(
  as.data.frame(diablo_design) %>%
    tibble::rownames_to_column("block"),
  file.path(tables_dir, "diablo_design_matrix.tsv")
)

# -----------------------------
# 3. Plot theme and colours
# -----------------------------

subtype_colors <- c(
  "LUAD" = "#1F77B4",
  "LUSC" = "#D62728"
)

layer_colors <- c(
  "mRNA" = "#1F77B4",
  "DNA methylation" = "#2CA02C",
  "RPPA" = "#9467BD"
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
    filename = file.path(figures_diablo_dir, filename),
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

fit_layer_preprocess_train <- function(x_train,
                                       layer_name,
                                       top_n,
                                       scale_features = TRUE) {
  feature_var <- matrixStats::colVars(x_train, na.rm = TRUE)
  names(feature_var) <- colnames(x_train)

  keep <- !is.na(feature_var) & is.finite(feature_var) & feature_var > 0
  x_train <- x_train[, keep, drop = FALSE]
  feature_var <- feature_var[keep]

  if (ncol(x_train) == 0) {
    stop("No usable features left after variance filtering for layer: ", layer_name)
  }

  n_select <- min(top_n, ncol(x_train))
  selected_features <- names(sort(feature_var, decreasing = TRUE))[seq_len(n_select)]
  x_train <- x_train[, selected_features, drop = FALSE]

  feature_medians <- matrixStats::colMedians(x_train, na.rm = TRUE)
  names(feature_medians) <- colnames(x_train)
  feature_medians[!is.finite(feature_medians)] <- 0

  x_train <- impute_with_medians(x_train, feature_medians)

  feature_sd_after_impute <- matrixStats::colSds(x_train, na.rm = TRUE)
  valid_after_impute <- is.finite(feature_sd_after_impute) & feature_sd_after_impute > 0

  x_train <- x_train[, valid_after_impute, drop = FALSE]
  feature_medians <- feature_medians[valid_after_impute]
  selected_features <- selected_features[valid_after_impute]

  if (scale_features) {
    feature_center <- colMeans(x_train)
    feature_scale <- matrixStats::colSds(x_train)
    feature_scale[!is.finite(feature_scale) | feature_scale == 0] <- 1

    x_train <- sweep(x_train, 2, feature_center, "-")
    x_train <- sweep(x_train, 2, feature_scale, "/")
  } else {
    feature_center <- rep(0, ncol(x_train))
    feature_scale <- rep(1, ncol(x_train))
    names(feature_center) <- colnames(x_train)
    names(feature_scale) <- colnames(x_train)
  }

  list(
    x_train = x_train,
    layer_name = layer_name,
    selected_features = selected_features,
    medians = feature_medians,
    center = feature_center,
    scale = feature_scale,
    scale_features = scale_features
  )
}

apply_layer_preprocess_test <- function(x_test, prep) {
  x_test <- x_test[, prep$selected_features, drop = FALSE]
  x_test <- impute_with_medians(x_test, prep$medians)

  if (isTRUE(prep$scale_features)) {
    x_test <- sweep(x_test, 2, prep$center, "-")
    x_test <- sweep(x_test, 2, prep$scale, "/")
  }

  x_test
}

fit_diablo_preprocess_train <- function(omics_list,
                                        train_idx,
                                        prefilter_plan) {
  prep_list <- list()
  train_blocks <- list()
  feature_counts <- list()

  for (layer_name in names(omics_list)) {
    top_n <- prefilter_plan %>%
      dplyr::filter(layer == layer_name) %>%
      dplyr::pull(top_n)

    x_train_layer <- omics_list[[layer_name]][train_idx, , drop = FALSE]

    prep <- fit_layer_preprocess_train(
      x_train = x_train_layer,
      layer_name = layer_name,
      top_n = top_n,
      scale_features = TRUE
    )

    prep_list[[layer_name]] <- prep
    train_blocks[[layer_name]] <- prep$x_train

    feature_counts[[layer_name]] <- tibble::tibble(
      layer = layer_name,
      n_features_used = ncol(prep$x_train)
    )
  }

  list(
    X_train = train_blocks,
    prep_list = prep_list,
    feature_counts = dplyr::bind_rows(feature_counts)
  )
}

apply_diablo_preprocess_test <- function(omics_list,
                                         test_idx,
                                         diablo_prep) {
  test_blocks <- list()

  for (layer_name in names(diablo_prep$prep_list)) {
    x_test_layer <- omics_list[[layer_name]][test_idx, , drop = FALSE]

    test_blocks[[layer_name]] <- apply_layer_preprocess_test(
      x_test = x_test_layer,
      prep = diablo_prep$prep_list[[layer_name]]
    )
  }

  test_blocks
}

limit_keepX_grid_to_features <- function(keepX_grid, X_train) {
  out <- keepX_grid

  for (layer_name in names(out)) {
    max_features <- ncol(X_train[[layer_name]])
    out[[layer_name]] <- unique(out[[layer_name]][out[[layer_name]] <= max_features])

    if (length(out[[layer_name]]) == 0) {
      out[[layer_name]] <- min(max_features, 5)
    }
  }

  out
}

limit_keepX_to_features <- function(keepX, X_train, ncomp) {
  out <- keepX

  for (layer_name in names(X_train)) {
    max_features <- ncol(X_train[[layer_name]])

    if (is.null(out[[layer_name]])) {
      out[[layer_name]] <- rep(min(25, max_features), ncomp)
    }

    out[[layer_name]] <- pmin(out[[layer_name]], max_features)

    if (length(out[[layer_name]]) < ncomp) {
      out[[layer_name]] <- rep(out[[layer_name]][1], ncomp)
    }

    out[[layer_name]] <- out[[layer_name]][seq_len(ncomp)]
  }

  out
}

tune_diablo_keepX <- function(X_train,
                              y_train,
                              design,
                              keepX_grid,
                              ncomp,
                              seed = 123) {
  keepX_grid <- limit_keepX_grid_to_features(keepX_grid, X_train)

  tune_fit <- mixOmics::tune.block.splsda(
    X = X_train,
    Y = y_train,
    ncomp = ncomp,
    test.keepX = keepX_grid,
    design = design,
    validation = "Mfold",
    folds = inner_k_diablo,
    nrepeat = 1,
    dist = prediction_distance,
    measure = "BER",
    weighted = TRUE,
    scale = FALSE,
    progressBar = FALSE,
    seed = seed
  )

  limit_keepX_to_features(tune_fit$choice.keepX, X_train, ncomp)
}

fit_diablo_model <- function(X_train,
                             y_train,
                             keepX,
                             design,
                             ncomp) {
  keepX <- limit_keepX_to_features(keepX, X_train, ncomp)

  mixOmics::block.splsda(
    X = X_train,
    Y = y_train,
    ncomp = ncomp,
    keepX = keepX,
    design = design,
    scale = FALSE
  )
}

extract_weighted_vote <- function(pred_obj,
                                  ncomp,
                                  distance = "centroids.dist") {
  candidates <- list(
    pred_obj$WeightedVote,
    pred_obj$MajorityVote,
    pred_obj$class
  )

  for (candidate in candidates) {
    if (is.null(candidate)) next

    if (is.list(candidate) && !is.null(candidate[[distance]])) {
      mat <- candidate[[distance]]
      if (is.matrix(mat) || is.data.frame(mat)) {
        return(as.character(mat[, ncomp]))
      }
    }

    if (is.matrix(candidate) || is.data.frame(candidate)) {
      return(as.character(candidate[, ncomp]))
    }
  }

  stop("Could not extract predicted classes from mixOmics predict() output.")
}

extract_lusc_score <- function(pred_obj,
                               ncomp,
                               predicted_classes) {
  candidates <- list(
    pred_obj$WeightedPredict,
    pred_obj$AveragedPredict,
    pred_obj$predict
  )

  for (candidate in candidates) {
    if (is.null(candidate)) next

    if (is.array(candidate) && length(dim(candidate)) == 3) {
      class_names <- dimnames(candidate)[[2]]
      if (!is.null(class_names) && "LUSC" %in% class_names) {
        return(as.numeric(candidate[, "LUSC", ncomp]))
      }
    }

    if (is.matrix(candidate) || is.data.frame(candidate)) {
      if ("LUSC" %in% colnames(candidate)) {
        return(as.numeric(candidate[, "LUSC"]))
      }
    }
  }

  # Fallback: binary score from predicted class.
  # This is less informative for ROC-AUC, but prevents the script from failing if
  # the installed mixOmics version changes the score object layout.
  as.numeric(predicted_classes == "LUSC")
}

compute_binary_metrics <- function(actual, predicted, score_lusc) {
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
      predictor = score_lusc,
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
    predictor = prediction_df$score_lusc,
    levels = c("LUAD", "LUSC"),
    direction = "<",
    quiet = TRUE
  )

  tibble::tibble(
    model_label = "DIABLO multiblock sPLS-DA",
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities
  ) %>%
    dplyr::arrange(FPR, TPR)
}

extract_diablo_loadings <- function(diablo_fit, fold_id = NA_integer_) {
  out <- list()

  for (layer_name in names(diablo_fit$loadings)) {
    loading_mat <- diablo_fit$loadings[[layer_name]]
    max_comp <- min(ncomp_diablo, ncol(loading_mat))

    for (comp_id in seq_len(max_comp)) {
      comp_loading <- loading_mat[, comp_id]

      out[[paste(layer_name, comp_id, sep = "_")]] <- tibble::tibble(
        layer = layer_name,
        component = comp_id,
        feature = rownames(loading_mat),
        loading = as.numeric(comp_loading),
        abs_loading = abs(as.numeric(comp_loading)),
        selected = abs(as.numeric(comp_loading)) > 0,
        fold = fold_id
      ) %>%
        dplyr::filter(selected) %>%
        dplyr::arrange(dplyr::desc(abs_loading))
    }
  }

  dplyr::bind_rows(out)
}

keepX_to_tibble <- function(keepX, fold_id) {
  out <- list()

  for (layer_name in names(keepX)) {
    out[[layer_name]] <- tibble::tibble(
      fold = fold_id,
      layer = layer_name,
      component = seq_along(keepX[[layer_name]]),
      keepX = as.numeric(keepX[[layer_name]])
    )
  }

  dplyr::bind_rows(out)
}

summarise_keepX_for_final_model <- function(keepX_by_fold,
                                            fixed_keepX,
                                            ncomp) {
  if (is.null(keepX_by_fold) || nrow(keepX_by_fold) == 0) {
    return(fixed_keepX)
  }

  out <- list()

  for (layer_name in unique(keepX_by_fold$layer)) {
    out[[layer_name]] <- numeric(ncomp)

    for (comp_id in seq_len(ncomp)) {
      vals <- keepX_by_fold %>%
        dplyr::filter(layer == layer_name, component == comp_id) %>%
        dplyr::pull(keepX)

      if (length(vals) == 0 || all(is.na(vals))) {
        out[[layer_name]][comp_id] <- fixed_keepX[[layer_name]][comp_id]
      } else {
        out[[layer_name]][comp_id] <- as.integer(round(stats::median(vals, na.rm = TRUE)))
      }
    }
  }

  out
}

# -----------------------------
# 5. Run outer-CV DIABLO
# -----------------------------

folds <- create_stratified_folds(
  y = labels$subtype,
  k = outer_k,
  seed = 123
)

all_predictions <- list()
all_fold_metrics <- list()
all_feature_counts <- list()
all_keepX <- list()
all_selected_features <- list()

y <- labels$subtype

for (fold_id in seq_along(folds)) {
  message("\n========================================")
  message("DIABLO | Fold: ", fold_id, "/", outer_k)
  message("========================================")

  test_idx <- folds[[fold_id]]
  train_idx <- setdiff(seq_len(nrow(labels)), test_idx)

  y_train <- y[train_idx]
  y_test <- y[test_idx]

  diablo_prep <- fit_diablo_preprocess_train(
    omics_list = omics_list,
    train_idx = train_idx,
    prefilter_plan = diablo_prefilter_plan
  )

  X_train <- diablo_prep$X_train
  X_test <- apply_diablo_preprocess_test(
    omics_list = omics_list,
    test_idx = test_idx,
    diablo_prep = diablo_prep
  )

  message(
    "Training blocks: ",
    paste(
      paste0(names(X_train), "=", vapply(X_train, ncol, integer(1))),
      collapse = ", "
    )
  )

  chosen_keepX <- diablo_fixed_keepX
  tuning_status <- "fixed_keepX"

  if (isTRUE(diablo_run_tuning)) {
    message("Tuning DIABLO keepX inside training fold...")

    tuning_attempt <- tryCatch(
      {
        tune_diablo_keepX(
          X_train = X_train,
          y_train = y_train,
          design = diablo_design,
          keepX_grid = diablo_keepX_grid,
          ncomp = ncomp_diablo,
          seed = 123 + fold_id
        )
      },
      error = function(e) {
        message("DIABLO tuning failed in fold ", fold_id, ": ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(tuning_attempt)) {
      chosen_keepX <- tuning_attempt
      tuning_status <- "inner_cv_tuned_keepX"
    } else {
      chosen_keepX <- diablo_fixed_keepX
      tuning_status <- "fallback_fixed_keepX"
    }
  }

  chosen_keepX <- limit_keepX_to_features(chosen_keepX, X_train, ncomp_diablo)

  diablo_fit <- fit_diablo_model(
    X_train = X_train,
    y_train = y_train,
    keepX = chosen_keepX,
    design = diablo_design,
    ncomp = ncomp_diablo
  )

  pred_obj <- predict(
    diablo_fit,
    newdata = X_test,
    dist = prediction_distance
  )

  predicted <- extract_weighted_vote(
    pred_obj = pred_obj,
    ncomp = ncomp_diablo,
    distance = prediction_distance
  )

  score_lusc <- extract_lusc_score(
    pred_obj = pred_obj,
    ncomp = ncomp_diablo,
    predicted_classes = predicted
  )

  prediction_df <- tibble::tibble(
    model_label = "DIABLO multiblock sPLS-DA",
    model = "block.splsda",
    integration = "DIABLO",
    fold = fold_id,
    sample_id = labels$sample_id[test_idx],
    actual = as.character(y_test),
    predicted = predicted,
    score_lusc = score_lusc,
    prediction_distance = prediction_distance,
    tuning_status = tuning_status
  )

  fold_metrics <- compute_binary_metrics(
    actual = prediction_df$actual,
    predicted = prediction_df$predicted,
    score_lusc = prediction_df$score_lusc
  ) %>%
    dplyr::mutate(
      model_label = "DIABLO multiblock sPLS-DA",
      model = "block.splsda",
      integration = "DIABLO",
      fold = fold_id,
      ncomp = ncomp_diablo,
      prediction_distance = prediction_distance,
      tuning_status = tuning_status,
      .before = 1
    )

  selected_features <- extract_diablo_loadings(
    diablo_fit = diablo_fit,
    fold_id = fold_id
  ) %>%
    dplyr::mutate(
      model_label = "DIABLO multiblock sPLS-DA",
      integration = "DIABLO"
    )

  feature_counts <- diablo_prep$feature_counts %>%
    dplyr::mutate(
      fold = fold_id,
      model_label = "DIABLO multiblock sPLS-DA"
    )

  all_predictions[[fold_id]] <- prediction_df
  all_fold_metrics[[fold_id]] <- fold_metrics
  all_feature_counts[[fold_id]] <- feature_counts
  all_keepX[[fold_id]] <- keepX_to_tibble(chosen_keepX, fold_id) %>%
    dplyr::mutate(tuning_status = tuning_status)
  all_selected_features[[fold_id]] <- selected_features
}

predictions_all <- dplyr::bind_rows(all_predictions) %>%
  dplyr::mutate(
    actual = factor(actual, levels = c("LUAD", "LUSC")),
    predicted = factor(predicted, levels = c("LUAD", "LUSC"))
  )

fold_metrics_all <- dplyr::bind_rows(all_fold_metrics)
feature_counts_all <- dplyr::bind_rows(all_feature_counts)
keepX_all <- dplyr::bind_rows(all_keepX)
selected_features_all <- dplyr::bind_rows(all_selected_features)

# -----------------------------
# 6. Metrics and summaries
# -----------------------------

overall_metrics <- compute_binary_metrics(
  actual = predictions_all$actual,
  predicted = predictions_all$predicted,
  score_lusc = predictions_all$score_lusc
) %>%
  dplyr::mutate(
    model_label = "DIABLO multiblock sPLS-DA",
    model = "block.splsda",
    integration = "DIABLO",
    n_samples = nrow(predictions_all),
    ncomp = ncomp_diablo,
    prediction_distance = prediction_distance,
    .before = 1
  )

metric_summary <- fold_metrics_all %>%
  dplyr::summarise(
    model_label = "DIABLO multiblock sPLS-DA",
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
    model_label = "DIABLO multiblock sPLS-DA",
    integration = "DIABLO",
    .before = 1
  )

selected_feature_summary <- selected_features_all %>%
  dplyr::group_by(layer, component, feature) %>%
  dplyr::summarise(
    selection_frequency = dplyr::n_distinct(fold),
    mean_abs_loading = mean(abs_loading, na.rm = TRUE),
    median_abs_loading = median(abs_loading, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(layer, component, dplyr::desc(selection_frequency), dplyr::desc(mean_abs_loading))

selected_feature_layer_summary <- selected_feature_summary %>%
  dplyr::group_by(layer, component) %>%
  dplyr::summarise(
    n_unique_selected_features = dplyr::n(),
    n_features_selected_in_all_folds = sum(selection_frequency == outer_k),
    .groups = "drop"
  )

# -----------------------------
# 7. Optional comparison with existing models
# -----------------------------

comparison_metrics <- overall_metrics %>%
  dplyr::select(
    model_label,
    integration,
    ROC_AUC,
    accuracy,
    balanced_accuracy,
    F1_LUAD,
    F1_LUSC,
    macro_F1
  )

single_path <- file.path(processed_dir, "04_single_omics_results.rds")
if (file.exists(single_path)) {
  single_results <- readRDS(single_path)

  if (!is.null(single_results$overall_metrics)) {
    single_compare <- single_results$overall_metrics %>%
      dplyr::mutate(
        model_label = paste(layer, model, sep = " + "),
        integration = "Single omics"
      ) %>%
      dplyr::select(
        model_label,
        integration,
        ROC_AUC,
        accuracy,
        balanced_accuracy,
        F1_LUAD,
        F1_LUSC,
        macro_F1
      )

    comparison_metrics <- dplyr::bind_rows(single_compare, comparison_metrics)
  }
}

early_path <- file.path(processed_dir, "05_early_fusion_results.rds")
if (file.exists(early_path)) {
  early_results <- readRDS(early_path)

  if (!is.null(early_results$overall_metrics)) {
    early_compare <- early_results$overall_metrics %>%
      dplyr::select(
        model_label,
        integration,
        ROC_AUC,
        accuracy,
        balanced_accuracy,
        F1_LUAD,
        F1_LUSC,
        macro_F1
      )

    comparison_metrics <- dplyr::bind_rows(early_compare, comparison_metrics)
  }
}

# -----------------------------
# 8. Final DIABLO model for visualization only
# -----------------------------

message("\nFitting final DIABLO model on all samples for visualization/interpretation only...")

full_idx <- seq_len(nrow(labels))
final_prep <- fit_diablo_preprocess_train(
  omics_list = omics_list,
  train_idx = full_idx,
  prefilter_plan = diablo_prefilter_plan
)

X_full <- final_prep$X_train
y_full <- labels$subtype

final_keepX <- summarise_keepX_for_final_model(
  keepX_by_fold = keepX_all,
  fixed_keepX = diablo_fixed_keepX,
  ncomp = ncomp_diablo
)

final_keepX <- limit_keepX_to_features(final_keepX, X_full, ncomp_diablo)

final_diablo_fit <- fit_diablo_model(
  X_train = X_full,
  y_train = y_full,
  keepX = final_keepX,
  design = diablo_design,
  ncomp = ncomp_diablo
)

final_selected_features <- extract_diablo_loadings(
  diablo_fit = final_diablo_fit,
  fold_id = NA_integer_
) %>%
  dplyr::mutate(
    model_label = "Final DIABLO multiblock sPLS-DA",
    integration = "DIABLO",
    source = "full_data_final_model_for_interpretation_only"
  )

final_scores_by_block <- purrr::imap_dfr(
  final_diablo_fit$variates,
  function(score_mat, layer_name) {
    tibble::tibble(
      sample_id = labels$sample_id,
      subtype = labels$subtype,
      layer = layer_name,
      comp1 = as.numeric(score_mat[, 1]),
      comp2 = if (ncol(score_mat) >= 2) as.numeric(score_mat[, 2]) else NA_real_
    )
  }
)

# Average sample scores across blocks for a compact final visualization.
final_average_scores <- final_scores_by_block %>%
  dplyr::group_by(sample_id, subtype) %>%
  dplyr::summarise(
    comp1 = mean(comp1, na.rm = TRUE),
    comp2 = mean(comp2, na.rm = TRUE),
    .groups = "drop"
  )

final_keepX_table <- keepX_to_tibble(final_keepX, fold_id = NA_integer_) %>%
  dplyr::mutate(source = "median_keepX_from_outer_folds_or_fixed_fallback")

# -----------------------------
# 9. Save tables
# -----------------------------

readr::write_tsv(predictions_all, file.path(tables_dir, "diablo_predictions.tsv"))
readr::write_tsv(fold_metrics_all, file.path(tables_dir, "diablo_fold_metrics.tsv"))
readr::write_tsv(overall_metrics, file.path(tables_dir, "diablo_overall_metrics.tsv"))
readr::write_tsv(metric_summary, file.path(tables_dir, "diablo_metric_summary.tsv"))
readr::write_tsv(confusion_counts, file.path(tables_dir, "diablo_confusion_matrix.tsv"))
readr::write_tsv(feature_counts_all, file.path(tables_dir, "diablo_feature_counts_by_fold.tsv"))
readr::write_tsv(keepX_all, file.path(tables_dir, "diablo_keepX_by_fold.tsv"))
readr::write_tsv(selected_features_all, file.path(tables_dir, "diablo_selected_features_by_fold.tsv"))
readr::write_tsv(selected_feature_summary, file.path(tables_dir, "diablo_selected_feature_summary.tsv"))
readr::write_tsv(selected_feature_layer_summary, file.path(tables_dir, "diablo_selected_feature_layer_summary.tsv"))
readr::write_tsv(final_selected_features, file.path(tables_dir, "diablo_final_selected_features.tsv"))
readr::write_tsv(final_scores_by_block, file.path(tables_dir, "diablo_final_scores_by_block.tsv"))
readr::write_tsv(final_average_scores, file.path(tables_dir, "diablo_final_average_scores.tsv"))
readr::write_tsv(final_keepX_table, file.path(tables_dir, "diablo_final_keepX.tsv"))
readr::write_tsv(comparison_metrics, file.path(tables_dir, "diablo_model_comparison_metrics.tsv"))

# -----------------------------
# 10. Figures
# -----------------------------

roc_df <- make_roc_data(predictions_all)

roc_plot <- ggplot2::ggplot(
  roc_df,
  ggplot2::aes(x = FPR, y = TPR)
) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "DIABLO ROC curve",
    subtitle = paste0(
      "Out-of-fold predictions; ROC-AUC = ",
      round(overall_metrics$ROC_AUC[1], 3)
    ),
    x = "False positive rate",
    y = "True positive rate",
    caption = "LUSC score is extracted from DIABLO prediction scores; predicted class uses weighted vote and centroid distance."
  ) +
  theme_model()

save_plot(roc_plot, "01_diablo_roc_curve.png", width = 7, height = 6)

metric_long <- overall_metrics %>%
  dplyr::select(ROC_AUC, balanced_accuracy, macro_F1, F1_LUAD, F1_LUSC) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = factor(
      metric,
      levels = c("ROC_AUC", "balanced_accuracy", "macro_F1", "F1_LUAD", "F1_LUSC"),
      labels = c("ROC-AUC", "Balanced accuracy", "Macro-F1", "F1 LUAD", "F1 LUSC")
    )
  )

metric_plot <- ggplot2::ggplot(
  metric_long,
  ggplot2::aes(x = metric, y = value)
) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(label = round(value, 3)),
    vjust = -0.35,
    size = 4
  ) +
  ggplot2::ylim(0, 1.05) +
  ggplot2::labs(
    title = "DIABLO cross-validated performance",
    subtitle = "Metrics calculated from outer-fold predictions",
    x = NULL,
    y = "Metric value"
  ) +
  theme_model() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))

save_plot(metric_plot, "02_diablo_metric_summary.png", width = 8, height = 5.5)

confusion_plot <- ggplot2::ggplot(
  confusion_counts,
  ggplot2::aes(x = predicted, y = actual, fill = n)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.8) +
  ggplot2::geom_text(
    ggplot2::aes(label = n),
    size = 6,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_gradient(low = "grey95", high = "grey35") +
  ggplot2::labs(
    title = "DIABLO confusion matrix",
    subtitle = "Out-of-fold predictions from 5-fold cross-validation",
    x = "Predicted subtype",
    y = "Actual subtype",
    fill = "Count"
  ) +
  theme_model()

save_plot(confusion_plot, "03_diablo_confusion_matrix.png", width = 6.5, height = 5.5)

selected_layer_plot_data <- selected_feature_layer_summary %>%
  dplyr::mutate(
    component = paste0("Component ", component)
  )

selected_layer_plot <- ggplot2::ggplot(
  selected_layer_plot_data,
  ggplot2::aes(x = layer, y = n_unique_selected_features, fill = layer)
) +
  ggplot2::geom_col(width = 0.65, show.legend = FALSE) +
  ggplot2::facet_wrap(~component) +
  ggplot2::scale_fill_manual(values = layer_colors) +
  ggplot2::geom_text(
    ggplot2::aes(label = scales::comma(n_unique_selected_features)),
    vjust = -0.35,
    size = 4
  ) +
  ggplot2::labs(
    title = "DIABLO selected features by omics layer",
    subtitle = "Unique selected features across outer CV folds",
    x = NULL,
    y = "Number of unique selected features"
  ) +
  theme_model() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))

save_plot(selected_layer_plot, "04_diablo_selected_features_by_layer.png", width = 8, height = 5.5)

top_features_plot_data <- selected_feature_summary %>%
  dplyr::group_by(layer) %>%
  dplyr::slice_max(
    order_by = selection_frequency * 100000 + mean_abs_loading,
    n = 10,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    feature_label = paste0(feature, " (C", component, ")")
  )

# Reorder manually to avoid adding another package dependency.
top_features_plot_data$feature_label <- factor(
  top_features_plot_data$feature_label,
  levels = unique(top_features_plot_data$feature_label[order(top_features_plot_data$mean_abs_loading)])
)

top_features_plot <- ggplot2::ggplot(
  top_features_plot_data,
  ggplot2::aes(x = feature_label, y = mean_abs_loading, fill = layer)
) +
  ggplot2::geom_col(show.legend = FALSE) +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(~layer, scales = "free_y") +
  ggplot2::scale_fill_manual(values = layer_colors) +
  ggplot2::labs(
    title = "Top DIABLO selected features",
    subtitle = "Ranked by cross-fold selection frequency and mean absolute loading",
    x = NULL,
    y = "Mean absolute loading"
  ) +
  theme_model()

save_plot(top_features_plot, "05_diablo_top_selected_features.png", width = 10, height = 8)

score_plot <- ggplot2::ggplot(
  final_average_scores,
  ggplot2::aes(x = comp1, y = comp2, color = subtype)
) +
  ggplot2::geom_point(size = 2.8, alpha = 0.85) +
  ggplot2::scale_color_manual(values = subtype_colors) +
  ggplot2::labs(
    title = "DIABLO final model component scores",
    subtitle = "Average scores across mRNA, methylation, and RPPA blocks; visualization model only",
    x = "DIABLO component 1",
    y = "DIABLO component 2",
    color = "Subtype"
  ) +
  theme_model()

save_plot(score_plot, "06_diablo_final_component_scores.png", width = 7.5, height = 6)

score_by_block_plot <- ggplot2::ggplot(
  final_scores_by_block,
  ggplot2::aes(x = comp1, y = comp2, color = subtype)
) +
  ggplot2::geom_point(size = 2.2, alpha = 0.75) +
  ggplot2::facet_wrap(~layer, scales = "free") +
  ggplot2::scale_color_manual(values = subtype_colors) +
  ggplot2::labs(
    title = "DIABLO block-specific component scores",
    subtitle = "Final model fitted on all samples for visualization only",
    x = "DIABLO component 1",
    y = "DIABLO component 2",
    color = "Subtype"
  ) +
  theme_model()

save_plot(score_by_block_plot, "07_diablo_block_specific_component_scores.png", width = 10, height = 6)

if (nrow(comparison_metrics) > 1) {
  comparison_plot_data <- comparison_metrics %>%
    dplyr::select(model_label, integration, ROC_AUC, balanced_accuracy, macro_F1) %>%
    tidyr::pivot_longer(
      cols = c(ROC_AUC, balanced_accuracy, macro_F1),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = factor(
        metric,
        levels = c("ROC_AUC", "balanced_accuracy", "macro_F1"),
        labels = c("ROC-AUC", "Balanced accuracy", "Macro-F1")
      ),
      model_label = factor(model_label, levels = unique(model_label[order(integration, model_label)]))
    )

  comparison_plot <- ggplot2::ggplot(
    comparison_plot_data,
    ggplot2::aes(x = model_label, y = value, fill = integration)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::facet_wrap(~metric, ncol = 1) +
    ggplot2::coord_flip() +
    ggplot2::ylim(0, 1.05) +
    ggplot2::labs(
      title = "DIABLO compared with previous models",
      subtitle = "Overall metrics from saved model result objects",
      x = NULL,
      y = "Metric value",
      fill = "Integration type"
    ) +
    theme_model()

  save_plot(comparison_plot, "08_diablo_model_comparison.png", width = 9, height = 9)
}

# -----------------------------
# 11. Save result object
# -----------------------------

diablo_results <- list(
  predictions = predictions_all,
  fold_metrics = fold_metrics_all,
  overall_metrics = overall_metrics,
  metric_summary = metric_summary,
  confusion_counts = confusion_counts,
  feature_counts_by_fold = feature_counts_all,
  keepX_by_fold = keepX_all,
  selected_features_by_fold = selected_features_all,
  selected_feature_summary = selected_feature_summary,
  selected_feature_layer_summary = selected_feature_layer_summary,
  final_selected_features = final_selected_features,
  final_scores_by_block = final_scores_by_block,
  final_average_scores = final_average_scores,
  final_keepX = final_keepX_table,
  comparison_metrics = comparison_metrics,
  modelling_plan = list(
    outer_k = outer_k,
    inner_k_diablo = inner_k_diablo,
    ncomp_diablo = ncomp_diablo,
    prediction_distance = prediction_distance,
    design_strength = design_strength,
    diablo_run_tuning = diablo_run_tuning,
    prefilter_plan = diablo_prefilter_plan,
    keepX_grid = diablo_keepX_grid,
    fixed_keepX = diablo_fixed_keepX,
    design = diablo_design
  ),
  interpretation = list(
    method_summary = paste(
      "DIABLO was implemented as supervised multiblock sparse PLS-DA using mRNA,",
      "DNA methylation and RPPA as separate omics blocks. Unlike early fusion,",
      "the blocks were not concatenated before modelling; DIABLO learned latent",
      "components and sparse feature sets jointly across blocks."
    ),
    leakage_control = paste(
      "Outer-fold performance was estimated with stratified 5-fold cross-validation.",
      "Within each training fold only, each omics block was variance-prefiltered,",
      "median-imputed, scaled, optionally tuned for keepX, and fitted with DIABLO.",
      "The held-out fold was transformed using only the training-fold parameters."
    ),
    batch_note = paste(
      "No site/source batch correction was applied before DIABLO because previous EDA showed",
      "that tissue-source-site variables were confounded with LUAD/LUSC subtype.",
      "Correcting for such variables could remove biological subtype signal."
    ),
    performance_note = paste(
      "Report predictive performance from outer-fold predictions only. The final DIABLO model",
      "fitted on all samples is included only for visualization and selected-feature interpretation."
    )
  )
)

saveRDS(
  diablo_results,
  file.path(processed_dir, "06_diablo_results.rds")
)

message("\nSaved DIABLO result object to: ", file.path(processed_dir, "06_diablo_results.rds"))
message("DIABLO figures saved in: ", figures_diablo_dir)
message("06_diablo_models.R completed successfully.")
