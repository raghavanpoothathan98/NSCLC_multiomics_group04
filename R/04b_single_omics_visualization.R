############################################################
# 04b_single_omics_visualization.R
# Make single-omics results easier to interpret
############################################################

# -----------------------------
# 1. Setup and load results
# -----------------------------

source("R/00_setup.R")

single_omics_results <- readRDS(
  file.path(processed_dir, "04_single_omics_results.rds")
)

metric_summary <- single_omics_results$metric_summary
overall_metrics <- single_omics_results$overall_metrics
confusion_matrices <- single_omics_results$confusion_matrices

message("Loaded single-omics results.")

# -----------------------------
# 2. Colors
# -----------------------------

model_colors <- c(
  "Elastic net" = "#1F77B4",
  "Random forest" = "#D55E00"
)

confusion_colors <- c(
  "Correct" = "#009E73",
  "Wrong" = "#D55E00"
)

# -----------------------------
# 3. Make confusion matrix more readable
# -----------------------------

confusion_readable <- confusion_matrices %>%
  dplyr::mutate(
    actual = as.character(actual),
    predicted = as.character(predicted),
    result = dplyr::if_else(actual == predicted, "Correct", "Wrong")
  ) %>%
  dplyr::group_by(layer, model, actual) %>%
  dplyr::mutate(
    actual_class_total = sum(Freq),
    percent_within_actual = round(Freq / actual_class_total * 100, 1)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    interpretation = dplyr::case_when(
      actual == "LUAD" & predicted == "LUAD" ~ "LUAD correctly classified",
      actual == "LUAD" & predicted == "LUSC" ~ "LUAD misclassified as LUSC",
      actual == "LUSC" & predicted == "LUAD" ~ "LUSC misclassified as LUAD",
      actual == "LUSC" & predicted == "LUSC" ~ "LUSC correctly classified",
      TRUE ~ "Other"
    ),
    label = paste0(Freq, "\n", percent_within_actual, "%")
  )

readr::write_tsv(
  confusion_readable,
  file.path(tables_dir, "single_omics_confusion_readable.tsv")
)

message("\nReadable confusion matrix:")
print(confusion_readable)

# -----------------------------
# 4. Confusion heatmap
# -----------------------------

p_confusion <- ggplot2::ggplot(
  confusion_readable,
  ggplot2::aes(x = predicted, y = actual, fill = result)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.8) +
  ggplot2::geom_text(
    ggplot2::aes(label = label),
    size = 4.2,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(values = confusion_colors) +
  ggplot2::facet_grid(layer ~ model) +
  ggplot2::labs(
    title = "Single-omics confusion matrices",
    subtitle = "Each cell shows count and percentage within the actual class",
    x = "Predicted subtype",
    y = "Actual subtype",
    fill = "Classification"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    axis.text = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  file.path(figures_dir, "single_omics_confusion_heatmaps.png"),
  p_confusion,
  width = 9,
  height = 8,
  dpi = 300
)

print(p_confusion)

# -----------------------------
# 5. Misclassification count summary
# -----------------------------

misclassification_summary <- confusion_readable %>%
  dplyr::filter(result == "Wrong") %>%
  dplyr::group_by(layer, model) %>%
  dplyr::summarise(
    total_errors = sum(Freq),
    .groups = "drop"
  ) %>%
  dplyr::arrange(total_errors)

readr::write_tsv(
  misclassification_summary,
  file.path(tables_dir, "single_omics_misclassification_summary.tsv")
)

message("\nMisclassification summary:")
print(misclassification_summary)

p_errors <- ggplot2::ggplot(
  misclassification_summary,
  ggplot2::aes(x = layer, y = total_errors, fill = model)
) +
  ggplot2::geom_col(position = "dodge", alpha = 0.85) +
  ggplot2::geom_text(
    ggplot2::aes(label = total_errors),
    position = ggplot2::position_dodge(width = 0.9),
    vjust = -0.3,
    size = 4
  ) +
  ggplot2::scale_fill_manual(values = model_colors) +
  ggplot2::labs(
    title = "Total misclassified samples by single-omics model",
    x = "Omics layer",
    y = "Number of incorrect predictions",
    fill = "Model"
  ) +
  ggplot2::theme_minimal(base_size = 12)

ggplot2::ggsave(
  file.path(figures_dir, "single_omics_misclassification_counts.png"),
  p_errors,
  width = 8,
  height = 5,
  dpi = 300
)

print(p_errors)

# -----------------------------
# 6. Make metric table easier to read
# -----------------------------

metrics_readable <- metric_summary %>%
  dplyr::mutate(
    ROC_AUC = paste0(round(mean_ROC_AUC, 3), " ± ", round(sd_ROC_AUC, 3)),
    balanced_accuracy = paste0(
      round(mean_balanced_accuracy, 3),
      " ± ",
      round(sd_balanced_accuracy, 3)
    ),
    F1_LUSC = paste0(round(mean_F1_LUSC, 3), " ± ", round(sd_F1_LUSC, 3))
  ) %>%
  dplyr::select(
    layer,
    model,
    ROC_AUC,
    balanced_accuracy,
    F1_LUSC
  ) %>%
  dplyr::arrange(dplyr::desc(ROC_AUC))

readr::write_tsv(
  metrics_readable,
  file.path(tables_dir, "single_omics_metrics_readable.tsv")
)

message("\nReadable performance table:")
print(metrics_readable)

# -----------------------------
# 7. Clean performance plot
# -----------------------------

metrics_long <- metric_summary %>%
  dplyr::select(
    layer,
    model,
    mean_ROC_AUC,
    mean_balanced_accuracy,
    mean_F1_LUSC
  ) %>%
  tidyr::pivot_longer(
    cols = c(mean_ROC_AUC, mean_balanced_accuracy, mean_F1_LUSC),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::case_when(
      metric == "mean_ROC_AUC" ~ "ROC-AUC",
      metric == "mean_balanced_accuracy" ~ "Balanced accuracy",
      metric == "mean_F1_LUSC" ~ "F1 score for LUSC",
      TRUE ~ metric
    )
  )

p_metrics <- ggplot2::ggplot(
  metrics_long,
  ggplot2::aes(x = layer, y = value, fill = model)
) +
  ggplot2::geom_col(position = "dodge", alpha = 0.85) +
  ggplot2::geom_text(
    ggplot2::aes(label = round(value, 3)),
    position = ggplot2::position_dodge(width = 0.9),
    vjust = -0.25,
    size = 3.5
  ) +
  ggplot2::facet_wrap(~ metric) +
  ggplot2::scale_fill_manual(values = model_colors) +
  ggplot2::coord_cartesian(ylim = c(0.85, 1.0)) +
  ggplot2::labs(
    title = "Single-omics model performance",
    subtitle = "Mean performance across outer cross-validation folds",
    x = "Omics layer",
    y = "Score",
    fill = "Model"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    strip.text = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  file.path(figures_dir, "single_omics_metrics_readable.png"),
  p_metrics,
  width = 10,
  height = 5.5,
  dpi = 300
)

print(p_metrics)

message("\n04b_single_omics_visualization.R completed successfully.")