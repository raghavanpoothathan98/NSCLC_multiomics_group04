############################################################
# 01_load_data.R
# Load already matched NSCLC multi-omics data
############################################################

# -----------------------------
# 1. Setup
# -----------------------------

source("R/00_setup.R")

# -----------------------------
# 2. Load labels
# -----------------------------

message("\nLoading sample labels from:")
message(paths$labels)

labels <- data.table::fread(
  paths$labels,
  data.table = FALSE,
  check.names = FALSE
)

message("\nLabel file columns:")
print(colnames(labels))

# Rename first column to sample_id if needed
if (!"sample_id" %in% colnames(labels)) {
  colnames(labels)[1] <- "sample_id"
}

# Detect subtype column
possible_subtype_cols <- c(
  "subtype",
  "Subtype",
  "cancer_type",
  "Cancer Type",
  "histological_type",
  "Histological Type",
  "label",
  "Label",
  "class",
  "Class",
  "tumor_type",
  "Tumor Type"
)

subtype_col <- intersect(possible_subtype_cols, colnames(labels))

if (length(subtype_col) == 0) {
  stop(
    "Could not detect subtype column. Rename the LUAD/LUSC column in sample_labels.tsv to 'subtype'."
  )
}

subtype_col <- subtype_col[1]

labels <- labels %>%
  dplyr::rename(subtype = dplyr::all_of(subtype_col)) %>%
  dplyr::select(sample_id, subtype)

labels$sample_id <- as.character(labels$sample_id)
labels$subtype <- as.character(labels$subtype)

# Standardize subtype names
labels$subtype <- dplyr::case_when(
  grepl("LUAD|adenocarcinoma", labels$subtype, ignore.case = TRUE) ~ "LUAD",
  grepl("LUSC|squamous", labels$subtype, ignore.case = TRUE) ~ "LUSC",
  TRUE ~ labels$subtype
)

labels$subtype <- factor(labels$subtype, levels = c("LUAD", "LUSC"))

message("\nLoaded labels:")
print(head(labels))

# -----------------------------
# 3. Helper function: read feature x sample matrix
# -----------------------------
# Your matched files have:
# rows = features
# columns = samples
#
# We transpose them to:
# rows = samples
# columns = features

read_feature_by_sample_matrix <- function(path, layer_name, reference_samples) {
  message("\nReading ", layer_name, " data from:")
  message(path)
  
  df <- data.table::fread(
    path,
    data.table = FALSE,
    check.names = FALSE,
    showProgress = TRUE
  )
  
  message(layer_name, " raw table dimensions:")
  print(dim(df))
  
  feature_ids <- as.character(df[[1]])
  sample_ids <- colnames(df)[-1]
  
  # Convert feature x sample table to numeric matrix
  mat <- as.matrix(df[, -1, drop = FALSE])
  mode(mat) <- "numeric"
  
  rownames(mat) <- make.unique(feature_ids)
  colnames(mat) <- sample_ids
  
  # Transpose: samples x features
  mat <- t(mat)
  
  message(layer_name, " matrix after transpose:")
  message(nrow(mat), " samples x ", ncol(mat), " features")
  
  # Check overlap with labels
  overlap <- intersect(rownames(mat), reference_samples)
  
  message(layer_name, " samples overlapping with labels: ", length(overlap))
  
  if (length(overlap) == 0) {
    stop(
      layer_name,
      " has zero sample overlap with labels. Check sample ID format."
    )
  }
  
  return(mat)
}

# -----------------------------
# 4. Load omics matrices
# -----------------------------

mrna <- read_feature_by_sample_matrix(
  paths$mrna,
  "mRNA",
  labels$sample_id
)

methylation <- read_feature_by_sample_matrix(
  paths$methylation,
  "DNA methylation",
  labels$sample_id
)

rppa <- read_feature_by_sample_matrix(
  paths$rppa,
  "RPPA",
  labels$sample_id
)

# -----------------------------
# 5. Keep common samples only
# -----------------------------

common_samples <- Reduce(
  intersect,
  list(
    labels$sample_id,
    rownames(mrna),
    rownames(methylation),
    rownames(rppa)
  )
)

message("\nCommon samples across labels and all omics layers:")
message(length(common_samples))

if (length(common_samples) == 0) {
  stop("No common samples found across labels, mRNA, methylation and RPPA.")
}

# Keep labels order
labels <- labels %>%
  dplyr::filter(sample_id %in% common_samples)

common_samples <- labels$sample_id

mrna <- mrna[common_samples, , drop = FALSE]
methylation <- methylation[common_samples, , drop = FALSE]
rppa <- rppa[common_samples, , drop = FALSE]

# Final order checks
stopifnot(identical(labels$sample_id, rownames(mrna)))
stopifnot(identical(labels$sample_id, rownames(methylation)))
stopifnot(identical(labels$sample_id, rownames(rppa)))

message("All omics layers and labels now have identical sample order.")

# -----------------------------
# 6. Dimension summary
# -----------------------------

dimension_summary <- tibble::tibble(
  layer = c("mRNA", "DNA methylation", "RPPA"),
  samples = c(nrow(mrna), nrow(methylation), nrow(rppa)),
  features = c(ncol(mrna), ncol(methylation), ncol(rppa))
)

message("\nDimension summary:")
print(dimension_summary)

readr::write_tsv(
  dimension_summary,
  file.path(tables_dir, "dimension_summary.tsv")
)

# -----------------------------
# 7. Class balance
# -----------------------------

class_balance <- labels %>%
  dplyr::count(subtype, name = "n") %>%
  dplyr::mutate(percent = round(n / sum(n) * 100, 2))

message("\nClass balance:")
print(class_balance)

readr::write_tsv(
  class_balance,
  file.path(tables_dir, "class_balance.tsv")
)

# -----------------------------
# 8. Basic sanity checks
# -----------------------------

if (any(is.na(labels$subtype))) {
  stop("Some subtype labels are NA. Check sample_labels.tsv.")
}

if (length(unique(labels$subtype)) != 2) {
  warning("Expected two classes: LUAD and LUSC. Found:")
  print(unique(labels$subtype))
}

if (nrow(labels) != nrow(mrna)) {
  stop("Number of labels does not match number of samples.")
}

# -----------------------------
# 9. Save loaded object
# -----------------------------
# This avoids reading the large methylation TSV again.

loaded_data <- list(
  labels = labels,
  mrna = mrna,
  methylation = methylation,
  rppa = rppa,
  dimension_summary = dimension_summary,
  class_balance = class_balance
)

saveRDS(
  loaded_data,
  file.path(processed_dir, "01_loaded_matched_data.rds"),
  compress = FALSE
)

message("\nSaved loaded data to:")
message(file.path(processed_dir, "01_loaded_matched_data.rds"))

message("\n01_load_data.R completed successfully.")