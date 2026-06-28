############################################################
# 01_load_data.R
# Load raw LUAD and LUSC multi-omics data and create matched data
#
# Input:
#   data/raw/LUAD/luad_tcga_pan_can_atlas_2018/
#   data/raw/LUSC/lusc_tcga_pan_can_atlas_2018/
#
# Output:
#   data/processed/01_loaded_matched_data.rds
#
# This script does:
#   - detect raw cBioPortal files
#   - force DNA methylation to HM450 only
#   - read clinical sample metadata
#   - create LUAD/LUSC labels
#   - read mRNA, DNA methylation and RPPA matrices
#   - standardize sample IDs
#   - match samples across labels, metadata and all omics layers
#   - save one official matched object
#
# This script does NOT do:
#   - imputation
#   - scaling
#   - feature selection
#   - batch correction
#   - model training
############################################################

# -----------------------------
# 1. Setup
# -----------------------------

source("R/00_setup.R")

# Keep FALSE to avoid confusion with old matched TSV files.
save_matched_tsv <- FALSE

loaded_rds_path <- file.path(processed_dir, "01_loaded_matched_data.rds")

# -----------------------------
# 2. Helper functions
# -----------------------------

list_raw_files <- function(folder) {
  files <- list.files(
    folder,
    recursive = TRUE,
    full.names = TRUE
  )
  
  files[file.exists(files) & !dir.exists(files)]
}

standardize_sample_id <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.", "-", x)
  x <- trimws(x)
  
  # TCGA sample-level barcode:
  # TCGA-XX-XXXX-XX
  is_tcga <- grepl("^TCGA-", x, ignore.case = TRUE)
  x[is_tcga] <- substr(x[is_tcga], 1, 15)
  
  x
}

read_cbio_table <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  
  on.exit(close(con), add = TRUE)
  
  lines <- readLines(con, n = 200, warn = FALSE)
  
  header_line <- which(
    !grepl("^#", lines) &
      nzchar(trimws(lines))
  )[1]
  
  if (is.na(header_line)) {
    stop("Could not find a valid header line in: ", path)
  }
  
  data.table::fread(
    path,
    data.table = FALSE,
    check.names = FALSE,
    skip = header_line - 1,
    showProgress = TRUE
  )
}

find_file <- function(files,
                      exact_names = character(0),
                      include_regex = NULL,
                      exclude_regex = NULL) {
  
  base <- basename(files)
  
  # Respect exact-name priority order.
  if (length(exact_names) > 0) {
    for (nm in exact_names) {
      exact_hit <- files[tolower(base) == tolower(nm)]
      if (length(exact_hit) > 0) {
        return(exact_hit[1])
      }
    }
  }
  
  keep <- rep(TRUE, length(files))
  
  if (!is.null(include_regex)) {
    keep <- keep & grepl(include_regex, base, ignore.case = TRUE)
  }
  
  if (!is.null(exclude_regex)) {
    keep <- keep & !grepl(exclude_regex, base, ignore.case = TRUE)
  }
  
  hit <- files[keep]
  
  if (length(hit) == 0) {
    return(NA_character_)
  }
  
  hit[1]
}

detect_raw_files_for_study <- function(folder, subtype_name) {
  files <- list_raw_files(folder)
  
  message("\nFiles found for ", subtype_name, ":")
  print(basename(files))
  
  clinical_file <- find_file(
    files,
    exact_names = c(
      "data_clinical_sample.txt",
      "data_clinical_sample.txt.gz"
    ),
    include_regex = "clinical.*sample"
  )
  
  mrna_file <- find_file(
    files,
    exact_names = c(
      "data_mrna_seq_v2_rsem.txt",
      "data_mrna_seq_v2_rsem.txt.gz",
      "data_RNA_Seq_v2_expression_median.txt",
      "data_RNA_Seq_v2_expression_median.txt.gz"
    ),
    include_regex = "mrna|rna_seq|rsem|expression",
    exclude_regex = "zscore|zscores|z-score|normal|clinical|mutation|mutations|cna|seg|maf|meta"
  )
  
  if (is.na(mrna_file)) {
    mrna_file <- find_file(
      files,
      include_regex = "mrna|rna_seq|rsem|expression",
      exclude_regex = "clinical|mutation|mutations|cna|seg|maf|meta"
    )
  }
  
  # IMPORTANT:
  # Force full HM450 methylation only.
  # Do NOT allow data_methylation_hm27_hm450_merged.txt,
  # because that creates the 670-sample / 22,601-probe workflow.
  methylation_file <- find_file(
    files,
    exact_names = c(
      "data_methylation_hm450.txt",
      "data_methylation_hm450.txt.gz"
    ),
    include_regex = "^data_methylation_hm450",
    exclude_regex = "merged|hm27|clinical|sample|meta"
  )
  
  rppa_file <- find_file(
    files,
    exact_names = c(
      "data_rppa.txt",
      "data_RPPA.txt",
      "data_rppa.txt.gz",
      "data_RPPA.txt.gz"
    ),
    include_regex = "rppa|protein",
    exclude_regex = "zscore|zscores|z-score|clinical|sample|meta"
  )
  
  if (is.na(rppa_file)) {
    rppa_file <- find_file(
      files,
      include_regex = "rppa|protein",
      exclude_regex = "clinical|sample|meta"
    )
  }
  
  tibble::tibble(
    subtype = subtype_name,
    layer = c("clinical_sample", "mRNA", "DNA_methylation", "RPPA"),
    file = c(clinical_file, mrna_file, methylation_file, rppa_file)
  )
}

detect_sample_column <- function(df, path_for_error = "") {
  possible_sample_cols <- c(
    "SAMPLE_ID",
    "Sample_ID",
    "sample_id",
    "Sample ID",
    "sample",
    "Sample",
    "PATIENT_ID",
    "Patient_ID",
    "patient_id"
  )
  
  sample_col <- intersect(possible_sample_cols, colnames(df))
  
  if (length(sample_col) == 0) {
    sample_col <- grep(
      "sample",
      colnames(df),
      ignore.case = TRUE,
      value = TRUE
    )
  }
  
  if (length(sample_col) == 0) {
    stop(
      "Could not detect sample ID column in clinical file:\n",
      path_for_error,
      "\nColumns found:\n",
      paste(colnames(df), collapse = ", ")
    )
  }
  
  sample_col[1]
}

derive_tcga_metadata <- function(sample_ids) {
  sample_parts <- strsplit(sample_ids, "-")
  
  is_tcga_barcode <- vapply(
    sample_parts,
    function(x) length(x) >= 4 && toupper(x[1]) == "TCGA",
    logical(1)
  )
  
  out <- tibble::tibble(sample_id = sample_ids)
  
  if (all(is_tcga_barcode)) {
    out <- out %>%
      dplyr::mutate(
        TCGA_TSS = vapply(sample_parts, function(x) x[2], character(1)),
        TCGA_participant = vapply(sample_parts, function(x) x[3], character(1)),
        TCGA_sample_type_code = substr(
          vapply(sample_parts, function(x) x[4], character(1)),
          1,
          2
        )
      )
  }
  
  out
}

read_clinical_metadata <- function(path, subtype_name) {
  message("\nReading clinical sample metadata for ", subtype_name, ":")
  message(path)
  
  clinical <- read_cbio_table(path)
  
  sample_col <- detect_sample_column(
    clinical,
    path_for_error = path
  )
  
  clinical <- clinical %>%
    dplyr::rename(sample_id = dplyr::all_of(sample_col)) %>%
    dplyr::mutate(
      sample_id = standardize_sample_id(sample_id),
      subtype = subtype_name,
      clinical_source = subtype_name
    ) %>%
    dplyr::filter(!is.na(sample_id), sample_id != "") %>%
    dplyr::distinct(sample_id, .keep_all = TRUE)
  
  barcode_metadata <- derive_tcga_metadata(clinical$sample_id)
  
  clinical <- clinical %>%
    dplyr::select(
      -dplyr::any_of(
        c("TCGA_TSS", "TCGA_participant", "TCGA_sample_type_code")
      )
    ) %>%
    dplyr::left_join(barcode_metadata, by = "sample_id")
  
  clinical
}

read_omics_matrix <- function(path, layer_name, reference_samples) {
  message("\nReading ", layer_name, " file:")
  message(path)
  
  df <- read_cbio_table(path)
  
  message(layer_name, " raw table dimensions:")
  print(dim(df))
  
  original_colnames <- colnames(df)
  standardized_colnames <- standardize_sample_id(original_colnames)
  
  sample_cols <- original_colnames[
    standardized_colnames %in% reference_samples
  ]
  
  if (length(sample_cols) == 0) {
    tcga_like <- grepl(
      "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}",
      standardized_colnames,
      ignore.case = TRUE
    )
    
    sample_cols <- original_colnames[tcga_like]
  }
  
  if (length(sample_cols) == 0) {
    stop(
      layer_name,
      " has no detectable sample columns.\n",
      "Check whether sample IDs in the omics file match clinical sample IDs.\n",
      "File: ",
      path
    )
  }
  
  non_sample_cols <- setdiff(original_colnames, sample_cols)
  
  possible_feature_cols <- c(
    "Hugo_Symbol",
    "Hugo Symbol",
    "GENE_SYMBOL",
    "Gene",
    "gene",
    "Composite.Element.REF",
    "Composite Element REF",
    "ID",
    "Name",
    "NAME"
  )
  
  feature_col <- intersect(possible_feature_cols, non_sample_cols)
  
  if (length(feature_col) == 0) {
    feature_col <- non_sample_cols[1]
  } else {
    feature_col <- feature_col[1]
  }
  
  feature_ids <- as.character(df[[feature_col]])
  
  missing_feature_id <- is.na(feature_ids) | feature_ids == ""
  
  if (any(missing_feature_id)) {
    feature_ids[missing_feature_id] <- paste0(
      "feature_",
      which(missing_feature_id)
    )
  }
  
  feature_ids <- make.unique(feature_ids)
  
  values <- df[, sample_cols, drop = FALSE]
  
  colnames(values) <- standardize_sample_id(colnames(values))
  
  if (any(duplicated(colnames(values)))) {
    warning(
      layer_name,
      ": duplicated sample columns after standardization. Keeping first occurrence."
    )
    
    values <- values[, !duplicated(colnames(values)), drop = FALSE]
  }
  
  values[] <- lapply(values, function(x) {
    suppressWarnings(as.numeric(x))
  })
  
  mat <- as.matrix(values)
  rownames(mat) <- feature_ids
  
  # Convert features × samples to samples × features
  mat <- t(mat)
  
  overlap <- intersect(rownames(mat), reference_samples)
  
  message(layer_name, " sample overlap with clinical labels: ", length(overlap))
  
  if (length(overlap) == 0) {
    stop(
      layer_name,
      " has zero overlap with clinical labels after sample ID standardization.\n",
      "Check sample ID format in clinical and omics files."
    )
  }
  
  mat <- mat[overlap, , drop = FALSE]
  
  message(layer_name, " final matrix:")
  message(nrow(mat), " samples x ", ncol(mat), " features")
  
  mat
}

combine_by_common_features <- function(mat1, mat2, layer_name) {
  common_features <- intersect(colnames(mat1), colnames(mat2))
  
  message("\n", layer_name, " common features between LUAD and LUSC:")
  message(length(common_features))
  
  if (length(common_features) == 0) {
    stop(layer_name, " has zero common features between LUAD and LUSC.")
  }
  
  mat1 <- mat1[, common_features, drop = FALSE]
  mat2 <- mat2[, common_features, drop = FALSE]
  
  combined <- rbind(mat1, mat2)
  
  if (any(duplicated(rownames(combined)))) {
    warning(
      layer_name,
      ": duplicated sample IDs after combining LUAD and LUSC. Keeping first occurrence."
    )
    
    combined <- combined[!duplicated(rownames(combined)), , drop = FALSE]
  }
  
  combined
}

summarize_missing <- function(mat, layer_name) {
  feature_missing <- colMeans(is.na(mat)) * 100
  
  tibble::tibble(
    layer = layer_name,
    samples = nrow(mat),
    features = ncol(mat),
    missing_values = sum(is.na(mat)),
    missing_percent = round(mean(is.na(mat)) * 100, 4),
    features_with_missing = sum(feature_missing > 0),
    max_feature_missing_percent = round(max(feature_missing, na.rm = TRUE), 4)
  )
}

make_duplicate_check <- function(labels, metadata, mrna, methylation, rppa) {
  tibble::tibble(
    object = c(
      "labels sample_id",
      "metadata sample_id",
      "mRNA sample rows",
      "DNA methylation sample rows",
      "RPPA sample rows",
      "mRNA features",
      "DNA methylation features",
      "RPPA features"
    ),
    duplicate_position = c(
      anyDuplicated(labels$sample_id),
      anyDuplicated(metadata$sample_id),
      anyDuplicated(rownames(mrna)),
      anyDuplicated(rownames(methylation)),
      anyDuplicated(rownames(rppa)),
      anyDuplicated(colnames(mrna)),
      anyDuplicated(colnames(methylation)),
      anyDuplicated(colnames(rppa))
    )
  )
}

# -----------------------------
# 3. Detect raw files
# -----------------------------

luad_detected <- detect_raw_files_for_study(
  luad_raw_dir,
  "LUAD"
)

lusc_detected <- detect_raw_files_for_study(
  lusc_raw_dir,
  "LUSC"
)

detected_files <- dplyr::bind_rows(
  luad_detected,
  lusc_detected
)

message("\nDetected raw files:")
print(detected_files)

readr::write_tsv(
  detected_files,
  file.path(tables_dir, "detected_raw_files.tsv")
)

if (any(is.na(detected_files$file))) {
  stop(
    "Some raw files could not be detected automatically.\n",
    "Check results/tables/detected_raw_files.tsv.\n",
    "Then adjust the filename patterns in 01_load_data.R."
  )
}

# Safety check: force HM450 file.
selected_methylation_files <- detected_files$file[
  detected_files$layer == "DNA_methylation"
]

if (!all(grepl("data_methylation_hm450\\.txt(\\.gz)?$", selected_methylation_files, ignore.case = TRUE))) {
  stop(
    "Wrong methylation file selected.\n",
    "Expected data_methylation_hm450.txt for both LUAD and LUSC.\n",
    "Selected files:\n",
    paste(selected_methylation_files, collapse = "\n")
  )
}

# -----------------------------
# 4. Read clinical metadata and create labels
# -----------------------------

luad_clinical_file <- luad_detected$file[
  luad_detected$layer == "clinical_sample"
]

lusc_clinical_file <- lusc_detected$file[
  lusc_detected$layer == "clinical_sample"
]

luad_metadata <- read_clinical_metadata(
  luad_clinical_file,
  "LUAD"
)

lusc_metadata <- read_clinical_metadata(
  lusc_clinical_file,
  "LUSC"
)

metadata <- dplyr::bind_rows(
  luad_metadata,
  lusc_metadata
) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)

labels <- metadata %>%
  dplyr::transmute(
    sample_id = sample_id,
    subtype = subtype
  ) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)

labels$subtype <- factor(labels$subtype, levels = c("LUAD", "LUSC"))

message("\nLabels from raw clinical metadata:")
print(dplyr::count(labels, subtype))

# -----------------------------
# 5. Read raw omics matrices
# -----------------------------

luad_reference_samples <- luad_metadata$sample_id
lusc_reference_samples <- lusc_metadata$sample_id

luad_mrna <- read_omics_matrix(
  luad_detected$file[luad_detected$layer == "mRNA"],
  "LUAD mRNA",
  luad_reference_samples
)

lusc_mrna <- read_omics_matrix(
  lusc_detected$file[lusc_detected$layer == "mRNA"],
  "LUSC mRNA",
  lusc_reference_samples
)

luad_methylation <- read_omics_matrix(
  luad_detected$file[luad_detected$layer == "DNA_methylation"],
  "LUAD DNA methylation",
  luad_reference_samples
)

lusc_methylation <- read_omics_matrix(
  lusc_detected$file[lusc_detected$layer == "DNA_methylation"],
  "LUSC DNA methylation",
  lusc_reference_samples
)

luad_rppa <- read_omics_matrix(
  luad_detected$file[luad_detected$layer == "RPPA"],
  "LUAD RPPA",
  luad_reference_samples
)

lusc_rppa <- read_omics_matrix(
  lusc_detected$file[lusc_detected$layer == "RPPA"],
  "LUSC RPPA",
  lusc_reference_samples
)

# -----------------------------
# 6. Combine LUAD and LUSC by common features
# -----------------------------

mrna <- combine_by_common_features(
  luad_mrna,
  lusc_mrna,
  "mRNA"
)

methylation <- combine_by_common_features(
  luad_methylation,
  lusc_methylation,
  "DNA methylation"
)

rppa <- combine_by_common_features(
  luad_rppa,
  lusc_rppa,
  "RPPA"
)

# -----------------------------
# 7. Match common samples across labels, metadata and omics
# -----------------------------

common_samples <- Reduce(
  intersect,
  list(
    labels$sample_id,
    metadata$sample_id,
    rownames(mrna),
    rownames(methylation),
    rownames(rppa)
  )
)

message("\nCommon samples across labels, metadata, mRNA, methylation and RPPA:")
message(length(common_samples))

if (length(common_samples) == 0) {
  stop(
    "No common samples found across labels, metadata and all omics layers.\n",
    "This usually means sample ID formats still do not match."
  )
}

labels <- labels %>%
  dplyr::filter(sample_id %in% common_samples) %>%
  dplyr::arrange(subtype, sample_id)

common_samples <- labels$sample_id

metadata <- metadata %>%
  dplyr::filter(sample_id %in% common_samples) %>%
  dplyr::arrange(match(sample_id, common_samples))

mrna <- mrna[common_samples, , drop = FALSE]
methylation <- methylation[common_samples, , drop = FALSE]
rppa <- rppa[common_samples, , drop = FALSE]

stopifnot(identical(labels$sample_id, metadata$sample_id))
stopifnot(identical(labels$sample_id, rownames(mrna)))
stopifnot(identical(labels$sample_id, rownames(methylation)))
stopifnot(identical(labels$sample_id, rownames(rppa)))

message("\nLabels, metadata and all omics layers now have identical sample order.")

# -----------------------------
# 8. Summaries
# -----------------------------

dimension_summary <- tibble::tibble(
  layer = c("mRNA", "DNA methylation", "RPPA"),
  samples = c(
    nrow(mrna),
    nrow(methylation),
    nrow(rppa)
  ),
  features = c(
    ncol(mrna),
    ncol(methylation),
    ncol(rppa)
  )
)

class_balance <- labels %>%
  dplyr::count(subtype, name = "n") %>%
  dplyr::mutate(percent = round(n / sum(n) * 100, 2))

missing_summary <- dplyr::bind_rows(
  summarize_missing(mrna, "mRNA"),
  summarize_missing(methylation, "DNA methylation"),
  summarize_missing(rppa, "RPPA")
)

duplicate_check <- make_duplicate_check(
  labels = labels,
  metadata = metadata,
  mrna = mrna,
  methylation = methylation,
  rppa = rppa
)

message("\nDimension summary:")
print(dimension_summary)

message("\nClass balance:")
print(class_balance)

message("\nMissing-value summary:")
print(missing_summary)

message("\nDuplicate check:")
print(duplicate_check)

if (any(duplicate_check$duplicate_position != 0)) {
  stop(
    "Duplicate IDs/features remain after loading. Check duplicate_check output."
  )
}

# Hard safety check for intended 567-sample HM450 dataset.
if (!all(dimension_summary$samples == 567)) {
  stop(
    "Unexpected matched sample count.\n",
    "Expected 567 matched samples for the HM450 workflow.\n",
    "Observed:\n",
    paste(capture.output(print(dimension_summary)), collapse = "\n")
  )
}

if (ncol(methylation) != 396065) {
  stop(
    "Unexpected methylation feature count.\n",
    "Expected 396,065 HM450 probes.\n",
    "Observed: ",
    ncol(methylation)
  )
}

expected_class_balance <- tibble::tibble(
  subtype = factor(c("LUAD", "LUSC"), levels = c("LUAD", "LUSC")),
  expected_n = c(316, 251)
)

observed_class_balance <- class_balance %>%
  dplyr::select(subtype, n)

class_check <- observed_class_balance %>%
  dplyr::left_join(expected_class_balance, by = "subtype")

if (!all(class_check$n == class_check$expected_n)) {
  stop(
    "Unexpected class balance.\n",
    "Expected LUAD = 316 and LUSC = 251.\n",
    "Observed:\n",
    paste(capture.output(print(class_balance)), collapse = "\n")
  )
}

readr::write_tsv(
  dimension_summary,
  file.path(tables_dir, "dimension_summary.tsv")
)

readr::write_tsv(
  class_balance,
  file.path(tables_dir, "class_balance.tsv")
)

readr::write_tsv(
  missing_summary,
  file.path(tables_dir, "missing_value_summary.tsv")
)

readr::write_tsv(
  duplicate_check,
  file.path(tables_dir, "duplicate_check_load_data.tsv")
)

readr::write_tsv(
  metadata,
  file.path(tables_dir, "matched_clinical_metadata.tsv")
)

# -----------------------------
# 9. Optional matched TSV outputs
# -----------------------------

if (isTRUE(save_matched_tsv)) {
  readr::write_tsv(
    labels,
    file.path(processed_dir, "sample_labels.tsv")
  )
  
  readr::write_tsv(
    metadata,
    file.path(processed_dir, "matched_clinical_metadata.tsv")
  )
  
  write_feature_by_sample_tsv <- function(mat, path) {
    out <- as.data.frame(t(mat), check.names = FALSE)
    out <- tibble::rownames_to_column(out, var = "feature_id")
    readr::write_tsv(out, path)
  }
  
  write_feature_by_sample_tsv(
    mrna,
    file.path(processed_dir, "matched_mrna.tsv")
  )
  
  write_feature_by_sample_tsv(
    methylation,
    file.path(processed_dir, "matched_methylation.tsv")
  )
  
  write_feature_by_sample_tsv(
    rppa,
    file.path(processed_dir, "matched_rppa.tsv")
  )
  
  message("\nOptional matched TSV files saved.")
}

# -----------------------------
# 10. Save official loaded matched object
# -----------------------------

loaded_data <- list(
  labels = labels,
  metadata = metadata,
  mrna = mrna,
  methylation = methylation,
  rppa = rppa,
  dimension_summary = dimension_summary,
  class_balance = class_balance,
  missing_summary = missing_summary,
  duplicate_check = duplicate_check,
  detected_files = detected_files,
  note = paste(
    "This object was generated directly from raw LUAD and LUSC folders.",
    "Samples are matched across labels, metadata, mRNA, DNA methylation, and RPPA.",
    "DNA methylation was forced to data_methylation_hm450.txt.",
    "No imputation, scaling, feature selection, batch correction, or modelling was applied."
  )
)

saveRDS(
  loaded_data,
  loaded_rds_path,
  compress = FALSE
)

message("\nSaved official loaded matched object to:")
message(loaded_rds_path)

message("\n01_load_data.R completed successfully.")