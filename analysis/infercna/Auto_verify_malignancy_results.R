project_dir <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline"
out_dir <- file.path(project_dir, "sn_outs")
setwd(project_dir)

manifest_path <- file.path(out_dir, "sample_manifest.csv")
if (!file.exists(manifest_path)) {
  stop("Missing sample manifest: ", manifest_path)
}

read_status <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  trimws(readLines(path, warn = FALSE, n = 1))
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE)
}

value_or_zero <- function(x, name) {
  if (!name %in% names(x)) {
    return(0L)
  }
  as.integer(x[[name]])
}

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)

status_rows <- do.call(rbind, lapply(manifest$sample, function(sample) {
  epi_f_path <- file.path(out_dir, "by_samples", sample, paste0(sample, "_epi_f.rds"))
  sample_summary_path <- file.path(out_dir, "by_samples", sample, paste0(sample, "_malignancy_summary.csv"))
  cluster_summary_path <- file.path(out_dir, "by_samples", sample, paste0(sample, "_malignancy_cluster_summary.csv"))

  data.frame(
    sample = sample,
    expr_filter_status = read_status(file.path(out_dir, "by_samples", sample, "expr_filter_status.txt")),
    infercna_status = read_status(file.path(out_dir, "by_samples", sample, "infercna_status.txt")),
    malignancy_status = read_status(file.path(out_dir, "by_samples", sample, "malignancy_status.txt")),
    epi_f_present = file.exists(epi_f_path),
    sample_summary_present = file.exists(sample_summary_path),
    cluster_summary_present = file.exists(cluster_summary_path),
    stringsAsFactors = FALSE
  )
}))

ok_samples <- status_rows$sample[
  status_rows$expr_filter_status == "ok" &
    status_rows$infercna_status == "ok" &
    status_rows$malignancy_status == "ok" &
    status_rows$epi_f_present &
    status_rows$sample_summary_present &
    status_rows$cluster_summary_present
]

if (length(ok_samples) == 0) {
  stop("No completed malignancy outputs found to verify.")
}

required_meta_cols <- c(
  "orig.ident",
  "reference_batch",
  "classification",
  "malignant_clus",
  "malignancy"
)

integrity_rows <- do.call(rbind, lapply(ok_samples, function(sample) {
  epi_f_path <- file.path(out_dir, "by_samples", sample, paste0(sample, "_epi_f.rds"))
  sample_summary_path <- file.path(out_dir, "by_samples", sample, paste0(sample, "_malignancy_summary.csv"))
  cluster_summary_path <- file.path(out_dir, "by_samples", sample, paste0(sample, "_malignancy_cluster_summary.csv"))

  epi <- readRDS(epi_f_path)
  meta <- methods::slot(epi, "meta.data")
  sample_summary <- safe_read_csv(sample_summary_path)
  cluster_summary <- safe_read_csv(cluster_summary_path)

  missing_meta <- setdiff(required_meta_cols, colnames(meta))

  class_counts <- table(meta$classification)
  malignancy_counts <- table(meta$malignancy)

  class_total <- sum(class_counts, na.rm = TRUE)
  malignancy_total <- sum(malignancy_counts, na.rm = TRUE)
  cluster_total <- if (!is.null(cluster_summary) && nrow(cluster_summary) > 0) {
    sum(cluster_summary$n_cells, na.rm = TRUE)
  } else {
    NA_integer_
  }

  summary_row <- if (!is.null(sample_summary) && nrow(sample_summary) > 0) {
    sample_summary[1, , drop = FALSE]
  } else {
    data.frame()
  }

  summary_epithelial_cells <- if (nrow(summary_row) > 0 && "epithelial_cells" %in% colnames(summary_row)) {
    as.integer(summary_row$epithelial_cells[1])
  } else {
    NA_integer_
  }

  summary_class_total <- if (nrow(summary_row) > 0) {
    value_or_zero(summary_row, "cna_malignant_cells") +
      value_or_zero(summary_row, "cna_unresolved_cells") +
      value_or_zero(summary_row, "cna_non_malignant_cells")
  } else {
    NA_integer_
  }

  summary_malignancy_total <- if (nrow(summary_row) > 0) {
    value_or_zero(summary_row, "malignant_level_1_cells") +
      value_or_zero(summary_row, "malignant_level_2_cells") +
      value_or_zero(summary_row, "malignant_unresolved_cells") +
      value_or_zero(summary_row, "non_malignant_level_1_cells") +
      value_or_zero(summary_row, "non_malignant_level_2_cells") +
      value_or_zero(summary_row, "non_malignant_unresolved_cells") +
      value_or_zero(summary_row, "unresolved_cells")
  } else {
    NA_integer_
  }

  check_map <- c(
    meta_columns_present = length(missing_meta) == 0,
    epi_matches_sample_summary = !is.na(summary_epithelial_cells) && nrow(meta) == summary_epithelial_cells,
    classification_total_matches_object = class_total == nrow(meta),
    classification_total_matches_summary = !is.na(summary_class_total) && class_total == summary_class_total,
    malignancy_total_matches_object = malignancy_total == nrow(meta),
    malignancy_total_matches_summary = !is.na(summary_malignancy_total) && malignancy_total == summary_malignancy_total,
    cluster_total_matches_object = !is.na(cluster_total) && cluster_total == nrow(meta),
    single_orig_ident = length(unique(meta$orig.ident)) == 1,
    single_reference_batch = length(unique(meta$reference_batch)) == 1
  )

  failed_checks <- names(check_map)[!check_map]
  if (length(missing_meta) > 0) {
    failed_checks <- c(failed_checks, paste0("missing_meta:", missing_meta))
  }

  data.frame(
    sample = sample,
    epithelial_cells = nrow(meta),
    malignant_total_cells = sum(meta$malignancy %in% c("malignant_level_1", "malignant_level_2"), na.rm = TRUE),
    unresolved_cells = sum(meta$malignancy == "unresolved", na.rm = TRUE),
    classification_total = class_total,
    malignancy_total = malignancy_total,
    cluster_total = cluster_total,
    summary_epithelial_cells = summary_epithelial_cells,
    reference_batch = paste(unique(meta$reference_batch), collapse = ";"),
    passed = length(failed_checks) == 0,
    failed_checks = paste(failed_checks, collapse = " | "),
    stringsAsFactors = FALSE
  )
}))

summary_row <- data.frame(
  total_manifest_samples = nrow(status_rows),
  expr_filter_ok_samples = sum(status_rows$expr_filter_status == "ok", na.rm = TRUE),
  infercna_ok_samples = sum(status_rows$expr_filter_status == "ok" & status_rows$infercna_status == "ok", na.rm = TRUE),
  malignancy_ok_samples = length(ok_samples),
  verified_samples = nrow(integrity_rows),
  passed_samples = sum(integrity_rows$passed, na.rm = TRUE),
  failed_samples = sum(!integrity_rows$passed, na.rm = TRUE),
  epithelial_cells_verified = sum(integrity_rows$epithelial_cells, na.rm = TRUE),
  malignant_cells_verified = sum(integrity_rows$malignant_total_cells, na.rm = TRUE),
  unresolved_cells_verified = sum(integrity_rows$unresolved_cells, na.rm = TRUE),
  stringsAsFactors = FALSE
)

write.csv(
  integrity_rows,
  file.path(out_dir, "Auto_malignancy_integrity_by_sample.csv"),
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  summary_row,
  file.path(out_dir, "Auto_malignancy_integrity_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)

if (any(!integrity_rows$passed)) {
  bad_samples <- integrity_rows$sample[!integrity_rows$passed]
  stop(
    "Malignancy integrity check failed for: ",
    paste(bad_samples, collapse = ", ")
  )
}

message(
  "Verified malignancy outputs for ",
  nrow(integrity_rows),
  " samples; malignant cells available for downstream merge: ",
  sum(integrity_rows$malignant_total_cells, na.rm = TRUE)
)
