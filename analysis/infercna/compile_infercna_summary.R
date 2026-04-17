suppressPackageStartupMessages({
  library("dplyr")
})

pipeline_root <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline"
sn_root <- file.path(pipeline_root, "sn_outs")
setwd(pipeline_root)

manifest_path <- file.path(sn_root, "sample_manifest.csv")
filtered_summary_path <- file.path(sn_root, "filtered_sample_summary.csv")

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
filtered_summary <- read.csv(filtered_summary_path, stringsAsFactors = FALSE)

read_status <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  trimws(readLines(path, warn = FALSE, n = 1))
}

read_optional_csv <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE)
}

status_rows <- lapply(manifest$sample, function(sample) {
  expr_status <- read_status(file.path(sn_root, "by_samples", sample, "expr_filter_status.txt"))
  infer_status <- read_status(file.path(sn_root, "by_samples", sample, "infercna_status.txt"))
  sample_summary_path <- file.path(sn_root, "by_samples", sample, paste0(sample, "_infercna_summary.csv"))
  cluster_summary_path <- file.path(sn_root, "by_samples", sample, paste0(sample, "_infercna_cluster_summary.csv"))
  epi_path <- file.path(sn_root, "by_samples", sample, paste0(sample, "_epi.rds"))
  signature_path <- file.path(sn_root, "by_samples", sample, paste0(sample, "_signatures.rds"))
  sample_summary_present <- file.exists(sample_summary_path)
  cluster_summary_present <- file.exists(cluster_summary_path)
  epi_present <- file.exists(epi_path)
  signature_present <- file.exists(signature_path)

  effective_infer_status <- infer_status
  if (identical(expr_status, "ok") && identical(infer_status, "ok") &&
      !(sample_summary_present && cluster_summary_present && epi_present && signature_present)) {
    effective_infer_status <- "ok_missing_outputs"
  }

  data.frame(
    sample = sample,
    expr_filter_status = expr_status,
    infercna_status = infer_status,
    effective_infercna_status = effective_infer_status,
    infercna_summary_present = sample_summary_present,
    infercna_cluster_summary_present = cluster_summary_present,
    infercna_epi_present = epi_present,
    infercna_signature_present = signature_present,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

status_rows <- status_rows %>%
  left_join(
    filtered_summary %>%
      select(sample, filtered_status = status, cells_before, cells_after),
    by = "sample"
  )

write.csv(
  status_rows,
  file.path(sn_root, "infercna_status_by_sample.csv"),
  row.names = FALSE,
  quote = TRUE
)

expr_ok_samples <- status_rows$sample[status_rows$expr_filter_status == "ok"]

sample_summaries <- lapply(expr_ok_samples, function(sample) {
  infer_status <- status_rows$effective_infercna_status[match(sample, status_rows$sample)]
  if (!identical(infer_status, "ok")) {
    return(NULL)
  }
  read_optional_csv(file.path(sn_root, "by_samples", sample, paste0(sample, "_infercna_summary.csv")))
}) %>% bind_rows()

cluster_summaries <- lapply(expr_ok_samples, function(sample) {
  infer_status <- status_rows$effective_infercna_status[match(sample, status_rows$sample)]
  if (!identical(infer_status, "ok")) {
    return(NULL)
  }
  out <- read_optional_csv(file.path(sn_root, "by_samples", sample, paste0(sample, "_infercna_cluster_summary.csv")))
  if (is.null(out)) {
    return(NULL)
  }
  out$seurat_clusters <- as.character(out$seurat_clusters)
  out$malignant_clus <- as.character(out$malignant_clus)
  out
}) %>% bind_rows()

signature_membership <- lapply(expr_ok_samples, function(sample) {
  infer_status <- status_rows$effective_infercna_status[match(sample, status_rows$sample)]
  if (!identical(infer_status, "ok")) {
    return(NULL)
  }
  sig_path <- file.path(sn_root, "by_samples", sample, paste0(sample, "_signatures.rds"))
  if (!file.exists(sig_path)) {
    return(NULL)
  }
  sig <- readRDS(sig_path)
  if (length(sig) == 0) {
    return(NULL)
  }
  data.frame(
    sample = sample,
    gene = unique(sig),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

if (nrow(sample_summaries) > 0) {
  sample_summaries <- sample_summaries %>%
    left_join(
      filtered_summary %>%
        select(sample, expr_filtered_cells = cells_after),
      by = "sample"
    ) %>%
    relocate(expr_filtered_cells, .after = filtered_cells)

  write.csv(
    sample_summaries,
    file.path(sn_root, "infercna_sample_summary.csv"),
    row.names = FALSE,
    quote = TRUE
  )

  overall_summary <- data.frame(
    expr_filter_ok_samples = sum(status_rows$expr_filter_status == "ok", na.rm = TRUE),
    infercna_ok_samples = sum(status_rows$expr_filter_status == "ok" & status_rows$effective_infercna_status == "ok", na.rm = TRUE),
    infercna_ok_missing_outputs = sum(status_rows$expr_filter_status == "ok" & status_rows$effective_infercna_status == "ok_missing_outputs", na.rm = TRUE),
    infercna_no_epi_samples = sum(status_rows$expr_filter_status == "ok" & status_rows$effective_infercna_status == "no_epi", na.rm = TRUE),
    infercna_no_cell_samples = sum(status_rows$expr_filter_status == "ok" & status_rows$effective_infercna_status == "no_cell", na.rm = TRUE),
    infercna_no_ref_samples = sum(status_rows$expr_filter_status == "ok" & status_rows$effective_infercna_status == "no_ref", na.rm = TRUE),
    total_filtered_cells_in_expr_ok_samples = sum(filtered_summary$cells_after[match(expr_ok_samples, filtered_summary$sample)], na.rm = TRUE),
    total_epithelial_cells = sum(sample_summaries$epithelial_cells, na.rm = TRUE),
    cna_malignant_cells = sum(sample_summaries$cna_malignant_cells, na.rm = TRUE),
    cna_unresolved_cells = sum(sample_summaries$cna_unresolved_cells, na.rm = TRUE),
    cna_non_malignant_cells = sum(sample_summaries$cna_non_malignant_cells, na.rm = TRUE),
    malignant_clus_cells = sum(sample_summaries$malignant_clus_cells, na.rm = TRUE),
    non_malignant_clus_cells = sum(sample_summaries$non_malignant_clus_cells, na.rm = TRUE),
    unresolved_clus_cells = sum(sample_summaries$unresolved_clus_cells, na.rm = TRUE),
    malignant_clus_n = sum(sample_summaries$malignant_clus_n, na.rm = TRUE),
    non_malignant_clus_n = sum(sample_summaries$non_malignant_clus_n, na.rm = TRUE),
    unresolved_clus_n = sum(sample_summaries$unresolved_clus_n, na.rm = TRUE),
    samples_with_signature_genes = sum(sample_summaries$signature_gene_count > 0, na.rm = TRUE),
    union_signature_genes = if (nrow(signature_membership) > 0) dplyr::n_distinct(signature_membership$gene) else 0L,
    stringsAsFactors = FALSE
  )

  if (overall_summary$total_epithelial_cells > 0) {
    overall_summary$pct_cna_malignant <- overall_summary$cna_malignant_cells / overall_summary$total_epithelial_cells
    overall_summary$pct_cna_unresolved <- overall_summary$cna_unresolved_cells / overall_summary$total_epithelial_cells
    overall_summary$pct_cna_non_malignant <- overall_summary$cna_non_malignant_cells / overall_summary$total_epithelial_cells
    overall_summary$pct_malignant_clus <- overall_summary$malignant_clus_cells / overall_summary$total_epithelial_cells
    overall_summary$pct_non_malignant_clus <- overall_summary$non_malignant_clus_cells / overall_summary$total_epithelial_cells
    overall_summary$pct_unresolved_clus <- overall_summary$unresolved_clus_cells / overall_summary$total_epithelial_cells
  }

  write.csv(
    overall_summary,
    file.path(sn_root, "infercna_overall_summary.csv"),
    row.names = FALSE,
    quote = TRUE
  )

  classification_overall <- data.frame(
    classification = c("cna_malignant", "cna_unresolved", "cna_non_malignant"),
    cells = c(
      sum(sample_summaries$cna_malignant_cells, na.rm = TRUE),
      sum(sample_summaries$cna_unresolved_cells, na.rm = TRUE),
      sum(sample_summaries$cna_non_malignant_cells, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  classification_overall$pct_of_epithelial <- if (sum(classification_overall$cells) > 0) {
    classification_overall$cells / sum(classification_overall$cells)
  } else {
    NA_real_
  }

  write.csv(
    classification_overall,
    file.path(sn_root, "infercna_classification_overall.csv"),
    row.names = FALSE,
    quote = TRUE
  )

  cluster_label_overall <- data.frame(
    malignant_clus = c("malignant_clus", "non_malignant_clus", "unresolved_clus"),
    cells = c(
      sum(sample_summaries$malignant_clus_cells, na.rm = TRUE),
      sum(sample_summaries$non_malignant_clus_cells, na.rm = TRUE),
      sum(sample_summaries$unresolved_clus_cells, na.rm = TRUE)
    ),
    n_clusters = c(
      sum(sample_summaries$malignant_clus_n, na.rm = TRUE),
      sum(sample_summaries$non_malignant_clus_n, na.rm = TRUE),
      sum(sample_summaries$unresolved_clus_n, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  cluster_label_overall$pct_of_epithelial <- if (sum(cluster_label_overall$cells) > 0) {
    cluster_label_overall$cells / sum(cluster_label_overall$cells)
  } else {
    NA_real_
  }

  write.csv(
    cluster_label_overall,
    file.path(sn_root, "infercna_cluster_label_overall.csv"),
    row.names = FALSE,
    quote = TRUE
  )
} else {
  write.csv(
    data.frame(),
    file.path(sn_root, "infercna_sample_summary.csv"),
    row.names = FALSE,
    quote = TRUE
  )
  write.csv(
    data.frame(),
    file.path(sn_root, "infercna_overall_summary.csv"),
    row.names = FALSE,
    quote = TRUE
  )
  write.csv(
    data.frame(),
    file.path(sn_root, "infercna_classification_overall.csv"),
    row.names = FALSE,
    quote = TRUE
  )
  write.csv(
    data.frame(),
    file.path(sn_root, "infercna_cluster_label_overall.csv"),
    row.names = FALSE,
    quote = TRUE
  )
}

if (nrow(cluster_summaries) > 0) {
  write.csv(
    cluster_summaries,
    file.path(sn_root, "infercna_cluster_summary.csv"),
    row.names = FALSE,
    quote = TRUE
  )
} else {
  write.csv(
    data.frame(),
    file.path(sn_root, "infercna_cluster_summary.csv"),
    row.names = FALSE,
    quote = TRUE
  )
}

if (nrow(signature_membership) > 0) {
  signature_summary <- signature_membership %>%
    group_by(gene) %>%
    summarise(
      n_samples = n_distinct(sample),
      samples = paste(sort(unique(sample)), collapse = "|"),
      .groups = "drop"
    ) %>%
    arrange(desc(n_samples), gene)

  write.csv(
    signature_membership,
    file.path(sn_root, "cancer_signatures_by_sample.csv"),
    row.names = FALSE,
    quote = TRUE
  )
  write.csv(
    signature_summary,
    file.path(sn_root, "cancer_signatures_summary.csv"),
    row.names = FALSE,
    quote = TRUE
  )
  write.table(
    signature_summary$gene,
    file.path(sn_root, "cancer_signatures.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
} else {
  write.csv(
    data.frame(sample = character(0), gene = character(0)),
    file.path(sn_root, "cancer_signatures_by_sample.csv"),
    row.names = FALSE,
    quote = TRUE
  )
  write.csv(
    data.frame(gene = character(0), n_samples = integer(0), samples = character(0)),
    file.path(sn_root, "cancer_signatures_summary.csv"),
    row.names = FALSE,
    quote = TRUE
  )
  write.table(
    character(0),
    file.path(sn_root, "cancer_signatures.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}
