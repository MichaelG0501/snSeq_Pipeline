suppressPackageStartupMessages({
  library("dplyr")
  library("Seurat")
})

set.seed(1)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
message(sample)

status_path <- file.path("by_samples", sample, "malignancy_status.txt")
sample_summary_path <- file.path("by_samples", sample, paste0(sample, "_malignancy_summary.csv"))
cluster_summary_path <- file.path("by_samples", sample, paste0(sample, "_malignancy_cluster_summary.csv"))

write_status <- function(status) {
  writeLines(status, status_path)
}

write_sample_summary <- function(df) {
  write.csv(df, sample_summary_path, row.names = FALSE, quote = TRUE)
}

write_cluster_summary <- function(df) {
  write.csv(df, cluster_summary_path, row.names = FALSE, quote = TRUE)
}

count_value <- function(x, value) {
  sum(x == value, na.rm = TRUE)
}

score_signature <- function(obj, genes) {
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(obj)))
  }
  colMeans(as.matrix(obj@assays$RNA$data[genes, , drop = FALSE]), na.rm = TRUE)
}

epi_path <- file.path("by_samples", sample, paste0(sample, "_epi.rds"))
write_status("running")
if (!file.exists(epi_path)) {
  write_sample_summary(data.frame(
    sample = sample,
    status = "missing_epi",
    stringsAsFactors = FALSE
  ))
  write_status("missing_epi")
  stop("Missing inferCNA epithelial object.")
}

epi <- readRDS(epi_path)
required_meta_cols <- c("classification", "malignant_clus", "seurat_clusters")
missing_meta_cols <- setdiff(required_meta_cols, colnames(epi@meta.data))
if (length(missing_meta_cols) > 0) {
  stop(
    "Missing metadata columns in ", sample, "_epi.rds: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

cell_cycle_genes <- read.csv(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)[, 1:3]

cc_genes <- cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1]
cc_genes <- intersect(cc_genes, rownames(epi))
if (length(cc_genes) > 0) {
  gene_means <- rowMeans(epi@assays$RNA$data[cc_genes, , drop = FALSE], na.rm = TRUE)
  cc_genes <- names(sort(gene_means, decreasing = TRUE))[seq_len(min(10, length(gene_means)))]
}

cc_score <- score_signature(epi, cc_genes)
epi$cc_score <- cc_score
cc_thr <- 1
cc_status <- rep("cc_unresolved", ncol(epi))
names(cc_status) <- colnames(epi)
cc_status[!is.na(cc_score) & cc_score >= cc_thr] <- "cc_malignant"
epi$cc_status <- cc_status

sig_file <- "cancer_signatures.txt"
if (file.exists(sig_file) && file.info(sig_file)$size > 0) {
  cs_genes <- read.table(sig_file, header = FALSE, stringsAsFactors = FALSE)$V1
  cs_genes <- intersect(cs_genes, rownames(epi))
  if (length(cs_genes) > 0) {
    gene_means <- rowMeans(epi@assays$RNA$data[cs_genes, , drop = FALSE], na.rm = TRUE)
    cs_genes <- names(sort(gene_means, decreasing = TRUE))[seq_len(min(50, length(gene_means)))]
  }
} else {
  cs_genes <- character(0)
}

cs_score <- score_signature(epi, cs_genes)
epi$cs_score <- cs_score
cs_thr <- 3
cs_status <- rep("cs_unresolved", ncol(epi))
names(cs_status) <- colnames(epi)
cs_status[!is.na(cs_score) & cs_score >= cs_thr] <- "cs_malignant"
epi$cs_status <- cs_status

epi$malignancy <- rep("unresolved", ncol(epi))
cluster_labels <- tapply(epi$malignant_clus, epi$seurat_clusters, function(x) {
  if (all(x == "malignant_clus")) {
    "malignant_clus"
  } else if (all(x == "non_malignant_clus")) {
    "non_malignant_clus"
  } else if (all(x == "unresolved_clus")) {
    "unresolved_clus"
  } else {
    NA_character_
  }
})

for (cl in names(cluster_labels)) {
  if (cluster_labels[cl] == "non_malignant_clus") {
    epi$malignancy[epi$seurat_clusters == cl & epi$classification == "cna_non_malignant"] <- "non_malignant_level_1"
    epi$malignancy[epi$seurat_clusters == cl & epi$classification == "cna_unresolved"] <- "non_malignant_level_2"
    epi$malignancy[epi$seurat_clusters == cl & epi$classification == "cna_malignant"] <- "non_malignant_unresolved"
  } else {
    malignant <- epi$seurat_clusters == cl & (
      epi$classification == "cna_malignant" |
        epi$cc_status == "cc_malignant" |
        epi$cs_status == "cs_malignant"
    )

    if (sum(malignant) / sum(epi$seurat_clusters == cl) >= 0.5) {
      epi$malignancy[epi$seurat_clusters == cl & epi$classification == "cna_malignant"] <- "malignant_level_1"
      epi$malignancy[epi$seurat_clusters == cl & epi$classification == "cna_unresolved"] <- "malignant_level_2"
      epi$malignancy[epi$seurat_clusters == cl & epi$classification == "cna_non_malignant"] <- "malignant_unresolved"
    } else {
      epi$malignancy[epi$seurat_clusters == cl] <- "unresolved"
    }
  }
}

cluster_summary <- lapply(names(cluster_labels), function(cl) {
  cl_cells <- epi$seurat_clusters == cl
  total <- sum(cl_cells, na.rm = TRUE)
  data.frame(
    sample = sample,
    seurat_clusters = as.character(cl),
    malignant_clus = as.character(cluster_labels[[cl]]),
    n_cells = total,
    cna_malignant = count_value(epi$classification[cl_cells], "cna_malignant"),
    cna_unresolved = count_value(epi$classification[cl_cells], "cna_unresolved"),
    cna_non_malignant = count_value(epi$classification[cl_cells], "cna_non_malignant"),
    cs_malignant = count_value(epi$cs_status[cl_cells], "cs_malignant"),
    cc_malignant = count_value(epi$cc_status[cl_cells], "cc_malignant"),
    malignant_level_1 = count_value(epi$malignancy[cl_cells], "malignant_level_1"),
    malignant_level_2 = count_value(epi$malignancy[cl_cells], "malignant_level_2"),
    malignant_unresolved = count_value(epi$malignancy[cl_cells], "malignant_unresolved"),
    non_malignant_level_1 = count_value(epi$malignancy[cl_cells], "non_malignant_level_1"),
    non_malignant_level_2 = count_value(epi$malignancy[cl_cells], "non_malignant_level_2"),
    non_malignant_unresolved = count_value(epi$malignancy[cl_cells], "non_malignant_unresolved"),
    unresolved = count_value(epi$malignancy[cl_cells], "unresolved"),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()
write_cluster_summary(cluster_summary)

sample_summary <- data.frame(
  sample = sample,
  status = "ok",
  reference_batch = unique(epi$reference_batch)[1],
  technology = unique(epi$technology)[1],
  orig_ident = unique(epi$orig.ident)[1],
  epithelial_cells = ncol(epi),
  cc_threshold = cc_thr,
  cs_threshold = cs_thr,
  cna_malignant_cells = count_value(epi$classification, "cna_malignant"),
  cna_unresolved_cells = count_value(epi$classification, "cna_unresolved"),
  cna_non_malignant_cells = count_value(epi$classification, "cna_non_malignant"),
  malignant_clus_cells = count_value(epi$malignant_clus, "malignant_clus"),
  non_malignant_clus_cells = count_value(epi$malignant_clus, "non_malignant_clus"),
  unresolved_clus_cells = count_value(epi$malignant_clus, "unresolved_clus"),
  cs_malignant_cells = count_value(epi$cs_status, "cs_malignant"),
  cc_malignant_cells = count_value(epi$cc_status, "cc_malignant"),
  malignant_level_1_cells = count_value(epi$malignancy, "malignant_level_1"),
  malignant_level_2_cells = count_value(epi$malignancy, "malignant_level_2"),
  malignant_unresolved_cells = count_value(epi$malignancy, "malignant_unresolved"),
  non_malignant_level_1_cells = count_value(epi$malignancy, "non_malignant_level_1"),
  non_malignant_level_2_cells = count_value(epi$malignancy, "non_malignant_level_2"),
  non_malignant_unresolved_cells = count_value(epi$malignancy, "non_malignant_unresolved"),
  unresolved_cells = count_value(epi$malignancy, "unresolved"),
  malignant_clus_n = sum(cluster_summary$malignant_clus == "malignant_clus", na.rm = TRUE),
  non_malignant_clus_n = sum(cluster_summary$malignant_clus == "non_malignant_clus", na.rm = TRUE),
  unresolved_clus_n = sum(cluster_summary$malignant_clus == "unresolved_clus", na.rm = TRUE),
  stringsAsFactors = FALSE
)

if (sample_summary$epithelial_cells > 0) {
  sample_summary$pct_malignant_level_1 <- sample_summary$malignant_level_1_cells / sample_summary$epithelial_cells
  sample_summary$pct_malignant_level_2 <- sample_summary$malignant_level_2_cells / sample_summary$epithelial_cells
  sample_summary$pct_malignant_total <- (sample_summary$malignant_level_1_cells + sample_summary$malignant_level_2_cells) / sample_summary$epithelial_cells
  sample_summary$pct_non_malignant_total <- (sample_summary$non_malignant_level_1_cells + sample_summary$non_malignant_level_2_cells) / sample_summary$epithelial_cells
  sample_summary$pct_unresolved <- sample_summary$unresolved_cells / sample_summary$epithelial_cells
} else {
  sample_summary$pct_malignant_level_1 <- NA_real_
  sample_summary$pct_malignant_level_2 <- NA_real_
  sample_summary$pct_malignant_total <- NA_real_
  sample_summary$pct_non_malignant_total <- NA_real_
  sample_summary$pct_unresolved <- NA_real_
}

write_sample_summary(sample_summary)
saveRDS(epi, file.path("by_samples", sample, paste0(sample, "_epi_f.rds")))
write_status("ok")
