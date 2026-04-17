suppressPackageStartupMessages({
  library("dplyr")
  library("Seurat")
  library("infercna")
})

set.seed(1)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
message(sample)

status_path <- file.path("by_samples", sample, "infercna_status.txt")
signature_path <- file.path("by_samples", sample, paste0(sample, "_signatures.rds"))
sample_summary_path <- file.path("by_samples", sample, paste0(sample, "_infercna_summary.csv"))
cluster_summary_path <- file.path("by_samples", sample, paste0(sample, "_infercna_cluster_summary.csv"))

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

get_layer <- function(obj, layer) {
  layer_data <- tryCatch(obj@assays$RNA@layers[[layer]], error = function(e) NULL)
  if (!is.null(layer_data)) {
    if (length(rownames(layer_data)) == 0) {
      rownames(layer_data) <- rownames(obj)
    }
    if (length(colnames(layer_data)) == 0) {
      colnames(layer_data) <- colnames(obj)
    }
    return(layer_data)
  }

  if (layer %in% c("counts", "data")) {
    layer_data <- GetAssayData(obj, slot = layer)
    if (length(rownames(layer_data)) == 0) {
      rownames(layer_data) <- rownames(obj)
    }
    if (length(colnames(layer_data)) == 0) {
      colnames(layer_data) <- colnames(obj)
    }
    return(layer_data)
  }

  if (layer == "CPM") {
    layer_data <- obj@assays$RNA$CPM
    if (length(rownames(layer_data)) == 0) {
      rownames(layer_data) <- rownames(obj)
    }
    if (length(colnames(layer_data)) == 0) {
      colnames(layer_data) <- colnames(obj)
    }
    return(layer_data)
  }

  stop("Layer not available: ", layer)
}

anno_path <- file.path("by_samples", sample, paste0(sample, "_anno.rds"))
write_status("running")
saveRDS(character(0), signature_path)

if (!file.exists(anno_path)) {
  write_sample_summary(data.frame(
    sample = sample,
    status = "missing_anno",
    reference_batch = NA_character_,
    stringsAsFactors = FALSE
  ))
  write_status("missing_anno")
  stop("Missing annotated object: ", anno_path)
}

tmdata <- readRDS(anno_path)
required_meta_cols <- c(
  "orig.ident", "celltype_update", "coexpression_loose", "marker_expression", "reference_batch"
)
missing_meta_cols <- setdiff(required_meta_cols, colnames(tmdata@meta.data))
if (length(missing_meta_cols) > 0) {
  stop(
    "Missing metadata columns in ", sample, ": ",
    paste(missing_meta_cols, collapse = ", "),
    ". Rerun earlier steps."
  )
}

reference_batch <- unique(tmdata$reference_batch)
reference_batch <- reference_batch[!is.na(reference_batch) & nzchar(reference_batch)]
if (length(reference_batch) != 1) {
  stop("Expected exactly one reference_batch in sample ", sample)
}

prefixed_cells <- paste0(tmdata$orig.ident, "_", colnames(tmdata))
tmdata <- RenameCells(tmdata, new.names = prefixed_cells)

if (sum(tmdata$coexpression_loose == "singlet" & tmdata$marker_expression == "good", na.rm = TRUE) > 0) {
  tmdata <- subset(tmdata, coexpression_loose == "singlet" & marker_expression == "good")
} else {
  write_sample_summary(data.frame(
    sample = sample,
    status = "no_cell",
    reference_batch = reference_batch,
    filtered_cells = 0L,
    epithelial_cells = 0L,
    stringsAsFactors = FALSE
  ))
  write_status("no_cell")
  saveRDS("No cells", file.path("by_samples", sample, "no_cell"))
  stop("No cells")
}

if (sum(tmdata$celltype_update == "epithelial", na.rm = TRUE) > 0) {
  data <- subset(tmdata, celltype_update %in% "epithelial")
} else {
  write_sample_summary(data.frame(
    sample = sample,
    status = "no_epi",
    reference_batch = reference_batch,
    filtered_cells = ncol(tmdata),
    epithelial_cells = 0L,
    stringsAsFactors = FALSE
  ))
  write_status("no_epi")
  saveRDS("No epithelial cells", file.path("by_samples", sample, "no_epi"))
  stop("No epithelial cells")
}

reference_path <- paste0(reference_batch, "_reference.rds")
if (!file.exists(reference_path)) {
  write_sample_summary(data.frame(
    sample = sample,
    status = "no_ref",
    reference_batch = reference_batch,
    filtered_cells = ncol(tmdata),
    epithelial_cells = ncol(data),
    stringsAsFactors = FALSE
  ))
  write_status("no_ref")
  saveRDS("No reference", file.path("by_samples", sample, "no_ref"))
  stop("Missing reference: ", reference_path)
}

reference <- readRDS(reference_path)
required_ref_fields <- c("matrix", "ref", "meta", "ref_ct")
missing_ref_fields <- setdiff(required_ref_fields, names(reference))
if (length(missing_ref_fields) > 0) {
  stop(
    "Reference file is missing fields: ",
    paste(missing_ref_fields, collapse = ", "),
    " in ", reference_path
  )
}

meta <- data@meta.data[, c("orig.ident", "celltype_update"), drop = FALSE]
matrix_data <- as.matrix(get_layer(data, "CPM"))
common_genes <- intersect(rownames(matrix_data), rownames(reference$matrix))
if (length(common_genes) == 0) {
  write_sample_summary(data.frame(
    sample = sample,
    status = "no_ref",
    reference_batch = reference_batch,
    filtered_cells = ncol(tmdata),
    epithelial_cells = ncol(data),
    shared_genes = 0L,
    stringsAsFactors = FALSE
  ))
  write_status("no_ref")
  saveRDS("No shared genes with reference", file.path("by_samples", sample, "no_ref"))
  stop("No shared genes with reference.")
}

matrix_data <- matrix_data[common_genes, , drop = FALSE]
reference_matrix <- as.matrix(reference$matrix[common_genes, , drop = FALSE])
meta <- rbind(meta, reference$meta)
matrix_data <- cbind(matrix_data, reference_matrix)
meta <- meta[colnames(matrix_data), , drop = FALSE]
ref <- reference$ref

outs <- infercna(matrix_data, refCells = ref)
saveRDS(outs, file.path("by_samples", sample, paste0(sample, "_outs.rds")))

sd_k_cor <- 2
sd_k_sig <- 2

coord <- cnaScatterPlot(outs)
meta <- meta %>%
  as.data.frame() %>%
  mutate(
    cna_cor = coord$cna.cor,
    cna_signal = coord$cna.signal
  )

meta$classification <- ifelse(
  meta$celltype_update != "epithelial",
  "non_epithelial_ref",
  NA_character_
)

ref_df <- meta %>% filter(celltype_update != "epithelial")
epi_df <- meta %>% filter(celltype_update == "epithelial")

mean_cor <- NA_real_
sd_cor <- NA_real_
mean_sig <- NA_real_
sd_sig <- NA_real_
thr_cor <- NA_real_
thr_sig <- NA_real_

if (nrow(epi_df) > 0 && nrow(ref_df) > 0) {
  mean_cor <- mean(ref_df$cna_cor, na.rm = TRUE)
  sd_cor <- sd(ref_df$cna_cor, na.rm = TRUE)
  mean_sig <- mean(ref_df$cna_signal, na.rm = TRUE)
  sd_sig <- sd(ref_df$cna_signal, na.rm = TRUE)

  thr_cor <- mean_cor + sd_k_cor * sd_cor
  thr_sig <- mean_sig + sd_k_sig * sd_sig

  meta$classification[meta$celltype_update == "epithelial"] <- ifelse(
    epi_df$cna_cor > thr_cor & epi_df$cna_signal > thr_sig,
    "cna_malignant",
    ifelse(
      epi_df$cna_cor > thr_cor | epi_df$cna_signal > thr_sig,
      "cna_unresolved",
      "cna_non_malignant"
    )
  )
}

tmdata@meta.data <- bind_cols(
  tmdata@meta.data,
  meta[rownames(tmdata@meta.data), c("cna_cor", "cna_signal", "classification"), drop = FALSE]
)

epi <- subset(tmdata, celltype_update == "epithelial")
ncells <- ncol(epi)

if (ncells > 30) {
  nfeat <- min(3000, max(2, nrow(epi)))
  epi <- FindVariableFeatures(epi, nfeatures = nfeat, verbose = FALSE)
  epi <- ScaleData(epi, features = VariableFeatures(epi), verbose = FALSE)
  max_pcs <- max(1, min(length(VariableFeatures(epi)), ncells - 1))
  npcs <- min(50, max(1, round(sqrt(ncells))))
  npcs <- min(npcs, max_pcs)
  epi <- RunPCA(epi, npcs = npcs, verbose = FALSE)
  pc_use <- 1:npcs

  if (ncells <= 10) {
    k_use <- max(1, ncells - 1)
  } else if (ncells <= 500) {
    k_use <- max(5, min(15, ncells - 1))
  } else if (ncells <= 5000) {
    k_use <- 30
  } else if (ncells <= 50000) {
    k_use <- 30
  } else {
    k_use <- min(50, round(sqrt(ncells)))
  }

  nn_method <- if (ncells > 20000) "annoy" else "rann"
  resolution <- if (ncells < 1000) 0.8 else if (ncells < 5000) 1 else 1.2

  epi <- FindNeighbors(
    epi,
    dims = pc_use,
    k.param = max(1, k_use),
    nn.method = nn_method,
    annoy.metric = "euclidean",
    verbose = FALSE
  )
  if (ncells >= 3) {
    epi <- FindClusters(epi, resolution = resolution, algorithm = 1, verbose = FALSE)
  }
  epi <- RunUMAP(
    epi,
    dims = pc_use,
    n.neighbors = max(2, k_use),
    min.dist = if (ncells > 5000) 0.15 else 0.3,
    spread = if (ncells > 5000) 1.5 else 1,
    verbose = FALSE
  )
} else {
  message("Not enough cells, assigning same X clusters")
  epi$seurat_clusters <- rep("X", ncol(epi))
  Idents(epi) <- "seurat_clusters"
}

tab <- table(epi$classification, epi$seurat_clusters)
clusters <- colnames(tab)
cluster_labels <- setNames(rep(NA_character_, length(clusters)), clusters)

for (cl in clusters) {
  malignant <- if ("cna_malignant" %in% rownames(tab)) tab["cna_malignant", cl] else 0
  non_malignant <- if ("cna_non_malignant" %in% rownames(tab)) tab["cna_non_malignant", cl] else 0
  total <- sum(tab[, cl])

  pct_non_malignant <- non_malignant / total
  pct_malignant <- malignant / total

  if (pct_non_malignant > 0.50) {
    cluster_labels[cl] <- "non_malignant_clus"
  } else if (pct_malignant > 0.50) {
    cluster_labels[cl] <- "malignant_clus"
  } else {
    cluster_labels[cl] <- "unresolved_clus"
  }
}

epi$malignant_clus <- as.character(cluster_labels[as.character(epi$seurat_clusters)])

malignant_ref <- colnames(epi)[
  epi$malignant_clus == "malignant_clus" & epi$classification == "cna_malignant"
]
nonmalignant_ref <- colnames(epi)[
  epi$malignant_clus == "non_malignant_clus" & epi$classification == "cna_non_malignant"
]

if (length(malignant_ref) > 50 && length(nonmalignant_ref) > 50) {
  cs_genes <- character(0)
  epi$signature <- rep("sig_unresolved", ncol(epi))
  epi$signature[malignant_ref] <- "sig_malignant_ref"
  epi$signature[nonmalignant_ref] <- "sig_non_malignant_ref"
  Idents(epi) <- "signature"

  dge <- FindMarkers(
    epi,
    ident.1 = "sig_malignant_ref",
    ident.2 = "sig_non_malignant_ref",
    assay = "RNA",
    logfc.threshold = 1,
    min.pct = 0.5,
    only.pos = TRUE
  )

  cs_genes <- rownames(dge[dge$p_val_adj < 0.01 & dge$pct.2 < 0.1, , drop = FALSE])
  if (length(cs_genes) > 0) {
    saveRDS(cs_genes, signature_path)
  } else {
    message("Not enough significant cancer signatures")
  }
} else {
  cs_genes <- character(0)
  message("Not enough reference cells for cancer signature scoring")
}

cluster_summary <- lapply(names(cluster_labels), function(cl) {
  cl_cells <- epi$seurat_clusters == cl
  total <- sum(cl_cells, na.rm = TRUE)
  data.frame(
    sample = sample,
    reference_batch = reference_batch,
    seurat_clusters = cl,
    malignant_clus = unname(cluster_labels[[cl]]),
    n_cells = total,
    cna_malignant = count_value(epi$classification[cl_cells], "cna_malignant"),
    cna_unresolved = count_value(epi$classification[cl_cells], "cna_unresolved"),
    cna_non_malignant = count_value(epi$classification[cl_cells], "cna_non_malignant"),
    pct_cna_malignant = if (total > 0) count_value(epi$classification[cl_cells], "cna_malignant") / total else NA_real_,
    pct_cna_unresolved = if (total > 0) count_value(epi$classification[cl_cells], "cna_unresolved") / total else NA_real_,
    pct_cna_non_malignant = if (total > 0) count_value(epi$classification[cl_cells], "cna_non_malignant") / total else NA_real_,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()
write_cluster_summary(cluster_summary)

malignant_clus_levels <- c("malignant_clus", "non_malignant_clus", "unresolved_clus")
sample_summary <- data.frame(
  sample = sample,
  status = "ok",
  reference_batch = reference_batch,
  reference_path = reference_path,
  reference_types = paste(reference$ref_ct, collapse = "|"),
  filtered_cells = ncol(tmdata),
  epithelial_cells = ncol(epi),
  shared_genes = length(common_genes),
  reference_cells_added = length(ref),
  mean_cor = mean_cor,
  sd_cor = sd_cor,
  thr_cor = thr_cor,
  mean_sig = mean_sig,
  sd_sig = sd_sig,
  thr_sig = thr_sig,
  cna_malignant_cells = count_value(epi$classification, "cna_malignant"),
  cna_unresolved_cells = count_value(epi$classification, "cna_unresolved"),
  cna_non_malignant_cells = count_value(epi$classification, "cna_non_malignant"),
  malignant_clus_cells = count_value(epi$malignant_clus, "malignant_clus"),
  non_malignant_clus_cells = count_value(epi$malignant_clus, "non_malignant_clus"),
  unresolved_clus_cells = count_value(epi$malignant_clus, "unresolved_clus"),
  malignant_clus_n = sum(unique(cluster_summary$seurat_clusters[cluster_summary$malignant_clus == "malignant_clus"]) != "", na.rm = TRUE),
  non_malignant_clus_n = sum(unique(cluster_summary$seurat_clusters[cluster_summary$malignant_clus == "non_malignant_clus"]) != "", na.rm = TRUE),
  unresolved_clus_n = sum(unique(cluster_summary$seurat_clusters[cluster_summary$malignant_clus == "unresolved_clus"]) != "", na.rm = TRUE),
  malignant_ref_cells = length(malignant_ref),
  nonmalignant_ref_cells = length(nonmalignant_ref),
  signature_gene_count = length(cs_genes),
  stringsAsFactors = FALSE
)

if (sample_summary$epithelial_cells > 0) {
  sample_summary$pct_cna_malignant <- sample_summary$cna_malignant_cells / sample_summary$epithelial_cells
  sample_summary$pct_cna_unresolved <- sample_summary$cna_unresolved_cells / sample_summary$epithelial_cells
  sample_summary$pct_cna_non_malignant <- sample_summary$cna_non_malignant_cells / sample_summary$epithelial_cells
  sample_summary$pct_malignant_clus <- sample_summary$malignant_clus_cells / sample_summary$epithelial_cells
  sample_summary$pct_non_malignant_clus <- sample_summary$non_malignant_clus_cells / sample_summary$epithelial_cells
  sample_summary$pct_unresolved_clus <- sample_summary$unresolved_clus_cells / sample_summary$epithelial_cells
} else {
  sample_summary$pct_cna_malignant <- NA_real_
  sample_summary$pct_cna_unresolved <- NA_real_
  sample_summary$pct_cna_non_malignant <- NA_real_
  sample_summary$pct_malignant_clus <- NA_real_
  sample_summary$pct_non_malignant_clus <- NA_real_
  sample_summary$pct_unresolved_clus <- NA_real_
}

write_sample_summary(sample_summary)
saveRDS(epi, file.path("by_samples", sample, paste0(sample, "_epi.rds")))
write_status("ok")
