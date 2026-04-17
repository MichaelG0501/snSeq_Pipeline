suppressPackageStartupMessages({
  library("dplyr")
  library("Matrix")
  library("parallel")
  library("Seurat")
  library("UCell")
})

set.seed(1)

project_dir <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline"
out_dir <- file.path(project_dir, "sn_outs")
ref_outs_dir <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs"
setwd(out_dir)

dir.create("Metaprogrammes_Results", recursive = TRUE, showWarnings = FALSE)

read_status <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  trimws(readLines(path, warn = FALSE, n = 1))
}

extract_layer_matrix_subset <- function(seurat_obj, genes, assay = "RNA", layer_prefix = "data") {
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  if (length(genes) == 0) {
    stop("No genes provided for layer extraction.")
  }

  assay_obj <- seurat_obj[[assay]]
  layer_names <- Layers(assay_obj)
  target_layers <- layer_names[grepl(paste0("^", layer_prefix, "(\\.|$)"), layer_names)]

  if (length(target_layers) == 0) {
    stop("No layers found with prefix '", layer_prefix, "' in assay ", assay)
  }

  layer_mats <- lapply(target_layers, function(layer_name) {
    layer_mat <- LayerData(seurat_obj, assay = assay, layer = layer_name)
    common_genes <- intersect(genes, rownames(layer_mat))
    out_mat <- Matrix(
      0,
      nrow = length(genes),
      ncol = ncol(layer_mat),
      sparse = TRUE,
      dimnames = list(genes, colnames(layer_mat))
    )
    if (length(common_genes) > 0) {
      out_mat[match(common_genes, genes), ] <- layer_mat[common_genes, , drop = FALSE]
    }
    out_mat
  })

  merged_mat <- do.call(cbind, layer_mats)
  merged_mat[, Cells(seurat_obj), drop = FALSE]
}

manifest_path <- file.path(out_dir, "sample_manifest.csv")
if (!file.exists(manifest_path)) {
  stop("Missing sample manifest: ", manifest_path)
}

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
sample_lookup <- manifest %>%
  select(sample, input_source, technology, reference_batch, cells)

status_rows <- lapply(manifest$sample, function(sample) {
  epi_f_path <- file.path("by_samples", sample, paste0(sample, "_epi_f.rds"))
  data.frame(
    sample = sample,
    expr_filter_status = read_status(file.path("by_samples", sample, "expr_filter_status.txt")),
    infercna_status = read_status(file.path("by_samples", sample, "infercna_status.txt")),
    malignancy_status = read_status(file.path("by_samples", sample, "malignancy_status.txt")),
    epi_f_present = file.exists(epi_f_path),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

eligible_samples <- status_rows$sample[
  status_rows$expr_filter_status == "ok" &
    status_rows$infercna_status == "ok" &
    status_rows$malignancy_status == "ok" &
    status_rows$epi_f_present
]

if (length(eligible_samples) == 0) {
  stop("No completed malignancy outputs available for malignant epithelial merge.")
}

local_mp_path <- file.path("Metaprogrammes_Results", "geneNMF_metaprograms_nMP_19.rds")
source_mp_path <- file.path(ref_outs_dir, "Metaprogrammes_Results", "geneNMF_metaprograms_nMP_19.rds")

if (file.exists(local_mp_path)) {
  geneNMF.metaprograms <- readRDS(local_mp_path)
} else {
  geneNMF.metaprograms <- readRDS(source_mp_path)
  saveRDS(geneNMF.metaprograms, local_mp_path)
}

writeLines(
  source_mp_path,
  con = file.path("Metaprogrammes_Results", "geneNMF_metaprograms_nMP_19_source.txt")
)

mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}

required_meta_cols <- c("orig.ident", "reference_batch", "technology", "classification", "malignancy")

tmdata_malignant <- list()
sample_rows <- list()

for (sample in eligible_samples) {
  sample_path <- file.path("by_samples", sample, paste0(sample, "_epi_f.rds"))
  tmdata <- readRDS(sample_path)
  meta <- tmdata@meta.data

  missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))
  if (length(missing_meta_cols) > 0) {
    stop(
      "Missing metadata in ", sample_path, ": ",
      paste(missing_meta_cols, collapse = ", ")
    )
  }
 
  malignant_keep <- meta$malignancy %in% c("malignant_level_1", "malignant_level_2")
  malignant_cells <- colnames(tmdata)[which(malignant_keep)]

  sample_rows[[sample]] <- data.frame(
    sample = sample,
    reference_batch = unique(meta$reference_batch)[1],
    technology = unique(meta$technology)[1],
    epithelial_cells = ncol(tmdata),
    malignant_cells = length(malignant_cells),
    malignant_level_1_cells = sum(meta$malignancy == "malignant_level_1", na.rm = TRUE),
    malignant_level_2_cells = sum(meta$malignancy == "malignant_level_2", na.rm = TRUE),
    cna_malignant_cells = sum(meta$classification == "cna_malignant", na.rm = TRUE),
    cna_unresolved_cells = sum(meta$classification == "cna_unresolved", na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  if (length(malignant_cells) == 0) {
    next
  }

  tmdata_malignant[[sample]] <- subset(tmdata, cells = malignant_cells)
}

sample_summary <- bind_rows(sample_rows) %>%
  left_join(sample_lookup, by = "sample")

write.csv(
  sample_summary,
  file.path(out_dir, "snSeq_malignant_epi_sample_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  status_rows,
  file.path(out_dir, "snSeq_malignant_epi_status_input.csv"),
  row.names = FALSE,
  quote = TRUE
)

if (length(tmdata_malignant) == 0) {
  stop("No malignant epithelial cells were found in completed step-6 outputs.")
}

all_cells <- unlist(lapply(tmdata_malignant, Cells), use.names = FALSE)
if (anyDuplicated(all_cells) > 0) {
  original_names <- names(tmdata_malignant)
  tmdata_malignant <- lapply(original_names, function(sample) {
    RenameCells(tmdata_malignant[[sample]], add.cell.id = sample)
  })
  names(tmdata_malignant) <- original_names
}

merged_obj <- tmdata_malignant[[1]]
if (length(tmdata_malignant) > 1) {
  merged_obj <- merge(
    x = tmdata_malignant[[1]],
    y = tmdata_malignant[-1],
    project = "snSeq_malignant_epi"
  )
}

saveRDS(merged_obj, file.path(out_dir, "snSeq_malignant_epi.rds"))
saveRDS(merged_obj@meta.data, file.path(out_dir, "snSeq_malignant_epi_meta.rds"))

mp.genes <- lapply(mp.genes, function(genes) {
  intersect(unique(genes), rownames(merged_obj))
})
mp.genes <- mp.genes[lengths(mp.genes) > 0]

if (length(mp.genes) == 0) {
  stop("No retained scRef metaprogram genes overlap the snRNA-seq malignant epithelial object.")
}

MP_df <- read.csv(
  "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv",
  check.names = FALSE
)
MP_list <- as.list(MP_df)
MP_list <- lapply(MP_list, function(genes) {
  unique(genes[genes != "" & !is.na(genes)])
})
names(MP_list) <- make.names(sub("^MP", "3CA_mp_", names(MP_list)))
MP_list <- lapply(MP_list, function(genes) {
  intersect(genes, rownames(merged_obj))
})
MP_list <- MP_list[lengths(MP_list) > 0]

if (length(MP_list) == 0) {
  stop("No 3CA metaprogram genes overlap the snRNA-seq malignant epithelial object.")
}

ncores_use <- min(8L, max(1L, parallel::detectCores()))
cell_cycle_genes <- read.csv(
  "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv",
  stringsAsFactors = FALSE
)[, 1:3]
cc_consensus_genes <- unique(cell_cycle_genes$Gene[cell_cycle_genes$Consensus == 1])
scoring_genes <- unique(c(unlist(mp.genes), unlist(MP_list), cc_consensus_genes))
scoring_expr <- extract_layer_matrix_subset(merged_obj, scoring_genes, assay = "RNA", layer_prefix = "data")

ucell_scores <- ScoreSignatures_UCell(
  scoring_expr,
  features = mp.genes,
  ncores = ncores_use,
  name = ""
)
saveRDS(ucell_scores, file.path("Metaprogrammes_Results", "UCell_nMP19_filtered.rds"))

write.csv(
  data.frame(
    mp = names(mp.genes),
    genes_in_signature = lengths(mp.genes),
    stringsAsFactors = FALSE
  ),
  file.path("Metaprogrammes_Results", "UCell_nMP19_filtered_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)

ucell_3ca <- ScoreSignatures_UCell(
  scoring_expr,
  features = MP_list,
  ncores = ncores_use,
  name = ""
)
saveRDS(ucell_3ca, file.path(out_dir, "UCell_3CA_MPs.rds"))

write.csv(
  data.frame(
    mp = names(MP_list),
    genes_in_signature = lengths(MP_list),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "UCell_3CA_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)

cc_consensus <- intersect(cc_consensus_genes, rownames(scoring_expr))
if (length(cc_consensus) == 0) {
  stop("No consensus cell-cycle genes overlap the snRNA-seq malignant epithelial object.")
}

cc_top50 <- names(sort(Matrix::rowMeans(scoring_expr[cc_consensus, , drop = FALSE]), decreasing = TRUE))[
  seq_len(min(50, length(cc_consensus)))
]
cc_score <- Matrix::colMeans(scoring_expr[cc_top50, , drop = FALSE])
names(cc_score) <- colnames(scoring_expr)

saveRDS(cc_top50, file.path(out_dir, "snSeq_malignant_epi_cc_top50.rds"))
saveRDS(cc_score, file.path(out_dir, "snSeq_malignant_epi_cc_score.rds"))

overall_summary <- data.frame(
  eligible_samples = length(eligible_samples),
  samples_with_malignant_cells = length(tmdata_malignant),
  malignant_cells_merged = ncol(merged_obj),
  malignant_level_1_cells = sum(merged_obj$malignancy == "malignant_level_1", na.rm = TRUE),
  malignant_level_2_cells = sum(merged_obj$malignancy == "malignant_level_2", na.rm = TRUE),
  batches_present = dplyr::n_distinct(merged_obj$reference_batch),
  retained_scRef_metaprograms = length(mp.genes),
  scored_3ca_metaprograms = length(MP_list),
  stringsAsFactors = FALSE
)

write.csv(
  overall_summary,
  file.path(out_dir, "snSeq_malignant_epi_overall_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)

message(
  "Saved malignant epithelial merge: ",
  ncol(merged_obj),
  " cells across ",
  length(tmdata_malignant),
  " samples."
)
