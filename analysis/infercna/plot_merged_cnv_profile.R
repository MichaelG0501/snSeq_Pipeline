####################
# Analysis registry:
#   Status: active terminal figure-generation
#   Script: analysis/infercna/plot_merged_cnv_profile.R
#   Methodology: analysis/methodology/infercna/plot_merged_cnv_profile_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Output tiers: sn_outs/infercna/cnv_profile/{intermediate,tables,figures,logs,reports}
####################

####################
# plot_merged_cnv_profile.R
#
# Plot a merged InferCNA profile from completed final epithelial malignancy
# objects. This script requires step-6 `_epi_f.rds` files with a `malignancy`
# metadata column; it does not recreate or fallback to older malignancy logic.
#
# Input:
#   sn_outs/by_samples/<sample>/<sample>_epi_f.rds
#   sn_outs/by_samples/<sample>/<sample>_outs.rds
#   /rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt
#
# Output:
#   sn_outs/infercna/cnv_profile/figures/snseq_merged_cnv_profile.pdf
#   sn_outs/infercna/cnv_profile/logs/plot_merged_cnv_profile.log
#
# Usage:
#   Rscript analysis/infercna/plot_merged_cnv_profile.R
####################
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(Seurat)
  library(infercna)
  library(ggplot2)
})

source("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/analysis/lib/config.R")
source(file.path(ANALYSIS_DIR, "lib", "logging.R"))

setwd(SN_OUTS_DIR)
output_dirs <- ensure_output_dirs("infercna/cnv_profile")
run_summary <- start_run_summary(
  script = "analysis/infercna/plot_merged_cnv_profile.R",
  inputs = c(
    file.path(SN_OUTS_DIR, "by_samples/<sample>/<sample>_epi_f.rds"),
    file.path(SN_OUTS_DIR, "by_samples/<sample>/<sample>_outs.rds"),
    "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"
  ),
  outputs = file.path(output_dirs["figures"], "snseq_merged_cnv_profile.pdf"),
  parameters = list(stage = "terminal_cnv_profile", requires_final_malignancy = TRUE)
)

sample_dirs <- list.dirs(path = "by_samples/", full.names = FALSE, recursive = FALSE)

all_meta <- list()
all_outs <- list()
all_ref_meta <- list()
all_ref_outs <- list()

# 10x Gene Order
gene_order <- read.table(
  "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt",
  header = FALSE, col.names = c("gene_id", "chromosome", "start", "end")
)

message("Collecting epithelial and reference cells from samples...")

for (sample in sample_dirs) {
  
  epi_file_f <- file.path("by_samples", sample, paste0(sample, "_epi_f.rds"))
  outs_file  <- file.path("by_samples", sample, paste0(sample, "_outs.rds"))
  
  if (!file.exists(outs_file)) next
  
  if (!file.exists(epi_file_f)) next
  epi <- readRDS(epi_file_f)
  
  outs <- readRDS(outs_file)
  
  # Epithelial subset
  meta_epi <- epi@meta.data
  common_cells <- intersect(colnames(outs), rownames(meta_epi))
  if (length(common_cells) == 0) next
  
  meta_epi <- meta_epi[common_cells, ]
  outs_epi <- outs[, common_cells, drop = FALSE]
  
  if (!"malignancy" %in% colnames(meta_epi)) {
    stop("Missing final `malignancy` metadata in ", epi_file_f, ". Run step 6 before this plotting script.")
  }
  
  # rename malignancy levels
  meta_epi$status <- plyr::revalue(
    meta_epi$malignancy,
    c(
      "malignant_level_1"        = "malignant",
      "malignant_level_2"        = "malignant",
      "non_malignant_level_1"    = "unresolved",
      "non_malignant_level_2"    = "unresolved",
      "malignant_unresolved"     = "unresolved",
      "non_malignant_unresolved" = "unresolved",
      "unresolved"               = "unresolved"
    ),
    warn_missing = FALSE
  )
  
  # Store study as reference_batch
  meta_epi$study <- meta_epi$reference_batch
  
  all_meta[[sample]] <- meta_epi
  all_outs[[sample]] <- outs_epi
  
  # Collect some reference cells (non-epithelial) if they exist in outs
  ref_cells <- setdiff(colnames(outs), rownames(meta_epi))
  if (length(ref_cells) > 0) {
    set.seed(42)
    ref_sub <- sample(ref_cells, min(length(ref_cells), 50))
    
    all_ref_outs[[sample]] <- outs[, ref_sub, drop = FALSE]
    all_ref_meta[[sample]] <- data.frame(
        status = "reference",
        orig.ident = unique(meta_epi$orig.ident)[1],
        cs_score = 0, 
        cc_score = 0, 
        study = unique(meta_epi$reference_batch)[1],
        row.names = ref_sub,
        stringsAsFactors = FALSE
    )
  }
  
  message("Finished processing sample: ", sample)
}

# 2. Merge all
message("Merging matrices...")
common_rows <- Reduce(intersect, lapply(c(all_outs, all_ref_outs), rownames))
all_outs_i <- lapply(all_outs, function(m) m[common_rows, , drop = FALSE])
all_ref_outs_i <- lapply(all_ref_outs, function(m) m[common_rows, , drop = FALSE])

common_cols <- c("orig.ident", "status", "cs_score", "cc_score", "study")
meta_epi_merged <- do.call(rbind, lapply(all_meta, function(m) m[, common_cols, drop = FALSE]))
meta_ref_merged <- do.call(rbind, all_ref_meta)

# Combine epithelial and reference
set.seed(1)
n_epi_target <- 5000
n_ref_target <- 800

keep_epi <- sample(seq_len(nrow(meta_epi_merged)), min(nrow(meta_epi_merged), n_epi_target))
keep_ref <- sample(seq_len(nrow(meta_ref_merged)), min(nrow(meta_ref_merged), n_ref_target))

meta <- rbind(meta_epi_merged[keep_epi, ], meta_ref_merged[keep_ref, ])
outs_merged <- cbind(
    do.call(cbind, all_outs_i)[, rownames(meta_epi_merged)[keep_epi]],
    do.call(cbind, all_ref_outs_i)[, rownames(meta_ref_merged)[keep_ref]]
)

# Fix study names (shorten if needed)
meta$study <- sub("_reference$", "", meta$study)

study_levels <- unique(meta$study)
study_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set3"))(length(study_levels)),
  study_levels
)

chrom_levels <- c(paste0("chr", 1:22), "chrX", "chrY")
common_genes <- intersect(rownames(outs_merged), gene_order$gene_id)

go <- gene_order %>%
  filter(gene_id %in% common_genes, chromosome %in% chrom_levels) %>%
  mutate(chromosome = factor(chromosome, levels = chrom_levels)) %>%
  arrange(chromosome, start)

outs_merged <- outs_merged[go$gene_id, , drop = FALSE]
stopifnot(identical(rownames(outs_merged), go$gene_id))

## --------------------------- BINNING ----------------------------- ##
message("Binning...")
bin_size <- 100L
go <- go %>%
  group_by(chromosome) %>%
  mutate(
    g_rank = row_number(),
    bin_in_chr = ((g_rank - 1L) %/% bin_size) + 1L,
    bin_key = paste(chromosome, bin_in_chr, sep = "_")
  ) %>%
  ungroup()

ordered_bin_keys <- unique(go$bin_key)
bins_idx <- split(seq_len(nrow(go)), factor(go$bin_key, levels = ordered_bin_keys))

binned_mat <- do.call(rbind, lapply(bins_idx, function(ix) colMeans(outs_merged[ix, , drop = FALSE])))
rownames(binned_mat) <- names(bins_idx)

row_chr_labels <- sub("_.*$", "", rownames(binned_mat))
row_chr <- factor(row_chr_labels, levels = unique(row_chr_labels))

top_ha <- HeatmapAnnotation(
  cancer_signature = as.numeric(meta$cs_score),
  cell_cycling     = as.numeric(meta$cc_score),
  batch            = factor(meta$study),
  
  col = list(
    cancer_signature = colorRamp2(c(0, 1, 1.5), c("white", "grey80", "black")),
    cell_cycling     = colorRamp2(c(0, 1, 1.5), c("white", "grey80", "black")),
    batch            = study_colors
  ),
  
  annotation_name_side   = "left",
  annotation_name_gp     = gpar(fontsize = 9),
  annotation_name_offset = unit(2, "mm"),
  annotation_height      = unit(c(4, 4, 4), "mm"),
  
  show_legend = c(cancer_signature = TRUE, cell_cycling = TRUE, batch = TRUE),
  annotation_legend_param = list(
    cancer_signature = list(title = "CS score"),
    cell_cycling     = list(title = "CC score"),
    batch            = list(title = "Batch")
  )
)

chr_used <- levels(droplevels(row_chr))
base_cols <- c(brewer.pal(12, "Paired"),
               brewer.pal(8,  "Dark2"),
               brewer.pal(9,  "Set1"),
               brewer.pal(12, "Set3"))
chr_cols <- setNames(base_cols[seq_along(chr_used)], chr_used)

left_chr_bar <- rowAnnotation(
  chr = row_chr,
  col = list(chr = chr_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  gp = gpar(col = NA),
  width = unit(4, "mm")
)

chr_bounds <- which(head(row_chr_labels, -1L) != tail(row_chr_labels, -1L))

class_order <- c("malignant", "unresolved", "reference")
meta$classification <- factor(meta$status, levels = class_order)
col_split <- meta$classification

ht <- Heatmap(
  binned_mat,
  name = "CNV",
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  column_split = col_split, 
  column_title_rot = 0, 
  cluster_column_slices = TRUE,
  show_column_dend = FALSE,
  column_gap = unit(2, "mm"),
  
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation  = top_ha,
  left_annotation = left_chr_bar,
  
  row_split = row_chr,
  row_gap = unit(0, "mm"),
  row_title_rot  = 0,
  rect_gp = gpar(col = NA),
  border = NA,
  
  layer_fun = function(j, i, x, y, w, h, fill) {
    hits <- intersect(i, chr_bounds)
    if (length(hits)) {
      id <- match(hits, i)
      yy <- y[id] - h[id]/2
      grid.segments(
        x0 = unit(0, "npc"), x1 = unit(1, "npc"),
        y0 = yy, y1 = yy,
        gp = gpar(col = "black", lwd = 1)
      )
    }
  }
)

pdf(file.path(output_dirs["figures"], "snseq_merged_cnv_profile.pdf"), width = 10, height = 8)
draw(ht)
dev.off()

run_summary <- finish_run_summary(run_summary, status = "ok")
write_run_summary(
  run_summary,
  file.path(output_dirs["logs"], "plot_merged_cnv_profile.log")
)

message("Done! Saved to sn_outs/infercna/cnv_profile/figures/snseq_merged_cnv_profile.pdf")
