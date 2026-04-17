####################
# Integrated into snSeq_Pipeline
# Adapted from scRef_Pipeline/analysis/cnv/cnv_subsetting.R
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

# Set working directory to sn_outs
setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs")

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
  epi_file   <- file.path("by_samples", sample, paste0(sample, "_epi.rds"))
  outs_file  <- file.path("by_samples", sample, paste0(sample, "_outs.rds"))
  
  if (!file.exists(outs_file)) next
  
  if (file.exists(epi_file_f)) {
    epi <- readRDS(epi_file_f)
  } else if (file.exists(epi_file)) {
    epi <- readRDS(epi_file)
  } else {
    next
  }
  
  outs <- readRDS(outs_file)
  
  # Epithelial subset
  meta_epi <- epi@meta.data
  common_cells <- intersect(colnames(outs), rownames(meta_epi))
  if (length(common_cells) == 0) next
  
  meta_epi <- meta_epi[common_cells, ]
  outs_epi <- outs[, common_cells, drop = FALSE]
  
  # Ensure cc_score and cs_score exist
  if (!"cc_score" %in% colnames(meta_epi)) meta_epi$cc_score <- 0
  if (!"cs_score" %in% colnames(meta_epi)) meta_epi$cs_score <- 0
  
  # Fallback malignancy logic if not present (logic from Malignancy.R)
  if (!"malignancy" %in% colnames(meta_epi)) {
    message("Malignancy column missing for ", sample, " - applying fallback logic...")
    meta_epi$malignancy <- "unresolved"
    
    # Requirement for fallback: classification and malignant_clus must exist (from Step 5)
    if ("classification" %in% colnames(meta_epi) && "malignant_clus" %in% colnames(meta_epi)) {
      cc_thr <- 1
      cs_thr <- 1
      
      cluster_col <- if ("seurat_clusters" %in% colnames(meta_epi)) "seurat_clusters" else "res.0.8"
      if (cluster_col %in% colnames(meta_epi)) {
        cluster_labels <- tapply(meta_epi$malignant_clus, meta_epi[[cluster_col]], function(x) {
          if (all(x == "malignant_clus", na.rm = TRUE)) "malignant_clus"
          else if (all(x == "non_malignant_clus", na.rm = TRUE)) "non_malignant_clus"
          else "unresolved_clus"
        })
        
        for (cl in names(cluster_labels)) {
          cells_in_cl <- meta_epi[[cluster_col]] == cl
          if (cluster_labels[cl] == "non_malignant_clus") {
             meta_epi$malignancy[cells_in_cl & meta_epi$classification == "cna_non_malignant"] <- "non_malignant_level_1"
             meta_epi$malignancy[cells_in_cl & meta_epi$classification == "cna_unresolved"] <- "non_malignant_level_2"
             meta_epi$malignancy[cells_in_cl & meta_epi$classification == "cna_malignant"] <- "non_malignant_unresolved"
          } else {
             is_malig <- cells_in_cl & (
               meta_epi$classification == "cna_malignant" | 
               meta_epi$cc_score >= cc_thr | 
               meta_epi$cs_score >= cs_thr
             )
             if (sum(is_malig) / sum(cells_in_cl) >= 0.5) {
               meta_epi$malignancy[cells_in_cl & meta_epi$classification == "cna_malignant"] <- "malignant_level_1"
               meta_epi$malignancy[cells_in_cl & meta_epi$classification == "cna_unresolved"] <- "malignant_level_2"
               meta_epi$malignancy[cells_in_cl & meta_epi$classification == "cna_non_malignant"] <- "malignant_unresolved"
             }
          }
        }
      }
    }
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

pdf("snseq_merged_cnv_profile.pdf", width = 10, height = 8)
draw(ht)
dev.off()

message("Done! Saved to sn_outs/snseq_merged_cnv_profile.pdf")
