####################
# Auto_states_unresolved_pan_cancer_noreg.R
# Unresolved-cell pan-cancer subclassification for snRNA-seq (noreg mode).
# Assignment based on top raw UCell 3CA score (for barplot).
# Unsupervised clustering for heatmap.
####################

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs")

# 1. Load data
message("Loading data...")
tmdata_all <- readRDS("snSeq_malignant_epi.rds")
state_noreg <- readRDS("Auto_topmp_v2_noreg_states_B.rds")
ucell_3ca <- readRDS("UCell_3CA_MPs.rds")
cc_score <- readRDS("snSeq_malignant_epi_cc_score.rds")

# 2. Align cells and define MPs
common_cells <- intersect(Cells(tmdata_all), rownames(ucell_3ca))
common_cells <- intersect(common_cells, names(state_noreg))
common_cells <- intersect(common_cells, names(cc_score))

tmdata_all <- tmdata_all[, common_cells]
ucell_3ca <- ucell_3ca[common_cells, , drop = FALSE]
state_noreg <- state_noreg[common_cells]
cc_score <- cc_score[common_cells]

CC_FIXED <- c(
  "X3CA_mp_1.Cell.Cycle...G2.M",
  "X3CA_mp_2.Cell.Cycle...G1.S",
  "X3CA_mp_3.Cell.Cylce.HMG.rich",
  "X3CA_mp_4.Chromatin",
  "X3CA_mp_5.Cell.cycle.single.nucleus"
)
pan_cancer_mps <- setdiff(colnames(ucell_3ca), CC_FIXED)

# 3. Logic for Unresolved cells
unresolved_cells <- names(state_noreg)[state_noreg == "Unresolved"]
if (length(unresolved_cells) == 0) {
  stop("No unresolved cells found.")
}

sub_scores_raw <- ucell_3ca[unresolved_cells, pan_cancer_mps, drop = FALSE]

# Top MP assigned for Barplot
top_mp_assigned <- pan_cancer_mps[max.col(sub_scores_raw, ties.method = "first")]
names(top_mp_assigned) <- unresolved_cells

# Clean names for plotting
clean_mp_name <- function(x) {
  x <- gsub("^X3CA_", "3CA_", x)
  gsub("\\.", " ", x)
}

# 4. Bar plot: Top MP assignment counts (No stacking)
plot_df <- data.frame(
  cell = unresolved_cells,
  top_mp = top_mp_assigned,
  stringsAsFactors = FALSE
) %>%
  mutate(top_mp_clean = clean_mp_name(top_mp))

# Order MPs by total abundance
mp_order <- plot_df %>%
  count(top_mp_clean) %>%
  arrange(desc(n)) %>%
  pull(top_mp_clean)

plot_df$top_mp_clean <- factor(plot_df$top_mp_clean, levels = mp_order)

p_bar <- ggplot(plot_df, aes(x = top_mp_clean)) +
  geom_bar(fill = "steelblue", color = "black", linewidth = 0.1) +
  labs(title = "Unresolved Cells: Pan-Cancer MP Assignment",
       subtitle = paste0("Top raw UCell score (Total cells: ", length(unresolved_cells), ")"),
       x = "Pan-Cancer MP",
       y = "Cell Count") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("Auto_topmp_v2_noreg_unresolved_pan_cancer_barplot.pdf", p_bar, width = 12, height = 8)

# 5. Unsupervised subclassification for Heatmap
z_normalise <- function(mat, sample_var, study_var) {
  clust_df <- as.data.frame(mat)
  clust_df$.cell <- rownames(mat)
  clust_df$.sample <- sample_var[rownames(mat)]
  clust_df$.study <- study_var[rownames(mat)]

  study_sd <- clust_df %>%
    group_by(.study) %>%
    summarise(across(all_of(colnames(mat)), ~ sd(.x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames(".study") %>%
    as.matrix()
  study_sd[is.na(study_sd) | study_sd == 0] <- 1

  clust_centered <- clust_df %>%
    group_by(.sample) %>%
    mutate(across(all_of(colnames(mat)), ~ .x - mean(.x, na.rm = TRUE))) %>%
    ungroup()

  mp_adj <- as.matrix(clust_centered[, colnames(mat), drop = FALSE])
  rownames(mp_adj) <- clust_centered$.cell
  for (mp in colnames(mp_adj)) {
    mp_adj[, mp] <- mp_adj[, mp] / study_sd[clust_centered$.study, mp]
  }
  mp_adj[!is.finite(mp_adj)] <- 0
  mp_adj
}

mp_adj <- z_normalise(as.matrix(sub_scores_raw), 
                      tmdata_all$orig.ident[unresolved_cells], 
                      tmdata_all$reference_batch[unresolved_cells])

# Clustering logic
mp_matrix <- t(mp_adj)
temp_obj <- CreateSeuratObject(counts = mp_matrix, assay = "MPs")
temp_obj <- ScaleData(temp_obj, assay = "MPs", layer = "counts", features = rownames(temp_obj), verbose = FALSE)
n_pcs <- min(30, nrow(mp_matrix) - 1)
temp_obj <- RunPCA(temp_obj, features = rownames(temp_obj), npcs = n_pcs, verbose = FALSE)
temp_obj <- FindNeighbors(temp_obj, reduction = "pca", dims = 1:n_pcs, graph.name = "MPs_snn", verbose = FALSE)
temp_obj <- FindClusters(temp_obj, graph.name = "MPs_snn", resolution = 1, verbose = FALSE)

subclass <- paste0("State_", as.numeric(as.character(temp_obj$seurat_clusters)) + 1)
names(subclass) <- unresolved_cells

# Sample cells for heatmap
set.seed(42)
max_cells_heatmap <- 6000
cells_to_plot <- if(length(unresolved_cells) > max_cells_heatmap) sample(unresolved_cells, max_cells_heatmap) else unresolved_cells

sub_scores_plot <- t(mp_adj[cells_to_plot, , drop = FALSE])
rownames(sub_scores_plot) <- clean_mp_name(rownames(sub_scores_plot))
split_vec <- factor(subclass[cells_to_plot])

study_vals <- tmdata_all$reference_batch[cells_to_plot]
study_levels <- unique(as.character(tmdata_all$reference_batch))
study_cols <- setNames(DiscretePalette(length(study_levels), palette = "polychrome"), study_levels)

subtype_cols <- setNames(hue_pal()(length(levels(split_vec))), levels(split_vec))

col_ann <- HeatmapAnnotation(
  Subclass = split_vec,
  CNA = as.character(tmdata_all$classification[cells_to_plot]),
  CC_score = cc_score[cells_to_plot],
  Study = study_vals,
  col = list(
    Subclass = subtype_cols,
    CNA = c(cna_malignant = "black", cna_unresolved = "grey70", cna_non_malignant = "white"),
    CC_score = colorRamp2(c(0, max(cc_score, na.rm = TRUE)), c("white", "darkgreen")),
    Study = study_cols
  ),
  annotation_name_side = "left",
  na_col = "white"
)

lim <- as.numeric(quantile(abs(sub_scores_plot), 0.98, na.rm = TRUE))
if (is.na(lim) || lim == 0) lim <- 1

ht <- Heatmap(
  sub_scores_plot,
  name = "Adj score",
  col = colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3")),
  top_annotation = col_ann,
  column_split = split_vec,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_rows = "ward.D2",
  show_row_dend = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8, fontface = "italic"),
  use_raster = TRUE,
  raster_quality = 5,
  column_title = "Unresolved Pan-Cancer Subclassification (Unsupervised Clustering)"
)

pdf("Auto_topmp_v2_noreg_unresolved_pan_cancer_heatmap.pdf", width = 18, height = 10, useDingbats = FALSE)
draw(ht, merge_legend = TRUE)
dev.off()

# Save summary
write.csv(plot_df %>% count(top_mp_clean, name = "cells"), "Auto_topmp_v2_noreg_unresolved_pan_cancer_summary.csv", row.names = FALSE)
saveRDS(data.frame(cell = unresolved_cells, subclass = subclass, top_mp = top_mp_assigned), "Auto_topmp_v2_noreg_unresolved_pan_cancer_results.rds")

message("Finished. Outputs saved to sn_outs/")
