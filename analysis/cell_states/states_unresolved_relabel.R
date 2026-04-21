####################
# Auto_unresolved_relabel.R
# Relabel unresolved noreg Approach-B cells using the scRef-retained 3CA MPs
# and regenerate the scRef-style proportion, boxplot, and heatmap outputs for
# the snRNA-seq malignant epithelial dataset.
####################

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(patchwork)
library(grid)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs")
task_prefix <- "task4"
out_dir <- file.path(getwd(), paste0(task_prefix, "_unresolved_states"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

####################
# constants
####################
state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive" = c("MP13", "MP12"),
  "SMG-like Metaplasia" = c("MP18", "MP16"),
  "Immune Infiltrating" = c("MP15")
)

group_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "Immune Infiltrating" = "#377EB8",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)

mp_descriptions <- c(
  "MP1"  = "G2M Cell Cycle",
  "MP9"  = "G1S Cell Cycle",
  "MP2"  = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP5"  = "Epithelial IFN Resp.",
  "MP10" = "Columnar Diff.",
  "MP8"  = "Intestinal Diff.",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP7"  = "DNA Damage Repair",
  "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP15" = "Immune Infiltration",
  "MP12" = "Neuro-responsive Epi"
)

cc_mps <- c("MP1", "MP7", "MP9")

retained_3ca_map <- c(
  "X3CA_mp_30.Respiration.1" = "Classic Proliferative",
  "X3CA_mp_12.Protein.maturation" = "3CA_EMT_and_Protein_maturation",
  "X3CA_mp_17.EMT.III" = "3CA_EMT_and_Protein_maturation"
)

####################
# helpers
####################
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

clean_3ca_name <- function(x) {
  x <- gsub("^X3CA_", "3CA_", x)
  gsub("\\.", " ", x)
}

make_cna_palette <- function(cna_vec) {
  cna_levels <- sort(unique(as.character(cna_vec[!is.na(cna_vec)])))
  if (length(cna_levels) == 0) {
    cna_levels <- "unknown"
  }
  cna_palette <- setNames(rep("grey70", length(cna_levels)), cna_levels)
  if ("cna_malignant" %in% names(cna_palette)) cna_palette["cna_malignant"] <- "black"
  if ("cna_unresolved" %in% names(cna_palette)) cna_palette["cna_unresolved"] <- "grey70"
  if ("cna_non_malignant" %in% names(cna_palette)) cna_palette["cna_non_malignant"] <- "white"
  cna_palette
}

####################
# load data
####################
tmdata_all <- readRDS("snSeq_malignant_epi.rds")
state_B <- readRDS("Auto_topmp_v2_noreg_states_B.rds")
mp_adj_noncc <- readRDS("Auto_topmp_v2_noreg_mp_adj.rds")
ucell_3ca <- readRDS("UCell_3CA_MPs.rds")
geneNMF.metaprograms <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds")
ucell_scores <- readRDS("Metaprogrammes_Results/UCell_nMP19_filtered.rds")
cc_score <- readRDS("snSeq_malignant_epi_cc_score.rds")

####################
# MP ordering
####################
mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}
retained_mps <- names(mp.genes)

state_ordered_mps <- unlist(state_groups, use.names = FALSE)
orig_tree_order <- geneNMF.metaprograms$programs.tree$order
orig_clusters <- geneNMF.metaprograms$programs.clusters[orig_tree_order]
orig_order <- paste0("MP", unique(orig_clusters))
orig_order <- orig_order[orig_order %in% retained_mps]
reordered_mps <- c(
  orig_order[orig_order %in% cc_mps],
  state_ordered_mps[state_ordered_mps %in% retained_mps]
)
reordered_mps <- unique(c(reordered_mps, orig_order))
mp_tree_order_names <- reordered_mps

group_order_pos <- sapply(state_groups, function(mps) {
  positions <- match(mps, mp_tree_order_names)
  if (all(is.na(positions))) return(Inf)
  min(positions, na.rm = TRUE)
})
ordered_group_names <- names(sort(group_order_pos))

####################
# align cells and relabel unresolved cells
####################
common_cells <- Reduce(intersect, list(
  Cells(tmdata_all),
  names(state_B),
  rownames(mp_adj_noncc),
  rownames(ucell_3ca),
  rownames(ucell_scores),
  names(cc_score)
))

tmdata_all <- tmdata_all[, common_cells]
state_B <- state_B[common_cells]
mp_adj_noncc <- mp_adj_noncc[common_cells, , drop = FALSE]
ucell_3ca <- ucell_3ca[common_cells, , drop = FALSE]
ucell_scores <- ucell_scores[common_cells, , drop = FALSE]
cc_score <- cc_score[common_cells]

sample_var <- tmdata_all$orig.ident
study_var <- tmdata_all$reference_batch
names(sample_var) <- Cells(tmdata_all)
names(study_var) <- Cells(tmdata_all)

unresolved_cells <- names(state_B)[state_B == "Unresolved"]
state_updated <- state_B

if (length(unresolved_cells) > 0) {
  CC_FIXED <- c(
    "X3CA_mp_1.Cell.Cycle...G2.M",
    "X3CA_mp_2.Cell.Cycle...G1.S",
    "X3CA_mp_3.Cell.Cylce.HMG.rich",
    "X3CA_mp_4.Chromatin",
    "X3CA_mp_5.Cell.cycle.single.nucleus"
  )

  unresolved_3ca <- ucell_3ca[unresolved_cells, setdiff(colnames(ucell_3ca), CC_FIXED), drop = FALSE]
  top_3ca_mp <- colnames(unresolved_3ca)[max.col(unresolved_3ca, ties.method = "first")]
  names(top_3ca_mp) <- unresolved_cells
  relabel_cells <- names(top_3ca_mp)[top_3ca_mp %in% names(retained_3ca_map)]
  state_updated[relabel_cells] <- retained_3ca_map[top_3ca_mp[relabel_cells]]

  coverage_df <- data.frame(
    cell = unresolved_cells,
    mp_label = top_3ca_mp,
    orig.ident = as.character(sample_var[unresolved_cells]),
    study = as.character(study_var[unresolved_cells]),
    stringsAsFactors = FALSE
  ) %>%
    group_by(mp_label) %>%
    summarise(
      n_cells = n(),
      n_samples = n_distinct(orig.ident),
      n_studies = n_distinct(study),
      pct_cells = 100 * n() / length(common_cells),
      retained_from_scref = first(mp_label %in% names(retained_3ca_map)),
      mapped_state = first(ifelse(mp_label %in% names(retained_3ca_map), retained_3ca_map[mp_label], "Unresolved")),
      .groups = "drop"
    ) %>%
    arrange(desc(n_cells))
} else {
  top_3ca_mp <- character(0)
  coverage_df <- data.frame()
}

write.csv(
  coverage_df,
  file.path(out_dir, "Auto_task4_unresolved_relabel_mp_coverage.csv"),
  row.names = FALSE
)

####################
# final state order and colors
####################
retained_3ca_available <- intersect(names(retained_3ca_map), colnames(ucell_3ca))
new_retained_3ca_names <- clean_3ca_name(retained_3ca_available)
new_retained_3ca_names <- setdiff(
  new_retained_3ca_names,
  c("3CA_mp_30 Respiration 1", "3CA_mp_12 Protein maturation", "3CA_mp_17 EMT III")
)
if (any(c("3CA_mp_12 Protein maturation", "3CA_mp_17 EMT III") %in% clean_3ca_name(retained_3ca_available))) {
  new_retained_3ca_names <- unique(c(new_retained_3ca_names, "3CA_EMT_and_Protein_maturation"))
}

state_level_order_updated <- c(
  ordered_group_names,
  new_retained_3ca_names,
  "Unresolved",
  "Hybrid"
)

custom_3ca_cols <- c(
  "3CA_EMT_and_Protein_maturation" = "#666666",
  "3CA_mp_1 Epithelial-1" = "#F781BF",
  "3CA_mp_5 Epithelial-5" = "#A65628",
  "3CA_mp_21 Epithelial-21" = "#FFFF33"
)
remaining_3ca <- setdiff(new_retained_3ca_names, names(custom_3ca_cols))
if (length(remaining_3ca) > 0) {
  new_state_cols <- c(
    custom_3ca_cols[intersect(names(custom_3ca_cols), new_retained_3ca_names)],
    setNames(scales::hue_pal()(length(remaining_3ca)), remaining_3ca)
  )
} else {
  new_state_cols <- custom_3ca_cols[intersect(names(custom_3ca_cols), new_retained_3ca_names)]
}
group_cols_updated <- c(group_cols, new_state_cols)

saveRDS(state_updated, file.path(out_dir, "Auto_task4_unresolved_relabel_states.rds"))
saveRDS(state_updated, "Auto_final_states.rds")

####################
# updated proportion plot
####################
prop_df <- data.frame(
  state = as.character(state_updated[common_cells]),
  study = as.character(tmdata_all$reference_batch[common_cells]),
  stringsAsFactors = FALSE
)
overall <- prop_df %>% count(state) %>% mutate(study = "Overall", pct = 100 * n / sum(n))
per_study <- prop_df %>%
  count(study, state) %>%
  group_by(study) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()
plot_df <- bind_rows(overall, per_study)
plot_df$state <- factor(plot_df$state, levels = state_level_order_updated)

p_bar <- ggplot(plot_df, aes(study, pct, fill = state)) +
  geom_col(color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), position = position_stack(vjust = 0.5), size = 2.2) +
  scale_fill_manual(values = group_cols_updated, drop = FALSE) +
  labs(title = "Updated state proportions (5 original + pan-cancer relabeled)", x = NULL, y = "% of cells", fill = "State") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pie_df <- overall %>%
  mutate(label = ifelse(pct < 5, "", paste0(state, "\n", sprintf("%.1f%%", pct))))

p_pie <- ggplot(pie_df, aes(x = "", y = pct, fill = state)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2.5) +
  scale_fill_manual(values = group_cols_updated, drop = FALSE) +
  labs(title = "Updated state proportions: overall pie", fill = "State") +
  theme_void(base_size = 11)

ggsave(
  file.path(out_dir, "Auto_task4_unresolved_relabel_proportion.pdf"),
  p_bar + p_pie + plot_layout(widths = c(2, 1)),
  width = 18,
  height = 8
)

####################
# cell-cycle boxplot
####################
cc_df <- data.frame(
  state = factor(as.character(state_updated[names(cc_score)]), levels = state_level_order_updated),
  cc_score = as.numeric(cc_score),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(state))

p_cc <- ggplot(cc_df, aes(state, cc_score, fill = state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 0.15, alpha = 0.2) +
  scale_fill_manual(values = group_cols_updated, drop = FALSE) +
  labs(title = "Finalised states: Cell-cycle score by state", x = NULL, y = "Cell-cycle score") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")

ggsave(
  file.path(out_dir, "Auto_task4_unresolved_relabel_cc_boxplot.pdf"),
  p_cc,
  width = 12,
  height = 7
)

####################
# per-cell heatmap
####################
cc_in_ucell <- intersect(cc_mps, colnames(ucell_scores))
cc_raw <- as.matrix(ucell_scores[common_cells, cc_in_ucell, drop = FALSE])
mp_adj_cc <- z_normalise(cc_raw, sample_var, study_var)
mp_adj_all <- cbind(mp_adj_noncc[common_cells, , drop = FALSE], mp_adj_cc)

z_3ca_all <- z_normalise(
  as.matrix(ucell_3ca[common_cells, retained_3ca_available, drop = FALSE]),
  sample_var,
  study_var
)
raw_3ca_all <- as.matrix(ucell_3ca[common_cells, retained_3ca_available, drop = FALSE])

cna_status <- as.character(tmdata_all@meta.data[common_cells, "classification"])
names(cna_status) <- common_cells
cna_palette <- make_cna_palette(cna_status)

set.seed(42)
MAX_CELLS_TOTAL <- 8000
state_counts <- table(state_updated)
state_fracs <- state_counts / sum(state_counts)
cells_per_state <- pmax(round(state_fracs * MAX_CELLS_TOTAL), 20)
state_cells <- split(names(state_updated), state_updated)
cells_to_plot <- unlist(
  mapply(
    function(cells, n) sample(cells, min(length(cells), n)),
    state_cells,
    cells_per_state[names(state_cells)],
    SIMPLIFY = FALSE
  ),
  use.names = FALSE
)
if (length(cells_to_plot) > MAX_CELLS_TOTAL) {
  cells_to_plot <- sample(cells_to_plot, MAX_CELLS_TOTAL)
}

sub_scores_orig <- t(mp_adj_all[cells_to_plot, , drop = FALSE])
cc_block_order <- cc_mps[cc_mps %in% rownames(sub_scores_orig)]
non_cc_block_order <- mp_tree_order_names[
  mp_tree_order_names %in% rownames(sub_scores_orig) & !(mp_tree_order_names %in% cc_mps)
]
mp_row_order <- c(cc_block_order, non_cc_block_order)
sub_scores <- sub_scores_orig[mp_row_order, , drop = FALSE]

if (length(retained_3ca_available) > 0) {
  pan_mat <- t(z_3ca_all[cells_to_plot, , drop = FALSE])
  pan_mat <- pan_mat[retained_3ca_available[retained_3ca_available %in% rownames(pan_mat)], , drop = FALSE]
  rownames(pan_mat) <- clean_3ca_name(rownames(pan_mat))

  raw_pan_mat <- t(raw_3ca_all[cells_to_plot, , drop = FALSE])
  raw_pan_mat <- raw_pan_mat[retained_3ca_available[retained_3ca_available %in% rownames(raw_pan_mat)], , drop = FALSE]
  rownames(raw_pan_mat) <- paste0(clean_3ca_name(rownames(raw_pan_mat)), " (raw)")

  sub_scores <- rbind(sub_scores, pan_mat, raw_pan_mat)
}

row_ids_raw <- rownames(sub_scores)
mp_label_map <- mp_descriptions
missing_mps <- setdiff(rownames(sub_scores), names(mp_label_map))
if (length(missing_mps) > 0) mp_label_map[missing_mps] <- missing_mps
rownames(sub_scores) <- mp_label_map[rownames(sub_scores)]

present_states <- intersect(state_level_order_updated, unique(as.character(state_updated[cells_to_plot])))
if (length(present_states) == 0) {
  present_states <- unique(as.character(state_updated[cells_to_plot]))
}
split_vec <- factor(as.character(state_updated[cells_to_plot]), levels = present_states)

state_df_full <- data.frame(
  cell = names(state_updated),
  state = as.character(state_updated),
  sample = as.character(tmdata_all@meta.data[names(state_updated), "orig.ident"]),
  study = as.character(tmdata_all@meta.data[names(state_updated), "reference_batch"]),
  stringsAsFactors = FALSE
)
total_samples <- length(unique(state_df_full$sample))
total_studies <- length(unique(state_df_full$study))
state_div_df <- state_df_full %>%
  group_by(state) %>%
  summarise(
    sample_cov = n_distinct(sample) / max(total_samples, 1),
    study_cov = n_distinct(study) / max(total_studies, 1),
    diversity_score = 0.5 * sample_cov + 0.5 * study_cov,
    .groups = "drop"
  )
state_div_map <- setNames(state_div_df$diversity_score, state_div_df$state)
state_div_vals <- state_div_map[as.character(split_vec)]
state_div_vals[is.na(state_div_vals)] <- 0
names(state_div_vals) <- cells_to_plot

local_group_cols <- group_cols_updated[names(group_cols_updated) %in% present_states]
for (st in present_states) {
  if (!st %in% names(local_group_cols)) local_group_cols[st] <- "grey80"
}
local_group_cols <- local_group_cols[present_states]

study_levels <- unique(as.character(tmdata_all$reference_batch))
study_vals <- tmdata_all@meta.data[cells_to_plot, "reference_batch"]
study_cols <- setNames(
  DiscretePalette(length(study_levels), palette = "polychrome"),
  study_levels
)

max_cc <- max(cc_score[cells_to_plot], na.rm = TRUE)
if (!is.finite(max_cc) || max_cc <= 0) max_cc <- 1

col_ann <- HeatmapAnnotation(
  State = split_vec,
  CNA = cna_status[cells_to_plot],
  CC_score = cc_score[cells_to_plot],
  Diversity = state_div_vals,
  Study = study_vals,
  col = list(
    State = local_group_cols,
    CNA = cna_palette,
    CC_score = colorRamp2(c(0, max_cc), c("white", "darkgreen")),
    Diversity = colorRamp2(c(0, 1), c("grey95", "purple4")),
    Study = study_cols
  ),
  annotation_name_side = "left",
  show_legend = TRUE,
  na_col = "white",
  annotation_legend_param = list(
    State = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
    CNA = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
    CC_score = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
    Diversity = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
    Study = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8))
  )
)

mp_to_group <- rep("Other", length(row_ids_raw))
names(mp_to_group) <- row_ids_raw
mp_to_group[cc_mps[cc_mps %in% names(mp_to_group)]] <- "Cell_cycle"
for (grp in names(state_groups)) {
  grp_mps <- intersect(state_groups[[grp]], names(mp_to_group))
  mp_to_group[grp_mps] <- grp
}
pan_rows <- clean_3ca_name(retained_3ca_available)
pan_rows_raw <- paste0(pan_rows, " (raw)")
mp_to_group[intersect(names(mp_to_group), pan_rows)] <- "PanCancer"
mp_to_group[intersect(names(mp_to_group), pan_rows_raw)] <- "PanCancerRaw"

mp_group_names <- c(
  "Cell_cycle" = "Cell cycle",
  "Classic Proliferative" = "Classic\nProlif",
  "Basal to Intestinal Metaplasia" = "Basal-IM",
  "Stress-adaptive" = "Stress\nadaptive",
  "SMG-like Metaplasia" = "SMG-like\nMeta",
  "Immune Infiltrating" = "Immune\nInfiltrating",
  "PanCancer" = "Pan-Cancer\n(norm)",
  "PanCancerRaw" = "Pan-Cancer\n(raw)",
  "Other" = "Other"
)

group_colors_row <- c(
  setNames(group_cols[names(state_groups)], mp_group_names[names(state_groups)]),
  "Cell cycle" = "gold",
  "Pan-Cancer\n(norm)" = "#4B0082",
  "Pan-Cancer\n(raw)" = "#7F0000",
  "Other" = "grey70"
)

mp_group_labels <- mp_group_names[mp_to_group[row_ids_raw]]
names(mp_group_labels) <- mp_label_map[row_ids_raw]

adj_rows <- rownames(sub_scores)[!grepl("\\(raw\\)$", rownames(sub_scores))]
raw_rows <- rownames(sub_scores)[grepl("\\(raw\\)$", rownames(sub_scores))]

row_ann_adj <- rowAnnotation(
  MP_group = factor(mp_group_labels[adj_rows], levels = mp_group_names),
  col = list(MP_group = group_colors_row),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_legend_param = list(
    MP_group = list(
      title = "MP Group",
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    )
  )
)

if (length(raw_rows) > 0) {
  row_ann_raw <- rowAnnotation(
    MP_group = factor(mp_group_labels[raw_rows], levels = mp_group_names),
    col = list(MP_group = group_colors_row),
    show_annotation_name = FALSE,
    show_legend = FALSE
  )
}

lim <- as.numeric(quantile(abs(sub_scores), 0.98, na.rm = TRUE))
if (!is.finite(lim) || lim == 0) lim <- 1
col_fun_sc <- colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3"))

cc_labels <- mp_descriptions[cc_mps]
pan_labels <- clean_3ca_name(retained_3ca_available)
row_split_adj <- factor(
  ifelse(
    adj_rows %in% pan_labels,
    "Pan-Cancer\n(norm)",
    ifelse(adj_rows %in% cc_labels, "Cell cycle\nMPs", "State\nMPs")
  ),
  levels = c("Cell cycle\nMPs", "State\nMPs", "Pan-Cancer\n(norm)")
)

col_order <- (function() {
  col_order_list <- lapply(levels(split_vec), function(lvl) {
    idx <- which(as.character(split_vec) == lvl)
    if (length(idx) <= 1) return(idx)
    mat_lvl <- sub_scores[adj_rows, idx, drop = FALSE]
    dcols <- dist(t(mat_lvl))
    hc <- hclust(dcols, method = "ward.D2")
    idx[hc$order]
  })
  full_ord <- unlist(col_order_list, use.names = FALSE)
  if (length(full_ord) != ncol(sub_scores) || !setequal(full_ord, seq_len(ncol(sub_scores)))) {
    return(seq_len(ncol(sub_scores)))
  }
  full_ord
})()

ht_adj <- Heatmap(
  sub_scores[adj_rows, , drop = FALSE],
  name = "Adj score",
  col = col_fun_sc,
  top_annotation = col_ann,
  left_annotation = row_ann_adj,
  column_split = split_vec,
  row_split = row_split_adj,
  column_order = col_order,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_names = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8, fontface = "italic"),
  column_title_rot = 30,
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_gap = unit(2, "mm"),
  column_gap = unit(1.5, "mm"),
  use_raster = TRUE,
  raster_quality = 3,
  border = FALSE,
  rect_gp = gpar(col = NA),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

if (length(raw_rows) > 0) {
  raw_lim <- as.numeric(quantile(sub_scores[raw_rows, ], 0.98, na.rm = TRUE))
  if (is.na(raw_lim) || raw_lim == 0) raw_lim <- 0.15
  col_fun_raw <- colorRamp2(c(0, raw_lim), c("white", "firebrick3"))

  ht_raw <- Heatmap(
    sub_scores[raw_rows, , drop = FALSE],
    name = "Raw score",
    col = col_fun_raw,
    left_annotation = row_ann_raw,
    column_split = split_vec,
    column_order = col_order,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8, fontface = "italic"),
    column_gap = unit(1.5, "mm"),
    use_raster = TRUE,
    raster_quality = 3,
    show_heatmap_legend = TRUE,
    column_title = NULL,
    border = FALSE,
    rect_gp = gpar(col = NA),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    )
  )

  ht <- ht_adj %v% ht_raw
} else {
  ht <- ht_adj
}

pdf(file.path(out_dir, "Auto_task4_unresolved_relabel_heatmap.pdf"), width = 18, height = 10, useDingbats = FALSE)
draw(ht, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
grid.text(
  "Unresolved cells subclass: 3CA-based relabeling (noreg)",
  x = unit(5, "mm"),
  y = unit(1, "npc") - unit(5, "mm"),
  just = c("left", "top"),
  gp = gpar(fontsize = 14, fontface = "bold")
)
dev.off()

write.csv(
  data.frame(
    n_cells_total = length(common_cells),
    n_unresolved_input = length(unresolved_cells),
    n_relabelled = sum(state_B == "Unresolved" & state_updated != "Unresolved"),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "Auto_task4_unresolved_relabel_summary.csv"),
  row.names = FALSE
)
