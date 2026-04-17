suppressPackageStartupMessages({
  library("dplyr")
  library("tidyr")
  library("ggplot2")
  library("scales")
  library("patchwork")
  library("ggrepel")
  library("colorspace")
  library("Seurat")
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/analysis/summary/summary.R")
}

analysis_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(analysis_dir, "..", ".."))
out_dir <- file.path(project_dir, "sn_outs")

required_files <- c(
  file.path(out_dir, "filtering_summary_sn.csv"),
  file.path(out_dir, "filtered_sample_summary.csv"),
  file.path(out_dir, "expression_filter_reason_overall.csv"),
  file.path(out_dir, "snSeq_merged.rds")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required summary inputs: ", paste(basename(missing_files), collapse = ", "))
}

before_df <- read.csv(
  file.path(out_dir, "filtering_summary_sn.csv"),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

filtered_df <- read.csv(
  file.path(out_dir, "filtered_sample_summary.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)

reason_df <- read.csv(
  file.path(out_dir, "expression_filter_reason_overall.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)

merged_obj <- readRDS(file.path(out_dir, "snSeq_merged.rds"))

legacy_celltypes <- intersect(
  c("lymph", "keratinocyte", "erythrocyte", "neutrophil"),
  unique(as.character(merged_obj$celltype_update))
)
if (length(legacy_celltypes) > 0) {
  stop(
    "Unexpected legacy cell types present in snSeq_merged.rds: ",
    paste(sort(legacy_celltypes), collapse = ", ")
  )
}

final_cells_from_summary <- sum(filtered_df$cells_after, na.rm = TRUE)
if (final_cells_from_summary != ncol(merged_obj)) {
  stop(
    "Merged-object cell count does not match filtered_sample_summary.csv: ",
    final_cells_from_summary, " vs ", ncol(merged_obj)
  )
}

before_df <- before_df %>%
  rename(
    orig.ident = sample,
    Raw = raw,
    Mitochondrial_Filter = `mito_DNA\npercentage < 15`,
    Gene_Count_Filter = `number of\ngenes > 200`,
    Housekeeping_Filter = `housekeeping\nexpression > 0.5`
  )

# Recover ground truth metadata for ALL samples from the manifest
manifest_path <- file.path(out_dir, "names_tmdata_sn.txt")
if (!file.exists(manifest_path)) {
  stop("Missing names_tmdata_sn.txt in sn_outs. Required for metadata recovery.")
}
sample_paths <- readLines(manifest_path)
meta_map <- data.frame(
  orig.ident = basename(sample_paths),
  batch_map = dirname(sample_paths),
  stringsAsFactors = FALSE
) %>%
  mutate(
    tech_map = case_when(
      batch_map == "gemx" ~ "GEMX",
      batch_map == "cynthia_sn" ~ "snSeq",
      TRUE ~ "Multiome"
    )
  )

summary_df <- before_df %>%
  left_join(
    filtered_df %>%
      transmute(
        orig.ident = sample,
        technology,
        reference_batch,
        input_source,
        Final = cells_after,
        pct_retained,
        status
      ),
    by = "orig.ident"
  ) %>%
  left_join(meta_map, by = "orig.ident") %>%
  mutate(
    technology = coalesce(technology, tech_map),
    reference_batch = coalesce(reference_batch, batch_map),
    input_source = coalesce(input_source, batch_map),
    Final = replace_na(Final, 0L),
    pct_retained = replace_na(pct_retained, 0),
    status = replace_na(status, "no_cell")
  ) %>%
  select(-tech_map, -batch_map)

tech_colors <- c(
  "snSeq" = "steelblue",
  "GEMX" = "darkorange",
  "Multiome" = "forestgreen"
)

####################
# scRef publication-matched celltype colours
####################
celltype_palette_base <- c(
  "b.cell" = "#E41A1C",
  "dendritic" = "#377EB8",
  "endothelial" = "#4DAF4A",
  "epithelial" = "#984EA3",
  "fibroblast" = "#FF7F00",
  "macrophage" = "#A65628",
  "mast" = "#F781BF",
  "nk.cell" = "#FFD700",
  "plasma" = "#999999",
  "t.cell" = "#00CED1"
)

bar_data <- summary_df %>%
  summarize(
    Raw = sum(Raw, na.rm = TRUE),
    `Mito_DNA < 15` = sum(Mitochondrial_Filter, na.rm = TRUE),
    `Min_genes > 200` = sum(Gene_Count_Filter, na.rm = TRUE),
    `HK_expr > 0.5` = sum(Housekeeping_Filter, na.rm = TRUE),
    `Final (Expr-filtered singlets)` = sum(Final, na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "filter_step", values_to = "cell_count") %>%
  mutate(
    filter_step = factor(
      filter_step,
      levels = c(
        "Raw",
        "Mito_DNA < 15",
        "Min_genes > 200",
        "HK_expr > 0.5",
        "Final (Expr-filtered singlets)"
      )
    )
  )

bar_plot <- ggplot(bar_data, aes(x = filter_step, y = cell_count, fill = filter_step)) +
  geom_col(width = 0.8, color = "gray20", linewidth = 0.25) +
  geom_text_repel(
    aes(label = comma(cell_count), y = cell_count),
    nudge_y = 0.04 * max(bar_data$cell_count, na.rm = TRUE),
    direction = "y",
    box.padding = 0.25,
    point.padding = 0.25,
    segment.color = NA,
    min.segment.length = 0,
    size = 3.2,
    max.overlaps = Inf
  ) +
  scale_fill_brewer(palette = "Blues", direction = -1, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14)), labels = comma) +
  labs(title = "Total Cells Remaining After Each Filter", x = NULL, y = "Number of Cells") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(5, 5, 5, 5)
  )

tech_pie_data <- summary_df %>%
  group_by(technology) %>%
  summarise(Final = sum(Final, na.rm = TRUE), .groups = "drop") %>%
  filter(Final > 0) %>%
  mutate(
    fraction = Final / sum(Final),
    label = paste0(technology, "\n", comma(Final), " (", percent(fraction, accuracy = 0.1), ")")
  ) %>%
  arrange(desc(technology)) %>%
  mutate(ypos = cumsum(fraction) - 0.5 * fraction)

tech_pie <- ggplot(tech_pie_data, aes(x = 1, y = fraction, fill = technology)) +
  geom_col(color = "white", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 1.95) +
  geom_segment(
    aes(x = 1.02, xend = 1.15, y = ypos, yend = ypos),
    inherit.aes = FALSE,
    color = "grey50",
    linewidth = 0.3
  ) +
  geom_label_repel(
    data = tech_pie_data,
    aes(x = 1.55, y = ypos, label = label),
    inherit.aes = FALSE,
    direction = "y",
    nudge_x = 0.1,
    box.padding = 0.45,
    point.padding = 0,
    segment.color = "grey50",
    min.segment.length = 0,
    size = 3.3,
    max.overlaps = Inf,
    label.r = grid::unit(0.1, "lines")
  ) +
  scale_fill_manual(values = tech_colors, guide = "none") +
  labs(title = "Final Distribution by Technology") +
  theme_void(base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.margin = margin(5, 12, 5, 5)
  )

pair1_combo <- (bar_plot + tech_pie) + plot_layout(ncol = 2, widths = c(1.6, 1))
ggsave(file.path(analysis_dir, "pair1_filter_overview.png"), bar_plot, width = 8, height = 4, dpi = 300)
ggsave(file.path(analysis_dir, "pair1_combo.png"), pair1_combo, width = 12, height = 6, dpi = 300)

reason_order <- c(
  "keep",
  "doublet",
  "marker_inconsistent",
  "marker_positive_unresolved",
  "poor_assigned",
  "poor_unresolved",
  "too_few_cells",
  "other"
)

reason_labels <- c(
  keep = "Kept",
  doublet = "Doublet",
  marker_inconsistent = "Marker mismatch",
  marker_positive_unresolved = "Marker-positive unresolved",
  poor_assigned = "Poor assigned marker signal",
  poor_unresolved = "Poor unresolved signal",
  too_few_cells = "Too few cells",
  other = "Other"
)

reason_palette <- c(
  keep = "#1b9e77",
  doublet = "#d95f02",
  marker_inconsistent = "#7570b3",
  marker_positive_unresolved = "#e6ab02",
  poor_assigned = "#666666",
  poor_unresolved = "#a6761d",
  too_few_cells = "#e7298a",
  other = "#bdbdbd"
)

reason_df <- reason_df %>%
  mutate(
    filter_reason = factor(filter_reason, levels = reason_order),
    label = reason_labels[as.character(filter_reason)]
  ) %>%
  arrange(filter_reason)

reason_plot <- ggplot(reason_df, aes(x = filter_reason, y = n_cells, fill = filter_reason)) +
  geom_col(width = 0.8, color = "gray20", linewidth = 0.25) +
  geom_text(
    aes(label = comma(n_cells)),
    vjust = -0.35,
    size = 3.2
  ) +
  scale_fill_manual(values = reason_palette, guide = "none") +
  scale_x_discrete(labels = reason_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = comma) +
  labs(title = "Why Cells Were Filtered", x = NULL, y = "Number of Cells") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(file.path(analysis_dir, "filter_reason_breakdown.png"), reason_plot, width = 9, height = 5, dpi = 300)

####################
# Pair-1 by batch (scRef per-study style, adapted to reference_batch)
####################
batch_summary <- summary_df %>%
  group_by(reference_batch) %>%
  summarise(
    Raw = sum(Raw, na.rm = TRUE),
    Mito_Filter = sum(Mitochondrial_Filter, na.rm = TRUE),
    Gene_Filter = sum(Gene_Count_Filter, na.rm = TRUE),
    HK_Filter = sum(Housekeeping_Filter, na.rm = TRUE),
    Final = sum(Final, na.rm = TRUE),
    .groups = "drop"
  )

batch_long <- batch_summary %>%
  pivot_longer(
    cols = c(Raw, Mito_Filter, Gene_Filter, HK_Filter, Final),
    names_to = "step",
    values_to = "cell_count"
  ) %>%
  mutate(
    step_order = case_when(
      step == "Raw" ~ 1L,
      step == "Mito_Filter" ~ 2L,
      step == "Gene_Filter" ~ 3L,
      step == "HK_Filter" ~ 4L,
      step == "Final" ~ 5L
    ),
    step_label = case_when(
      step == "Raw" ~ "Raw",
      step == "Mito_Filter" ~ "Mito_DNA < 15",
      step == "Gene_Filter" ~ "Min_genes > 200",
      step == "HK_Filter" ~ "HK_expr > 0.5",
      step == "Final" ~ "Final singlets",
      TRUE ~ step
    ),
    step_f = factor(
      step_order,
      levels = c(1, 2, 3, 4, 5),
      labels = c("Raw", "Mito", "Gene", "HK", "Singlets")
    )
  )

pair1_by_batch <- ggplot(batch_long, aes(x = step_f, y = cell_count, fill = step_f)) +
  geom_col(width = 0.65, color = "gray20", linewidth = 0.25) +
  geom_text(
    aes(label = comma(cell_count), y = cell_count),
    vjust = -0.4,
    size = 3
  ) +
  facet_wrap(~ reference_batch, scales = "free_y") +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0.12, 0.14))
  ) +
  scale_fill_brewer(
    palette = "Blues",
    direction = -1,
    guide = "none"
  ) +
  labs(
    title = "Cells Remaining After Each Filter (by reference batch)",
    x = NULL,
    y = "Number of Cells"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.margin = margin(10, 10, 20, 10),
    panel.spacing = unit(1, "lines")
  ) +
  coord_cartesian(clip = "off") +
  geom_text(
    data = batch_long,
    aes(x = step_f, y = 0, label = step_label),
    inherit.aes = FALSE,
    vjust = 1.6,
    size = 2.7
  )

ggsave(file.path(analysis_dir, "pair1_by_batch.png"), pair1_by_batch, width = 16, height = 5, dpi = 300)

celltype_df <- as.data.frame(table(merged_obj$celltype_update), stringsAsFactors = FALSE) %>%
  rename(celltype_update = Var1, n = Freq) %>%
  filter(!is.na(celltype_update), nzchar(celltype_update), n > 0) %>%
  arrange(desc(n)) %>%
  mutate(
    percent = 100 * n / sum(n),
    display_label = paste0("N = ", comma(n), "\n(", sprintf("%.1f", percent), "%)")
  )

if (sum(celltype_df$n) != ncol(merged_obj)) {
  stop("Cell type counts do not sum to the merged-object cell count.")
}

extra_celltypes <- setdiff(as.character(celltype_df$celltype_update), names(celltype_palette_base))
celltype_palette <- celltype_palette_base
if (length(extra_celltypes) > 0) {
  celltype_palette <- c(
    celltype_palette,
    setNames(qualitative_hcl(length(extra_celltypes), palette = "Dark 3"), extra_celltypes)
  )
}

celltype_df <- celltype_df %>%
  filter(celltype_update %in% names(celltype_palette)) %>%
  arrange(n) %>%
  mutate(celltype_update = factor(celltype_update, levels = celltype_update))

bar_pad <- max(celltype_df$n) * 0.02

celltype_bar <- ggplot(celltype_df, aes(x = celltype_update, y = -n, fill = celltype_update)) +
  geom_col(width = 0.8, color = "black") +
  geom_text(
    aes(label = display_label),
    hjust = 1,
    nudge_y = -bar_pad * 1.8,
    lineheight = 0.95,
    size = 4
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = celltype_palette, drop = FALSE) +
  scale_y_continuous(
    labels = function(x) comma(abs(x)),
    expand = expansion(mult = c(0.20, 0.02))
  ) +
  labs(title = "Cell Type Counts", x = NULL, y = NULL) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 12, face = "bold", colour = "black"),
    axis.text.x = element_text(size = 10, colour = "black"),
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  )

ct_pie_data <- celltype_df %>%
  mutate(
    fraction = n / sum(n),
    label = paste0(as.character(celltype_update), "\n", comma(n), " (", percent(fraction, accuracy = 0.1), ")"),
    pie_label = ifelse(fraction >= 0.03, label, "")
  ) %>%
  arrange(desc(celltype_update)) %>%
  mutate(ypos = cumsum(fraction) - 0.5 * fraction)

celltype_pie <- ggplot(ct_pie_data, aes(x = 1, y = fraction, fill = celltype_update)) +
  geom_col(color = "white", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.1) +
  geom_segment(
    aes(x = 1.02, xend = 1.12, y = ypos, yend = ypos),
    inherit.aes = FALSE,
    color = "grey70",
    linewidth = 0.25
  ) +
  geom_label_repel(
    data = dplyr::filter(ct_pie_data, pie_label != ""),
    aes(x = 1.7, y = ypos, label = pie_label),
    inherit.aes = FALSE,
    direction = "y",
    nudge_x = 0.15,
    box.padding = 0.55,
    point.padding = 0,
    segment.color = "grey60",
    min.segment.length = 0,
    size = 3.1,
    max.overlaps = Inf,
    label.r = grid::unit(0.1, "lines")
  ) +
  scale_fill_manual(values = celltype_palette, guide = "none") +
  labs(title = "Cell Type Composition") +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.margin = margin(5, 12, 5, 5)
  )

pair2_combo <- (celltype_bar + celltype_pie) + plot_layout(ncol = 2, widths = c(1.7, 1.1))
ggsave(file.path(analysis_dir, "pair2_celltype_composition.png"), pair2_combo, width = 13, height = 6, dpi = 300)

####################
# Pair-3 UMAPs (publication-style final celltype + technology)
####################
if (!"umap" %in% names(merged_obj@reductions)) {
  stop("Missing UMAP reduction in snSeq_merged.rds. Rerun step 4 to regenerate the merged object.")
}

umap_mat <- Embeddings(merged_obj, "umap")
umap_df <- data.frame(
  UMAP_1 = umap_mat[, 1],
  UMAP_2 = umap_mat[, 2],
  celltype_update = as.character(merged_obj$celltype_update),
  technology = as.character(merged_obj$technology),
  stringsAsFactors = FALSE
)

set.seed(0)
umap_df <- umap_df[sample(nrow(umap_df)), , drop = FALSE]

umap_celltype_levels <- c(
  intersect(names(celltype_palette_base), unique(umap_df$celltype_update)),
  sort(setdiff(unique(umap_df$celltype_update), names(celltype_palette_base)))
)

umap_df <- umap_df %>%
  mutate(
    celltype_update = factor(celltype_update, levels = umap_celltype_levels),
    technology = factor(technology, levels = names(tech_colors))
  )

umap_celltype <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = celltype_update)) +
  geom_point(size = 0.03, alpha = 0.8) +
  scale_color_manual(values = celltype_palette, drop = FALSE) +
  labs(
    title = "Final Cell Type",
    subtitle = paste0("Total filtered cells: ", comma(nrow(umap_df))),
    color = "Celltype",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right",
    legend.text = element_text(size = 9),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    aspect.ratio = 1
  )

umap_technology <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = technology)) +
  geom_point(size = 0.03, alpha = 0.8) +
  scale_color_manual(values = tech_colors, drop = FALSE) +
  labs(
    title = "Technology",
    subtitle = "Filtered cohort",
    color = "Technology",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right",
    legend.text = element_text(size = 9),
    aspect.ratio = 1
  )

pair3_umap <- (umap_celltype + umap_technology) + plot_layout(ncol = 2, widths = c(1.1, 1))
ggsave(file.path(analysis_dir, "pair3_umap_celltype_technology.png"), pair3_umap, width = 13, height = 6, dpi = 300)

counts <- summary_df %>%
  transmute(
    orig.ident,
    before = Raw,
    after = Final,
    technology
  )

counts_long <- counts %>%
  pivot_longer(
    cols = c("before", "after"),
    names_to = "filter_status",
    values_to = "cell_count"
  ) %>%
  mutate(
    filter_status = factor(filter_status, levels = c("before", "after")),
    orig.ident = factor(orig.ident)
  ) %>%
  mutate(
    orig.ident = factor(
      orig.ident,
      levels = counts %>%
        arrange(factor(technology, levels = names(tech_colors)), orig.ident) %>%
        pull(orig.ident)
    )
  )

axis_label_colors <- ifelse(
  levels(counts_long$orig.ident) %in% counts$orig.ident[counts$after == 0],
  "red",
  "black"
)

sample_plot <- ggplot(
  counts_long,
  aes(x = orig.ident, y = cell_count + 1, fill = technology, alpha = filter_status)
) +
  geom_col(
    width = 0.75,
    position = position_dodge(width = 0.85),
    color = "gray20",
    linewidth = 0.15
  ) +
  scale_y_log10(
    breaks = c(1, 11, 101, 1001, 10001, 100001),
    labels = c("0", "10", "100", "1k", "10k", "100k")
  ) +
  scale_fill_manual(values = tech_colors, drop = FALSE) +
  scale_alpha_manual(
    name = "",
    values = c("before" = 1.0, "after" = 0.6),
    labels = c("Raw", "Final")
  ) +
  labs(
    x = NULL,
    y = "Number of Cells (log10 scale)",
    fill = "Technology"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = axis_label_colors),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, face = "italic"),
    legend.position = "top",
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed")
  )

ggsave(file.path(analysis_dir, "sample_retention.png"), sample_plot, width = 16, height = 7, dpi = 300)

retention_df <- summary_df %>%
  transmute(
    orig.ident,
    technology,
    pct_retained,
    status
  ) %>%
  mutate(
    orig.ident = factor(
      orig.ident,
      levels = counts %>%
        arrange(factor(technology, levels = names(tech_colors)), orig.ident) %>%
        pull(orig.ident)
    )
  )

retention_plot <- ggplot(retention_df, aes(x = orig.ident, y = pct_retained, fill = technology, alpha = status)) +
  geom_col(color = "gray20", linewidth = 0.2) +
  scale_fill_manual(values = tech_colors, drop = FALSE) +
  scale_alpha_manual(values = c("ok" = 1, "no_cell" = 0.55), guide = "none") +
  labs(
    title = "Step-4 Retention Rate By Sample",
    x = NULL,
    y = "Cells retained after expression filtering (%)",
    fill = "Technology"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = axis_label_colors),
    legend.position = "top"
  )

ggsave(file.path(analysis_dir, "sample_retention_pct.png"), retention_plot, width = 16, height = 6, dpi = 300)

pdf(file.path(analysis_dir, "summary_overview.pdf"), width = 13, height = 8)
print(pair1_combo)
print(pair1_by_batch)
print(reason_plot)
print(pair2_combo)
print(pair3_umap)
print(sample_plot)
print(retention_plot)
dev.off()
