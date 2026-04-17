suppressPackageStartupMessages({
  library("Seurat")
  library("dplyr")
  library("purrr")
  library("tibble")
  library("readxl")
  library("ggplot2")
  library("patchwork")
  library("Matrix")
  library("grid")
})

set.seed(1)

out_dir <- "sn_outs"
manifest_path <- "sn_outs/sample_manifest.csv"
marker_workbook <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Marker_Genes.xlsx"
annotation_marker_path <- "sn_outs/annotation_score_markers.csv"
expression_filter_marker_path <- "sn_outs/expression_filter_markers.csv"

required_files <- c(
  manifest_path,
  marker_workbook,
  annotation_marker_path,
  expression_filter_marker_path
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required annotation plotting inputs: ", paste(basename(missing_files), collapse = ", "))
}

sample_manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
sample_ids <- unique(sample_manifest$sample)

## ---- Parameters ----
## Available choices for sample_mode:
## 1. "all"    - Process all samples in the manifest.
## 2. "subset" - Process a random subset (defined by n_samples and sample_seed).
## 3. "<ID>"   - Enter a specific sample ID (e.g., "patient_A_pre") to plot only that sample.
sample_mode   <- "E_post_T1_biopsy" 

## If sample_mode is "subset", how many samples should be randomly picked?
n_samples     <- 10       

## Random seed for reproducibility when using "subset" mode.
sample_seed   <- 1        

## ---- Process Selection ----

if (sample_mode == "all") {
  selected_samples <- sample_ids
  selection_label  <- "all_samples"
} else if (sample_mode == "subset") {
  set.seed(sample_seed)
  selected_samples <- sort(sample(sample_ids, size = min(n_samples, length(sample_ids)), replace = FALSE))
  selection_label  <- paste0("subset_", length(selected_samples), "_seed_", sample_seed)
} else {
  # Direct sample selection
  if (!sample_mode %in% sample_ids) {
    stop("Sample mode '", sample_mode, "' is not 'all', 'subset', or a valid sample ID from the manifest.")
  }
  selected_samples <- sample_mode
  selection_label  <- paste0("single_", sample_mode)
}

pdf_path      <- file.path(out_dir, paste0("snseq_annotation_featureplots_", selection_label, ".pdf"))
status_path   <- file.path(out_dir, paste0("annotation_featureplots_status_", selection_label, ".csv"))
selected_path <- file.path(out_dir, paste0("selected_samples_", selection_label, ".csv"))

## ---- Cell type definitions ----

celltype_map <- c(
  "Fibroblast" = "fibroblast", "Macrophage" = "macrophage", "Mast" = "mast",
  "B cell" = "b.cell", "T cell" = "t.cell", "Dendritic" = "dendritic",
  "Endothelial" = "endothelial", "Epithelial" = "epithelial",
  "NK cell" = "nk.cell", "Plasma" = "plasma"
)

celltype_order <- c(
  "epithelial", "endothelial", "fibroblast", "macrophage",
  "b.cell", "plasma", "t.cell", "nk.cell", "mast", "dendritic"
)

celltype_display <- c(
  epithelial = "Epithelial", endothelial = "Endothelial",
  fibroblast = "Fibroblast", macrophage = "Macrophage",
  "b.cell" = "B cell", plasma = "Plasma",
  "t.cell" = "T cell", "nk.cell" = "NK cell",
  mast = "Mast", dendritic = "Dendritic"
)

celltype_colors <- c(
  epithelial = "#D73027", endothelial = "#2C7BB6",
  fibroblast = "#A6761D", macrophage = "#1B9E77",
  "b.cell" = "#7570B3", plasma = "#E7298A",
  "t.cell" = "#66A61E", "nk.cell" = "#1F78B4",
  mast = "#E6AB02", dendritic = "#6A3D9A"
)

canonical_markers <- list(
  fibroblast = c("COL3A1", "COL1A2", "LUM", "COL1A1", "COL6A3", "DCN"),
  macrophage = c("CSF1R", "TYROBP", "CD14", "CD163", "AIF1", "CD68"),
  mast = c("MS4A2", "CPA3", "TPSB2", "TPSAB1"),
  epithelial = c("KRT7", "MUC1", "KRT19", "EPCAM"),
  t.cell = c("CD3E", "CD3D", "CD2", "CD3G"),
  b.cell = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1"),
  nk.cell = c("GNLY", "NKG7", "PRF1", "GZMB", "KLRB1"),
  plasma = c("MZB1", "JCHAIN", "DERL3"),
  dendritic = c("CLEC10A", "CCR7", "CD86"),
  endothelial = c("ENG", "CLEC14A", "CLDN5", "VWF", "CDH5")
)

## Filter status colors
filter_status_colors <- c(
  "Kept" = "#2CA02C",
  "Low Expr" = "#FF7F0E",
  "Others" = "#D62728"
)

## ---- Marker loading ----

combine_marker_scores <- function(df, w_specificity = 0.2, w_sensitivity = 0.8) {
  pr <- function(x) {
    r <- rank(x, ties.method = "average", na.last = "keep")
    r / (sum(!is.na(x)) + 1)
  }
  combined <- (w_specificity * pr(df$specificity) + w_sensitivity * pr(df$sensitivity)) /
    (w_specificity + w_sensitivity)
  df %>% mutate(Combined = combined) %>% arrange(desc(Combined))
}

markers_3ca <- read_excel(marker_workbook, sheet = 1) %>%
  filter(specificity > 0.2, cell_type != "Malignant") %>%
  mutate(cell_type = recode(cell_type, !!!celltype_map))

markers_ranked <- markers_3ca %>%
  split(.$cell_type) %>%
  lapply(combine_marker_scores)

initial_top50 <- imap(markers_ranked, ~ {
  .x %>% arrange(desc(Combined)) %>% slice_head(n = 50) %>% pull(gene)
})
initial_top50 <- initial_top50[celltype_order]

annotation_marker_tbl <- read.csv(annotation_marker_path, stringsAsFactors = FALSE)
expression_filter_marker_tbl <- read.csv(expression_filter_marker_path, stringsAsFactors = FALSE) %>%
  rename(celltype = celltype) %>%
  mutate(score_weight = ifelse(is.na(score_weight), 0.1, pmax(score_weight, 0.1)))

## ---- Helper functions ----

get_expr_layer <- function(obj, layer = "data") {
  layer_data <- tryCatch(obj@assays$RNA@layers[[layer]], error = function(e) NULL)
  if (!is.null(layer_data)) {
    if (length(rownames(layer_data)) == 0) rownames(layer_data) <- rownames(obj)
    if (length(colnames(layer_data)) == 0) colnames(layer_data) <- colnames(obj)
    return(layer_data)
  }
  GetAssayData(obj, slot = layer)
}

row_zscore <- function(mat) {
  out <- t(scale(t(mat)))
  out[!is.finite(out)] <- 0
  out
}

add_score_matrix <- function(obj, score_mat, prefix) {
  for (ct in celltype_order) {
    obj[[paste0(prefix, ct)]] <- score_mat[colnames(obj), ct]
  }
  obj
}

## ---- Score builders (all produce per-cell scores, not fractions) ----

## Page 1: Initial 3CA top-50 Z-scores (same as Clustering.R)
build_initial_score_matrix <- function(obj, marker_sets) {
  expr <- get_expr_layer(obj, "data")
  marker_sets <- marker_sets[celltype_order]
  gene_union <- intersect(unique(unlist(marker_sets, use.names = FALSE)), rownames(expr))
  score_mat <- matrix(0, nrow = ncol(obj), ncol = length(celltype_order),
                      dimnames = list(colnames(obj), celltype_order))
  if (length(gene_union) == 0) return(score_mat)
  expr_mat <- as.matrix(expr[gene_union, , drop = FALSE])
  expr_z <- row_zscore(expr_mat)
  for (ct in celltype_order) {
    genes <- intersect(marker_sets[[ct]], rownames(expr_z))
    if (length(genes) == 0) next
    score_mat[, ct] <- Matrix::colMeans(expr_z[genes, , drop = FALSE])
  }
  row_zscore(score_mat)
}

## Page 2: Refined weighted Z-scores (same as Annotation.R / Expr_filtering.R)
build_refined_score_matrix <- function(obj, marker_tbl) {
  marker_tbl <- marker_tbl %>%
    filter(celltype %in% celltype_order) %>%
    mutate(gene = as.character(gene))
  score_mat <- matrix(0, nrow = ncol(obj), ncol = length(celltype_order),
                      dimnames = list(colnames(obj), celltype_order))
  gene_union <- intersect(unique(marker_tbl$gene), rownames(obj))
  if (length(gene_union) == 0) return(score_mat)
  expr_mat <- as.matrix(get_expr_layer(obj, "data")[gene_union, , drop = FALSE])
  expr_z <- row_zscore(expr_mat)
  for (ct in celltype_order) {
    sub_tbl <- marker_tbl %>% filter(celltype == ct, gene %in% gene_union)
    if (nrow(sub_tbl) == 0) next
    weights <- sub_tbl$score_weight
    names(weights) <- sub_tbl$gene
    weights <- weights / sum(weights)
    score_mat[, ct] <- as.numeric(crossprod(weights, expr_z[sub_tbl$gene, , drop = FALSE]))
  }
  row_zscore(score_mat)
}

## Page 3: Canonical marker mean expression (log-normalized data layer)
build_canonical_score_matrix <- function(obj) {
  expr <- get_expr_layer(obj, "data")
  score_mat <- matrix(0, nrow = ncol(obj), ncol = length(celltype_order),
                      dimnames = list(colnames(obj), celltype_order))
  for (ct in celltype_order) {
    genes <- intersect(canonical_markers[[ct]], rownames(expr))
    if (length(genes) == 0) next
    score_mat[, ct] <- as.numeric(Matrix::colMeans(expr[genes, , drop = FALSE]))
  }
  score_mat
}

## Page 4: Expression-filter panel mean expression (top-4 approach from Expr_filtering.R)
build_filter_score_matrix <- function(obj, marker_tbl, expr_cut = 0.5) {
  expr <- get_expr_layer(obj, "data")
  score_mat <- matrix(0, nrow = ncol(obj), ncol = length(celltype_order),
                      dimnames = list(colnames(obj), celltype_order))

  for (ct in celltype_order) {
    genes <- marker_tbl %>% filter(celltype == ct) %>% pull(gene)
    genes <- intersect(genes, rownames(expr))
    if (length(genes) == 0) next

    ## Same logic as Expr_filtering.R: mean of the top up-to-4 expressed markers above threshold
    sub_expr <- as.matrix(expr[genes, , drop = FALSE])
    scores <- apply(sub_expr, 2, function(cell_vals) {
      above <- cell_vals[cell_vals > expr_cut]
      if (length(above) == 0) return(0)
      mean(sort(above, decreasing = TRUE)[1:min(4, length(above))])
    })
    score_mat[, ct] <- scores
  }
  score_mat
}

## ---- Fixed scale limits per page type ----
page_limits <- list(
  initial_score_ = c(-2.5, 2.5),
  refined_score_ = c(-2.5, 2.5),
  canonical_score_ = NULL,   # auto from data
  filter_score_ = NULL       # auto from data
)

pretty_sample_title <- function(sample, n_cells, page_label) {
  paste0(sample, " (n = ", format(n_cells, big.mark = ","), ") | ", page_label)
}

## ---- Adaptive point size ----
get_point_size <- function(n_cells) {
  if (n_cells > 40000) 0.35
  else if (n_cells > 20000) 0.8
  else if (n_cells > 10000) 0.9
  else if (n_cells > 5000) 1.0
  else if (n_cells > 2000) 1.2
  else 1.4
}

## ---- Plot functions ----

base_dimplot <- function(obj, group_by, title_text, label = FALSE,
                         palette = NULL, show_legend = FALSE) {
  point_size <- get_point_size(ncol(obj))
  p <- DimPlot(object = obj, group.by = group_by, label = label, repel = label,
               raster = FALSE, pt.size = point_size,
               label.size = 4, combine = FALSE)[[1]]
  if (!is.null(palette)) {
    present_levels <- unique(as.character(obj@meta.data[[group_by]]))
    present_levels <- present_levels[!is.na(present_levels)]
    p <- p + scale_color_manual(
      values = palette,
      breaks = intersect(names(palette), present_levels),
      labels = if (!is.null(names(palette))) {
        vals <- intersect(names(palette), present_levels)
        if (all(vals %in% names(celltype_display))) celltype_display[vals] else vals
      },
      drop = TRUE
    )
  }
  p + ggtitle(title_text) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width = unit(0.35, "cm")
    )
}

## Diverging blue-white-red for Z-scores (has negative values)
score_featureplot_diverging <- function(obj, feature_name, title_text, limits) {
  point_size <- get_point_size(ncol(obj))
  coords <- as.data.frame(Embeddings(obj, reduction = "umap")[, 1:2, drop = FALSE])
  colnames(coords) <- c("UMAP_1", "UMAP_2")
  values <- as.numeric(obj@meta.data[[feature_name]])
  plot_df <- cbind(coords, score = values)
  plot_df$score_plot <- pmin(pmax(plot_df$score, limits[1]), limits[2])

  ## Sort by absolute value so extreme values plot on top
  plot_df <- plot_df[order(abs(plot_df$score_plot)), , drop = FALSE]

  ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(color = "grey88", size = point_size + 0.1, stroke = 0) +
    geom_point(aes(color = score_plot), size = point_size, stroke = 0, alpha = 0.95) +
    scale_color_gradient2(
      low = "#2166AC", mid = "white", high = "#D73027", midpoint = 0,
      limits = limits, oob = scales::squish, name = NULL,
      guide = guide_colorbar(barwidth = unit(0.35, "cm"), barheight = unit(1.8, "cm"))
    ) +
    ggtitle(title_text) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      axis.text = element_blank(), axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 6),
      legend.margin = margin(0, 0, 0, -4)
    )
}

## Sequential grey-to-red for expression scores (non-negative values)
score_featureplot_sequential <- function(obj, feature_name, title_text, limits) {
  point_size <- get_point_size(ncol(obj))
  coords <- as.data.frame(Embeddings(obj, reduction = "umap")[, 1:2, drop = FALSE])
  colnames(coords) <- c("UMAP_1", "UMAP_2")
  values <- as.numeric(obj@meta.data[[feature_name]])
  plot_df <- cbind(coords, score = values)
  plot_df$score_plot <- pmin(pmax(plot_df$score, limits[1]), limits[2])
  plot_df <- plot_df[order(plot_df$score_plot), , drop = FALSE]

  ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(color = "grey88", size = point_size + 0.1, stroke = 0) +
    geom_point(aes(color = score_plot), size = point_size, stroke = 0, alpha = 0.95) +
    scale_color_gradient(
      low = "grey88", high = "#D73027",
      limits = limits, oob = scales::squish, name = NULL,
      guide = guide_colorbar(barwidth = unit(0.35, "cm"), barheight = unit(1.8, "cm"))
    ) +
    ggtitle(title_text) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      axis.text = element_blank(), axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 6),
      legend.margin = margin(0, 0, 0, -4)
    )
}

## ---- Page builders ----

## Z-score pages (initial and refined): diverging blue-white-red
build_zscore_page <- function(obj, score_prefix, page_title,
                              label_column, label_title) {
  limits <- page_limits[[score_prefix]]
  plots <- list(base_dimplot(obj, "seurat_clusters", "Clusters", label = TRUE))
  for (ct in celltype_order) {
    plots[[length(plots) + 1]] <- score_featureplot_diverging(
      obj, paste0(score_prefix, ct), celltype_display[[ct]], limits = limits
    )
  }
  plots[[length(plots) + 1]] <- base_dimplot(
    obj, label_column, label_title,
    label = FALSE, palette = celltype_colors, show_legend = TRUE
  )
  wrap_plots(plots, ncol = 4, nrow = 3) +
    plot_annotation(title = page_title)
}

## Expression score pages (canonical and filter): sequential grey-to-red
## Computes a shared scale across all celltypes for comparability
build_expression_page <- function(obj, score_prefix, page_title,
                                  label_column, label_title,
                                  label_palette = celltype_colors) {
  ## Compute shared scale across all celltypes
  all_vals <- unlist(lapply(celltype_order, function(ct) {
    obj@meta.data[[paste0(score_prefix, ct)]]
  }))
  shared_max <- max(all_vals, na.rm = TRUE)
  shared_max <- max(shared_max, 0.1)  # ensure non-zero
  limits <- c(0, ceiling(shared_max * 10) / 10)  # round up to nearest 0.1

  plots <- list(base_dimplot(obj, "seurat_clusters", "Clusters", label = TRUE))
  for (ct in celltype_order) {
    plots[[length(plots) + 1]] <- score_featureplot_sequential(
      obj, paste0(score_prefix, ct), celltype_display[[ct]], limits = limits
    )
  }
  plots[[length(plots) + 1]] <- base_dimplot(
    obj, label_column, label_title,
    label = FALSE, palette = label_palette, show_legend = TRUE
  )
  wrap_plots(plots, ncol = 4, nrow = 3) +
    plot_annotation(title = page_title)
}

## ---- Filter status helper ----

add_filter_status <- function(obj) {
  if (!"filter_reason" %in% colnames(obj@meta.data)) {
    obj$filter_status <- "No Filter Data"
    return(obj)
  }
  obj$filter_status <- case_when(
    obj$filter_reason == "keep" ~ "Kept",
    obj$filter_reason %in% c("poor_assigned", "poor_unresolved",
                             "marker_positive_unresolved") ~ "Low Expr",
    TRUE ~ "Others"
  )
  obj$filter_status <- factor(obj$filter_status, levels = c("Kept", "Low Expr", "Others"))
  obj
}

## ---- Main loop ----

write.csv(
  data.frame(sample = selected_samples, selection_mode = sample_mode,
             sample_seed = sample_seed, stringsAsFactors = FALSE),
  selected_path, row.names = FALSE
)

plot_status <- vector("list", length(selected_samples))
names(plot_status) <- selected_samples

grDevices::cairo_pdf(pdf_path, width = 18, height = 13, onefile = TRUE)

for (idx in seq_along(selected_samples)) {
  sample <- selected_samples[idx]
  anno_path <- file.path(out_dir, "by_samples", sample, paste0(sample, "_anno.rds"))

  if (!file.exists(anno_path)) {
    plot_status[[sample]] <- data.frame(
      sample = sample, status = "missing_anno_rds", stringsAsFactors = FALSE
    )
    next
  }

  obj <- readRDS(anno_path)
  if (!"umap" %in% names(obj@reductions)) {
    stop("Missing UMAP reduction in annotated object: ", sample)
  }

  ## Compute all score matrices
  obj <- add_score_matrix(obj, build_initial_score_matrix(obj, initial_top50),
                          prefix = "initial_score_")
  obj <- add_score_matrix(obj, build_refined_score_matrix(obj, expression_filter_marker_tbl),
                          prefix = "refined_score_")
  obj <- add_score_matrix(obj, build_canonical_score_matrix(obj),
                          prefix = "canonical_score_")
  obj <- add_score_matrix(obj, build_filter_score_matrix(obj, expression_filter_marker_tbl),
                          prefix = "filter_score_")
  obj <- add_filter_status(obj)

  ## Page 1: Initial 3CA Top-50 Scores (Z-score, diverging blue-white-red)
  print(build_zscore_page(
    obj = obj, score_prefix = "initial_score_",
    page_title = pretty_sample_title(sample, ncol(obj), "Initial 3CA Top-50 Scores"),
    label_column = "celltype_initial", label_title = "Celltype Initial"
  ))

  ## Page 2: Refined Weighted Scores (Z-score, diverging blue-white-red)
  print(build_zscore_page(
    obj = obj, score_prefix = "refined_score_",
    page_title = pretty_sample_title(sample, ncol(obj), "Refined Weighted Scores"),
    label_column = "celltype_update", label_title = "Celltype Update"
  ))

  ## Page 3: Canonical Marker Expression (mean log-expression, grey-to-red)
  print(build_expression_page(
    obj = obj, score_prefix = "canonical_score_",
    page_title = pretty_sample_title(sample, ncol(obj), "Canonical Marker Expression"),
    label_column = "celltype_update", label_title = "Celltype Update"
  ))

  ## Page 4: Expression-Filter Panel Scores (top-4 mean expression, grey-to-red)
  ## Last panel shows filter status instead of celltype
  print(build_expression_page(
    obj = obj, score_prefix = "filter_score_",
    page_title = pretty_sample_title(sample, ncol(obj), "Expression-Filter Panel Scores"),
    label_column = "filter_status", label_title = "Filter Status",
    label_palette = filter_status_colors
  ))

  if (idx < length(selected_samples)) grid.newpage()

  plot_status[[sample]] <- data.frame(
    sample = sample, status = "written",
    n_cells = ncol(obj), stringsAsFactors = FALSE
  )
  rm(obj); gc(verbose = FALSE)
}

dev.off()
write.csv(bind_rows(plot_status), status_path, row.names = FALSE)
