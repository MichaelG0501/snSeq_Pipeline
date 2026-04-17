suppressPackageStartupMessages({
  library("Seurat")
  library("dplyr")
  library("Matrix")
  library("ComplexHeatmap")
  library("circlize")
  library("gridExtra")
  library("grid")
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/analysis/plotting/qc_heatmap.R")
}

analysis_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(analysis_dir, "..", ".."))
out_dir <- file.path(project_dir, "sn_outs")
setwd(out_dir)

required_files <- c(
  "sample_manifest.csv",
  "filtered_sample_summary.csv",
  "expression_filter_markers.csv"
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required QC heatmap inputs: ", paste(missing_files, collapse = ", "))
}

max_mt <- 15
min_ngenes <- 200
min_hk_expr <- 0.5

ct_reorder <- c(
  "b.cell", "dendritic", "endothelial", "epithelial",
  "fibroblast", "macrophage", "mast", "nk.cell", "nk.cell|t.cell",
  "plasma", "t.cell", "t.cell|nk.cell", "unresolved",
  "unresolved_inconsistent"
)

housekeeping <- c(
  "ACTB", "GAPDH", "RPS11", "RPS13", "RPS14", "RPS15", "RPS16", "RPS18",
  "RPS19", "RPS20", "RPL10", "RPL13", "RPL15", "RPL18"
)

load_marker_panel <- function(marker_path) {
  marker_tbl <- read.csv(marker_path, stringsAsFactors = FALSE)
  required_cols <- c("celltype", "gene", "rank")
  missing_cols <- setdiff(required_cols, colnames(marker_tbl))
  if (length(missing_cols) > 0) {
    stop("expression_filter_markers.csv is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  marker_tbl <- marker_tbl %>%
    mutate(
      celltype = as.character(celltype),
      gene = as.character(gene)
    ) %>%
    arrange(celltype, rank, gene)

  ordered_celltypes <- intersect(ct_reorder, unique(marker_tbl$celltype))
  marker_tbl <- marker_tbl %>%
    filter(celltype %in% ordered_celltypes)

  marker_list <- lapply(ordered_celltypes, function(ct) {
    unique(marker_tbl$gene[marker_tbl$celltype == ct])
  })
  names(marker_list) <- ordered_celltypes
  marker_list$housekeeping <- housekeeping

  marker_vec <- unique(c(unlist(marker_list, use.names = FALSE), housekeeping))

  list(marker_tbl = marker_tbl, markers_list = marker_list, markers = marker_vec)
}

marker_panel <- load_marker_panel("expression_filter_markers.csv")
markers_list <- marker_panel$markers_list
markers <- marker_panel$markers

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
    layer_data <- suppressWarnings(GetAssayData(obj, slot = layer))
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

pad_matrix <- function(mat, all_genes) {
  missing_genes <- setdiff(all_genes, rownames(mat))
  if (length(missing_genes) > 0) {
    zero_mat <- Matrix::Matrix(
      0,
      nrow = length(missing_genes),
      ncol = ncol(mat),
      dimnames = list(missing_genes, colnames(mat)),
      sparse = TRUE
    )
    mat <- rbind(mat, zero_mat)
  }
  mat[all_genes, , drop = FALSE]
}

normalize_label <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "unresolved"
  x
}

extract_stage_meta <- function(obj, stage_name, identity_col, sample_name = NULL) {
  meta <- obj@meta.data
  if (!identity_col %in% colnames(meta)) {
    stop("Missing identity column ", identity_col, " in stage object for ", if (!is.null(sample_name)) sample_name else stage_name)
  }
  meta$cell_id <- rownames(meta)
  meta$sample <- if (!is.null(sample_name)) sample_name else meta$orig.ident
  meta$plot_cell_id <- if (!is.null(sample_name)) {
    paste(sample_name, meta$cell_id, sep = "___")
  } else {
    meta$cell_id
  }
  meta$stage_name <- stage_name
  meta$celltype_plot <- normalize_label(meta[[identity_col]])
  meta[, c(
    "plot_cell_id", "cell_id", "sample", "orig.ident", "nFeature_RNA",
    "percent.mt", "celltype_plot", "stage_name"
  ), drop = FALSE]
}

combine_stage_expression <- function(stage_meta, sample_loader) {
  expr_parts <- list()
  meta_parts <- list()

  split_meta <- split(stage_meta, stage_meta$sample)
  for (sample_name in names(split_meta)) {
    obj <- sample_loader(sample_name)
    meta_sub <- split_meta[[sample_name]]
    expr <- get_layer(obj, "data")
    expr <- expr[intersect(markers, rownames(expr)), meta_sub$cell_id, drop = FALSE]
    expr <- pad_matrix(expr, unique(markers))
    colnames(expr) <- meta_sub$plot_cell_id
    rownames(meta_sub) <- meta_sub$plot_cell_id
    expr_parts[[sample_name]] <- expr
    meta_parts[[sample_name]] <- meta_sub
    rm(obj, expr, meta_sub)
    gc()
  }

  expr_mat <- do.call(cbind, expr_parts)
  meta_df <- bind_rows(meta_parts)
  rownames(meta_df) <- meta_df$plot_cell_id
  list(expr = expr_mat[, rownames(meta_df), drop = FALSE], meta = meta_df)
}

build_stage_from_anno <- function(stage_name, identity_col, cell_filter = NULL) {
  manifest <- read.csv("sample_manifest.csv", stringsAsFactors = FALSE)
  sample_ids <- unique(manifest$sample)

  meta_parts <- lapply(sample_ids, function(sample_name) {
    anno_path <- file.path("by_samples", sample_name, paste0(sample_name, "_anno.rds"))
    if (!file.exists(anno_path)) {
      return(NULL)
    }
    obj <- readRDS(anno_path)
    if (!is.null(cell_filter)) {
      keep_idx <- cell_filter(obj@meta.data)
      keep_idx[is.na(keep_idx)] <- FALSE
      if (!any(keep_idx)) {
        rm(obj)
        gc()
        return(NULL)
      }
      obj <- subset(obj, cells = rownames(obj@meta.data)[keep_idx])
    }
    meta <- extract_stage_meta(
      obj = obj,
      stage_name = stage_name,
      identity_col = identity_col,
      sample_name = sample_name
    )
    rm(obj)
    gc()
    meta
  })

  meta_all <- bind_rows(meta_parts)
  if (nrow(meta_all) == 0) {
    stop("No cells available for stage ", stage_name, ".")
  }

  combine_stage_expression(
    stage_meta = meta_all,
    sample_loader = function(sample_name) {
      readRDS(file.path("by_samples", sample_name, paste0(sample_name, "_anno.rds")))
    }
  )
}

build_prefilter_stage <- function() {
  build_stage_from_anno(
    stage_name = "post_qc_pre_expression_filter",
    identity_col = "celltype_update_step3"
  )
}

build_final_stage <- function() {
  sample_summary <- read.csv("filtered_sample_summary.csv", stringsAsFactors = FALSE)
  sample_ids <- unique(sample_summary$sample[sample_summary$status == "ok"])

  meta_parts <- lapply(sample_ids, function(sample_name) {
    filtered_path <- file.path("by_samples", sample_name, paste0(sample_name, "_filtered.rds"))
    if (!file.exists(filtered_path)) {
      return(NULL)
    }
    obj <- readRDS(filtered_path)
    meta <- extract_stage_meta(
      obj = obj,
      stage_name = "final_filtered_singlets",
      identity_col = "celltype_update",
      sample_name = sample_name
    )
    rm(obj)
    gc()
    meta
  })

  meta_all <- bind_rows(meta_parts)
  if (nrow(meta_all) == 0) {
    stop("No cells available for stage final_filtered_singlets.")
  }

  combine_stage_expression(
    stage_meta = meta_all,
    sample_loader = function(sample_name) {
      readRDS(file.path("by_samples", sample_name, paste0(sample_name, "_filtered.rds")))
    }
  )
}

plot_heatmap <- function(stage_obj, sampleid, reorder = TRUE, reorder_levels = ct_reorder) {
  expr_data <- as.data.frame(t(as.matrix(stage_obj$expr)))
  missing_markers <- setdiff(markers, colnames(expr_data))
  for (gene in missing_markers) {
    expr_data[[gene]] <- 0
  }
  expr_data <- expr_data[, markers, drop = FALSE]

  local_markers_list <- markers_list
  for (name in names(local_markers_list)) {
    local_markers_list[[name]] <- local_markers_list[[name]][local_markers_list[[name]] %in% colnames(expr_data)]
  }

  expr_data$celltype <- stage_obj$meta[rownames(expr_data), "celltype_plot", drop = TRUE]
  expr_data$celltype <- sapply(expr_data$celltype, function(x) gsub(" ", "\n", x))
  desired_order <- if (isTRUE(reorder)) reorder_levels else unique(expr_data$celltype)
  expr_data$celltype <- factor(expr_data$celltype, levels = c(desired_order, setdiff(unique(expr_data$celltype), desired_order)))
  expr_data <- expr_data[order(expr_data$celltype), , drop = FALSE]

  hk_avg <- rowMeans(expr_data[, colnames(expr_data) %in% local_markers_list$housekeeping, drop = FALSE])
  hk_avg <- matrix(hk_avg, nrow = 1, dimnames = list("avg_hk", names(hk_avg)))
  nCounts <- stage_obj$meta[rownames(expr_data), "nFeature_RNA", drop = TRUE]
  names(nCounts) <- rownames(expr_data)
  nCounts <- matrix(nCounts, nrow = 1, dimnames = list("nGenes", names(nCounts)))

  heatplot <- t(as.matrix(expr_data[, 1:(ncol(expr_data) - 1), drop = FALSE]))
  expr_max <- max(heatplot, na.rm = TRUE)
  if (!is.finite(expr_max) || expr_max <= 0) {
    expr_max <- 1
  }

  heatmap_grobs <- list()
  len <- length(unique(as.character(expr_data$celltype)))
  present_celltypes <- intersect(desired_order, unique(as.character(expr_data$celltype)))
  extra_celltypes <- setdiff(unique(as.character(expr_data$celltype)), present_celltypes)
  present_celltypes <- c(present_celltypes, sort(extra_celltypes))
  len <- length(present_celltypes)

  for (i in seq_along(local_markers_list)) {
    marker <- local_markers_list[[i]]
    marker_names <- rownames(heatplot[marker, , drop = FALSE])
    temp <- list()
    for (j in seq_along(present_celltypes)) {
      celltype <- present_celltypes[j]
      cells <- intersect(rownames(expr_data)[expr_data$celltype == celltype], colnames(heatplot))
      if (length(cells) == 0) {
        stop("No plotting columns found for cell type ", celltype, " in ", sampleid)
      }

      ht <- Heatmap(
        heatplot[marker, cells, drop = FALSE],
        col = colorRamp2(c(0, round(0.6 * expr_max, 1), ceiling(expr_max)),
                         c("#D0D0D0", "red4", "red4")),
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 40),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_heatmap_legend = FALSE,
        use_raster = FALSE
      )

      ht_grob <- grid.grabExpr(draw(ht, newpage = FALSE, padding = unit(c(2, 1, 2, 1), "mm")))
      temp[[j]] <- ht_grob
    }

    gene_label_grobs <- lapply(marker_names, function(name) {
      textGrob(label = name, x = unit(1, "npc"), just = "right", gp = gpar(fontsize = 40))
    })
    gene_label_col <- arrangeGrob(grobs = gene_label_grobs, ncol = 1)
    temp <- c(list(gene_label_col), temp)

    text_grob <- textGrob(
      names(local_markers_list)[i],
      gp = gpar(fontsize = 50, fontface = "bold")
    )
    rect_grob <- rectGrob(gp = gpar(fill = "grey", col = NA))
    merged_grob <- gTree(children = gList(rect_grob, text_grob))
    temp[[length(present_celltypes) + 2]] <- merged_grob

    combined_grob <- do.call(
      arrangeGrob,
      c(temp, list(ncol = len + 2, widths = c((len + 1) / 23, rep(1, len), (len + 1) / 12)))
    )
    heatmap_grobs[[length(heatmap_grobs) + 1]] <- combined_grob
  }

  stats_grobs <- list()
  for (i in seq_along(list(nCounts, hk_avg))) {
    data <- list(nCounts, hk_avg)[[i]]
    data_names <- rownames(data)
    temp <- list()

    for (j in seq_along(present_celltypes)) {
      celltype <- present_celltypes[j]
      cells <- intersect(rownames(expr_data)[expr_data$celltype == celltype], colnames(data))
      if (length(cells) == 0) {
        stop("No stats columns found for cell type ", celltype, " in ", sampleid)
      }

      if (i == 1) {
        gene_max <- max(nCounts, na.rm = TRUE)
        gene_high <- ceiling(gene_max / 1000) * 1000
        if (!is.finite(gene_high) || gene_high <= min_ngenes) {
          gene_high <- min_ngenes + 1000
        }
        ht_add <- Heatmap(
          data[, cells, drop = FALSE],
          name = "Number of\ngenes detected",
          col = colorRamp2(c(min_ngenes, round(0.7 * gene_max, -3), gene_high),
                           c("#D0D0D0", "blue3", "blue3")),
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = FALSE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 40),
          show_column_names = FALSE,
          show_heatmap_legend = FALSE,
          use_raster = FALSE
        )
      } else {
        hk_high <- ceiling(max(hk_avg, na.rm = TRUE) * 10) / 10
        hk_high <- max(hk_high, min_hk_expr)
        ht_add <- Heatmap(
          data[, cells, drop = FALSE],
          name = "Average\nhousekeeping\nexpression",
          col = colorRamp2(c(0, min_hk_expr, hk_high), c("#D0D0D0", "#A0B8E6", "blue3")),
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = FALSE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 40),
          show_column_names = FALSE,
          show_heatmap_legend = FALSE,
          use_raster = FALSE
        )
      }

      ht_grob <- grid.grabExpr(draw(ht_add, newpage = FALSE, padding = unit(c(6, 1.5, 6, 1.5), "mm")))
      temp[[j]] <- ht_grob
    }

    gene_label_grobs <- lapply(data_names, function(name) {
      textGrob(label = name, just = "center", gp = gpar(fontsize = 40))
    })
    gene_label_col <- arrangeGrob(grobs = gene_label_grobs, ncol = 1)
    temp <- c(list(gene_label_col), temp)

    text_grob <- textGrob(
      c("Number\nof genes", "Average\nhousekeeping\nexpression")[i],
      gp = gpar(fontsize = 40, fontface = "bold")
    )
    rect_grob <- rectGrob(gp = gpar(fill = "white", col = NA))
    merged_grob <- gTree(children = gList(rect_grob, text_grob))
    temp[[length(present_celltypes) + 2]] <- merged_grob

    combined_grob <- do.call(
      arrangeGrob,
      c(temp, list(ncol = len + 2, widths = c((len + 1) / 23, rep(1, len), (len + 1) / 12)))
    )
    stats_grobs[[length(stats_grobs) + 1]] <- combined_grob
  }

  celltype_labels <- c("", as.character(present_celltypes), "", "")
  text_grobs <- lapply(celltype_labels, function(celltype) {
    textGrob(celltype, gp = gpar(fontsize = 50, fontface = "bold"), just = "center", rot = 0)
  })
  text_grob_row <- do.call(
    arrangeGrob,
    c(text_grobs, list(ncol = len + 2, widths = c((len + 1) / 23, rep(1, len), (len + 1) / 12)))
  )

  celltype_num <- c("", as.numeric(table(expr_data$celltype)[present_celltypes]), "Markers", "")
  num_grobs <- lapply(celltype_num, function(celltype_num) {
    textGrob(celltype_num, gp = gpar(fontsize = 40, fontface = "bold"), just = "center")
  })
  num_grob_row <- do.call(
    arrangeGrob,
    c(num_grobs, list(ncol = len + 2, widths = c((len + 1) / 23, rep(1, len), (len + 1) / 12)))
  )

  title_grob <- textGrob(
    paste0(
      "Mitochondria DNA < ", max_mt,
      "     &&     Minimum number of genes > ", min_ngenes,
      "     &&     Minimum HK expression > ", min_hk_expr
    ),
    gp = gpar(fontsize = 40),
    just = "center"
  )

  gene_max <- max(nCounts, na.rm = TRUE)
  gene_high <- ceiling(gene_max / 1000) * 1000
  gene_high <- max(gene_high, min_ngenes + 1000)
  hk_high <- ceiling(max(hk_avg, na.rm = TRUE) * 10) / 10
  hk_high <- max(hk_high, min_hk_expr)

  legend_obj1 <- Legend(
    title = "\nExpression\nlevels (E)\n",
    at = pretty(c(0, ceiling(expr_max)), n = 5),
    labels_gp = gpar(fontsize = 50),
    title_gp = gpar(fontsize = 50, fontface = "bold"),
    grid_height = unit(220, "mm"),
    legend_width = unit(40, "mm"),
    legend_height = unit(220, "mm"),
    grid_width = unit(40, "mm"),
    title_position = "topcenter",
    col_fun = colorRamp2(c(0, round(0.6 * expr_max, 1), ceiling(expr_max)),
                         c("#D0D0D0", "red4", "red4"))
  )

  legend_obj2 <- Legend(
    title = "\nNumber\nof genes\n",
    at = round(seq(min_ngenes, gene_high, length.out = 4), -3),
    labels_gp = gpar(fontsize = 50),
    title_gp = gpar(fontsize = 50, fontface = "bold"),
    grid_height = unit(220, "mm"),
    legend_width = unit(40, "mm"),
    legend_height = unit(220, "mm"),
    grid_width = unit(40, "mm"),
    title_position = "topcenter",
    col_fun = colorRamp2(c(min_ngenes, round(0.7 * gene_max, -3), gene_high),
                         c("#D0D0D0", "blue3", "blue3"))
  )

  legend_obj3 <- Legend(
    title = "\nAverage\nhousekeeping\nexpression\n",
    at = seq(floor(min_hk_expr), hk_high, length.out = 4) %>% round(1),
    labels_gp = gpar(fontsize = 50),
    title_gp = gpar(fontsize = 50, fontface = "bold"),
    grid_height = unit(220, "mm"),
    legend_width = unit(40, "mm"),
    legend_height = unit(220, "mm"),
    grid_width = unit(40, "mm"),
    title_position = "topcenter",
    col_fun = colorRamp2(c(0, min_hk_expr, hk_high), c("#D0D0D0", "#A0B8E6", "blue3"))
  )

  main_content <- arrangeGrob(
    grobs = c(list(title_grob), list(text_grob_row), list(num_grob_row), heatmap_grobs, stats_grobs),
    ncol = 1,
    heights = c(2.5, 2, 1.5, rep(4, length(heatmap_grobs)), 4, 4)
  )

  legend_grob <- textGrob(label = sampleid, gp = gpar(fontsize = 35, fontface = "bold"))
  legend_grob0 <- textGrob(
    label = paste0("Total cell count: ", sum(expr_data$celltype %in% present_celltypes)),
    gp = gpar(fontsize = 35)
  )
  legend_grob1 <- grid.grabExpr(draw(legend_obj1))
  legend_grob2 <- grid.grabExpr(draw(legend_obj2))
  legend_grob3 <- grid.grabExpr(draw(legend_obj3))
  legend_column <- arrangeGrob(
    legend_grob, legend_grob0, legend_grob1, legend_grob2, legend_grob3,
    ncol = 1, heights = c(0.3, 0.2, 2, 2, 2)
  )

  grid.arrange(main_content, legend_column, ncol = 2, widths = c(9.3, 1))
}

prefilter_stage <- build_prefilter_stage()
final_stage <- build_final_stage()

summary_tbl <- bind_rows(
  prefilter_stage$meta %>%
    count(stage_name, celltype_plot, name = "n_plotted") %>%
    mutate(total_cells_stage = nrow(prefilter_stage$meta)),
  final_stage$meta %>%
    count(stage_name, celltype_plot, name = "n_plotted") %>%
    mutate(total_cells_stage = nrow(final_stage$meta))
)

write.csv(summary_tbl, "Auto_qc_heatmap_stage_summary.csv", row.names = FALSE)

png("Auto_QC_snSeq_prefilter.png", width = 80, height = 50, units = "in", res = 400)
plot_heatmap(
  stage_obj = prefilter_stage,
  sampleid = "snSeq cohort after QC, before expression filtering",
  reorder = TRUE
)
dev.off()

png("Auto_QC_snSeq_final.png", width = 80, height = 50, units = "in", res = 400)
plot_heatmap(
  stage_obj = final_stage,
  sampleid = "snSeq cohort final singlets with good expression",
  reorder = TRUE
)
dev.off()

pdf("Auto_QC_snSeq_heatmaps.pdf", width = 80, height = 50, onefile = TRUE, useDingbats = FALSE)
plot_heatmap(
  stage_obj = prefilter_stage,
  sampleid = "snSeq cohort after QC, before expression filtering",
  reorder = TRUE
)
grid.newpage()
plot_heatmap(
  stage_obj = final_stage,
  sampleid = "snSeq cohort final singlets with good expression",
  reorder = TRUE
)
dev.off()
