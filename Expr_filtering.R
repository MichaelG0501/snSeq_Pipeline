suppressPackageStartupMessages({
  library("Seurat")
  library("dplyr")
  library("tidyr")
  library("parallel")
  library("Matrix")
  library("ggplot2")
  library("patchwork")
})

reticulate::use_condaenv(
  "dmtcp",
  conda = "/rds/general/user/sg3723/home/anaconda3/bin/conda",
  required = TRUE
)

set.seed(1)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs")

target_celltypes <- c(
  "epithelial", "endothelial", "fibroblast", "macrophage", "b.cell",
  "plasma", "t.cell", "nk.cell", "mast", "dendritic"
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

canonical_backfill_priority <- list(
  epithelial = c("EPCAM", "CLDN7", "KRT19", "KRT7", "PIGR", "PLPP2", "LCN2"),
  endothelial = c("CDH5", "VWF", "CLDN5", "ENG", "CLEC14A"),
  fibroblast = c("DCN", "LUM", "COL1A1", "COL1A2", "COL3A1"),
  macrophage = c("CSF1R", "CD14", "AIF1", "TYROBP", "CD68"),
  b.cell = c("CD19", "CD22", "CD79A", "CD79B", "MS4A1"),
  plasma = c("DERL3", "MZB1", "JCHAIN", "FCRL5"),
  t.cell = c("CD2", "CD3E", "CD3G", "CD3D"),
  nk.cell = c("GNLY", "KLRB1", "NKG7", "PRF1", "GZMB"),
  mast = c("CPA3", "MS4A2", "HPGDS", "TPSB2", "TPSAB1"),
  dendritic = c("CLEC10A", "CCR7", "CD86")
)

target_panel_sizes <- c(
  epithelial = 7,
  endothelial = 8,
  fibroblast = 8,
  macrophage = 8,
  b.cell = 6,
  plasma = 6,
  t.cell = 6,
  nk.cell = 6,
  mast = 4,
  dendritic = 3
)

allowed_pairs <- list(
  c("t.cell", "nk.cell"),
  c("nk.cell", "t.cell"),
  c("b.cell", "plasma"),
  c("plasma", "b.cell")
)

pair_labels <- vapply(allowed_pairs, paste, collapse = "|", FUN.VALUE = character(1))

non_exclusive_groups <- list(
  c("macrophage", "mast", "dendritic"),
  c("t.cell", "nk.cell", "dendritic"),
  c("b.cell", "plasma")
)

gap_cut <- 0.8
n_workers <- 8
marker_top_n <- 4
expr_detect_cut <- 0.5
doublet_score_cut <- 1.0
doublet_detect_cut <- 2L
good_score_cut <- 1.25
good_detect_cut <- 2L
single_marker_score_cut <- 5
single_marker_detect_cut <- 1L
cluster_inconsistent_frac_cut <- 0.1
cluster_unresolved_count_cut <- 50L
cluster_min_cells_cut <- 100L
strong_weight_cut <- 5

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

normalize_pair_label <- function(ct) {
  ifelse(
    grepl("plasma\\|", ct, ignore.case = TRUE), "plasma",
    ifelse(
      grepl("b\\.cell\\|", ct, ignore.case = TRUE), "b.cell",
      ifelse(
        grepl("nk\\.cell\\|", ct, ignore.case = TRUE), "nk.cell",
        ifelse(grepl("t\\.cell\\|", ct, ignore.case = TRUE), "t.cell", ct)
      )
    )
  )
}

resolve_marker_col <- function(celltype_label, score_row) {
  if (is.na(celltype_label) || !nzchar(celltype_label)) {
    return(NA_character_)
  }

  pair <- strsplit(celltype_label, "\\|")[[1]]
  pair <- pair[nzchar(pair)]
  pair_key <- paste(pair, collapse = "|")

  if (length(pair) == 2 && pair_key %in% pair_labels) {
    pair <- intersect(pair, names(score_row))
    if (length(pair) == 0) {
      return(NA_character_)
    }
    return(names(which.max(score_row[pair])))
  }

  celltype_label
}

safe_scale <- function(x) {
  if (length(x) <= 1 || isTRUE(all.equal(stats::sd(x, na.rm = TRUE), 0))) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

build_expression_filter_markers <- function(marker_path) {
  if (!file.exists(marker_path)) {
    stop("Missing annotation_score_markers.csv. Run step 3 first.")
  }

  marker_tbl <- read.csv(marker_path, stringsAsFactors = FALSE)
  required_cols <- c("gene", "home_celltype", "home_pct", "off_pct_max", "score_weight")
  missing_cols <- setdiff(required_cols, colnames(marker_tbl))
  if (length(missing_cols) > 0) {
    stop(
      "annotation_score_markers.csv is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  marker_tbl <- marker_tbl %>%
    mutate(
      gene = as.character(gene),
      home_celltype = as.character(home_celltype)
    )

  panels <- vector("list", length(target_celltypes))
  names(panels) <- target_celltypes

  for (ct in target_celltypes) {
    derived <- marker_tbl %>%
      filter(home_celltype == ct) %>%
      arrange(desc(score_weight), desc(home_pct), off_pct_max, gene)

    strong_genes <- derived$gene[derived$score_weight >= strong_weight_cut]
    if (length(strong_genes) == 0 && nrow(derived) > 0) {
      strong_genes <- head(derived$gene, min(3, nrow(derived)))
    }

    backfill <- canonical_backfill_priority[[ct]]
    if (is.null(backfill)) {
      backfill <- canonical_markers[[ct]]
    }

    target_size <- unname(target_panel_sizes[ct])
    selected <- unique(c(strong_genes, backfill, derived$gene))
    selected <- head(selected, target_size)

    if (length(selected) < 3) {
      stop("Unable to build a >=3 marker expression-filter panel for cell type: ", ct)
    }

    panel_tbl <- lapply(seq_along(selected), function(rank_id) {
      gene <- selected[rank_id]
      derived_row <- derived[derived$gene == gene, , drop = FALSE]

      source <- if (nrow(derived_row) > 0 && gene %in% strong_genes) {
        "derived"
      } else if (gene %in% backfill && nrow(derived_row) > 0) {
        "backfill_promoted_derived"
      } else if (gene %in% backfill) {
        "canonical_backfill"
      } else {
        "derived_additional"
      }

      data.frame(
        celltype = ct,
        gene = gene,
        rank = rank_id,
        source = source,
        target_panel_size = target_size,
        home_pct = if (nrow(derived_row) > 0) derived_row$home_pct[1] else NA_real_,
        off_pct_max = if (nrow(derived_row) > 0) derived_row$off_pct_max[1] else NA_real_,
        score_weight = if (nrow(derived_row) > 0) derived_row$score_weight[1] else NA_real_,
        stringsAsFactors = FALSE
      )
    })

    panels[[ct]] <- bind_rows(panel_tbl)
  }

  bind_rows(panels)
}

score_markers <- function(tmdata, markers_list, top_n = marker_top_n, expr_cut = expr_detect_cut) {
  expr <- get_layer(tmdata, "data")
  all_markers <- unique(unlist(markers_list, use.names = FALSE))
  all_markers <- all_markers[all_markers %in% rownames(expr)]

  if (length(all_markers) == 0) {
    stop("None of the expression-filter markers are present in the assay.")
  }

  expr_sub <- expr[all_markers, , drop = FALSE]

  score_mtx <- matrix(
    0,
    nrow = ncol(expr_sub),
    ncol = length(markers_list),
    dimnames = list(colnames(expr_sub), names(markers_list))
  )

  detect_mtx <- matrix(
    0L,
    nrow = ncol(expr_sub),
    ncol = length(markers_list),
    dimnames = list(colnames(expr_sub), names(markers_list))
  )

  for (celltype in names(markers_list)) {
    marker_genes <- intersect(markers_list[[celltype]], rownames(expr_sub))
    if (length(marker_genes) == 0) {
      next
    }

    marker_expr <- as.matrix(expr_sub[marker_genes, , drop = FALSE])
    detect_mtx[, celltype] <- colSums(marker_expr > expr_cut, na.rm = TRUE)
    score_mtx[, celltype] <- apply(marker_expr, 2, function(v) {
      expressed <- sum(v > expr_cut, na.rm = TRUE)
      if (expressed == 0) {
        return(0)
      }
      k <- min(top_n, expressed, length(v))
      mean(sort(v, decreasing = TRUE)[seq_len(k)])
    })
  }

  list(score = score_mtx, detected = detect_mtx)
}

compute_top_two_diff <- function(score_row) {
  x <- as.numeric(score_row)
  names(x) <- names(score_row)

  if (!length(x) || all(x <= 0)) {
    return(Inf)
  }

  ord <- order(x, decreasing = TRUE)
  top_markers <- names(x)[ord]
  top1 <- top_markers[1]

  conflicting_groups <- non_exclusive_groups[vapply(
    non_exclusive_groups,
    function(g) top1 %in% g,
    logical(1)
  )]

  if (length(conflicting_groups) > 0) {
    conflicting_markers <- unique(unlist(conflicting_groups))
    other_idx <- which(!top_markers %in% conflicting_markers & x[top_markers] > 0)[1]
    if (is.na(other_idx)) {
      return(Inf)
    }
    return(x[top1] - x[top_markers[other_idx]])
  }

  positive_idx <- which(x[top_markers] > 0)
  if (length(positive_idx) < 2) {
    return(Inf)
  }

  x[top_markers[1]] - x[top_markers[2]]
}

passes_marker_gate <- function(score, detected) {
  if (is.na(score) || is.na(detected)) {
    return(FALSE)
  }

  (score >= good_score_cut && detected >= good_detect_cut) ||
    (score >= single_marker_score_cut && detected >= single_marker_detect_cut)
}

classify_filter_status <- function(celltype_labels, score_mtx, detect_mtx) {
  n_cells <- nrow(score_mtx)
  marker_expression <- rep(NA_character_, n_cells)
  assigned_score <- rep(NA_real_, n_cells)
  assigned_detected <- rep(NA_integer_, n_cells)
  max_filter_score <- rep(0, n_cells)
  max_filter_celltype <- rep(NA_character_, n_cells)
  max_filter_detected <- rep(0L, n_cells)

  for (i in seq_len(n_cells)) {
    row_vals <- score_mtx[i, , drop = TRUE]
    row_detect <- detect_mtx[i, , drop = TRUE]
    best_ct <- names(which.max(row_vals))

    max_filter_score[i] <- if (length(best_ct) > 0) row_vals[best_ct] else 0
    max_filter_celltype[i] <- if (length(best_ct) > 0) best_ct else NA_character_
    max_filter_detected[i] <- if (length(best_ct) > 0) row_detect[best_ct] else 0L

    celltype <- as.character(celltype_labels[i])
    is_resolved <- !is.na(celltype) && nzchar(celltype) && !(celltype %in% c("unresolved", "unknown"))

    if (is_resolved) {
      marker_col <- resolve_marker_col(celltype, row_vals)
      if (!is.na(marker_col) && marker_col %in% colnames(score_mtx)) {
        assigned_score[i] <- row_vals[marker_col]
        assigned_detected[i] <- row_detect[marker_col]
      }

      if (passes_marker_gate(assigned_score[i], assigned_detected[i])) {
        marker_expression[i] <- "good"
      } else if (passes_marker_gate(max_filter_score[i], max_filter_detected[i])) {
        marker_expression[i] <- "good_inconsistent"
      } else {
        marker_expression[i] <- "poor"
      }
    } else {
      marker_expression[i] <- if (passes_marker_gate(max_filter_score[i], max_filter_detected[i])) {
        "good_unresolved"
      } else {
        "poor_unresolved"
      }
    }
  }

  data.frame(
    marker_expression = marker_expression,
    marker_match_score = assigned_score,
    marker_match_detected = assigned_detected,
    max_filter_score = max_filter_score,
    max_filter_celltype = max_filter_celltype,
    max_filter_detected = max_filter_detected,
    stringsAsFactors = FALSE
  )
}

build_filter_reason <- function(coexpression_loose, marker_expression) {
  ifelse(
    coexpression_loose == "doublet", "doublet",
    ifelse(
      marker_expression == "good", "keep",
      ifelse(
        marker_expression == "good_inconsistent", "marker_inconsistent",
        ifelse(
          marker_expression == "good_unresolved", "marker_positive_unresolved",
          ifelse(
            marker_expression == "poor", "poor_assigned",
            ifelse(
              marker_expression == "poor_unresolved",
              "poor_unresolved",
              "other"
            )
          )
        )
      )
    )
  )
}

build_cluster_scores_long <- function(cluster_ids, score_mtx) {
  score_df <- as.data.frame(score_mtx, stringsAsFactors = FALSE)
  score_df$cluster <- as.character(cluster_ids)

  score_df %>%
    group_by(cluster) %>%
    summarize(across(all_of(colnames(score_mtx)), median, na.rm = TRUE), .groups = "drop") %>%
    pivot_longer(-cluster, names_to = "cell_type", values_to = "score") %>%
    group_by(cluster) %>%
    mutate(z = safe_scale(score)) %>%
    ungroup()
}

call_cluster_from_scores <- function(scores_long) {
  cluster_ids <- unique(scores_long$cluster)

  bind_rows(lapply(cluster_ids, function(cluster_id) {
    sub <- scores_long[scores_long$cluster == cluster_id, , drop = FALSE]
    az <- sub$z
    act <- sub$cell_type
    max_idx <- which.max(az)
    top_ct <- act[max_idx]
    top_z <- az[max_idx]
    other_idx <- setdiff(seq_along(az), max_idx)

    if (length(other_idx) == 0) {
      step2 <- top_ct
    } else {
      margins <- top_z - az[other_idx]
      close_idx <- other_idx[margins < gap_cut]

      if (length(close_idx) == 0) {
        step2 <- top_ct
      } else if (length(close_idx) == 1) {
        other_ct <- act[close_idx]
        pair_vec <- c(top_ct, other_ct)
        is_allowed <- any(vapply(allowed_pairs, function(p) identical(p, pair_vec), logical(1)))
        step2 <- if (is_allowed) paste0(top_ct, "|", other_ct) else "unresolved"
      } else {
        step2 <- "unresolved"
      }
    }

    data.frame(cluster = cluster_id, step2 = step2, stringsAsFactors = FALSE)
  }))
}

manifest_path <- "sample_manifest.csv"
if (!file.exists(manifest_path)) {
  stop("Missing sample_manifest.csv. Run step 1 first.")
}

filter_marker_tbl <- build_expression_filter_markers("annotation_score_markers.csv")
write.csv(filter_marker_tbl, "expression_filter_markers.csv", row.names = FALSE)
write.csv(
  filter_marker_tbl %>%
    count(celltype, name = "n_markers"),
  "expression_filter_marker_counts.csv",
  row.names = FALSE
)

markers_list <- split(filter_marker_tbl$gene, filter_marker_tbl$celltype)
marker_columns <- names(markers_list)

sample_manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
required_manifest_cols <- c("sample", "reference_batch", "technology", "input_source")
missing_manifest_cols <- setdiff(required_manifest_cols, colnames(sample_manifest))
if (length(missing_manifest_cols) > 0) {
  stop(
    "sample_manifest.csv is missing required columns: ",
    paste(missing_manifest_cols, collapse = ", "),
    ". Rerun step 1."
  )
}

sample_dirs <- unique(sample_manifest$sample)
tmdata_annotated <- list()

for (sample in sample_dirs) {
  anno_path <- file.path("by_samples", sample, paste0(sample, "_anno.rds"))
  status_path <- file.path("by_samples", sample, "expr_filter_status.txt")

  if (!file.exists(anno_path)) {
    writeLines("missing_anno", status_path)
    next
  }

  tmdata <- readRDS(anno_path)
  required_meta_cols <- c(
    "orig.ident", "celltype_update", "celltype_initial",
    "reference_batch", "technology", "input_source"
  )
  missing_meta_cols <- setdiff(required_meta_cols, colnames(tmdata@meta.data))
  if (length(missing_meta_cols) > 0) {
    stop(
      "Missing metadata columns in ", sample, ": ",
      paste(missing_meta_cols, collapse = ", "),
      ". Rerun earlier steps."
    )
  }

  drop_cols <- intersect(
    c(
      marker_columns,
      "coexpression", "coexpression_loose", "marker_expression",
      "marker_match_score", "marker_match_detected",
      "max_filter_score", "max_filter_celltype", "max_filter_detected",
      "filter_keep", "filter_reason", "expr_filter_reclustered",
      "celltype_update_step3"
    ),
    colnames(tmdata@meta.data)
  )

  if (length(drop_cols) > 0) {
    tmdata@meta.data[, drop_cols] <- NULL
  }

  tmdata$celltype_update <- normalize_pair_label(as.character(tmdata$celltype_update))
  tmdata$celltype_update_step3 <- tmdata$celltype_update
  tmdata$expr_filter_reclustered <- FALSE
  tmdata_annotated[[sample]] <- tmdata
}

for (name in names(tmdata_annotated)) {
  scored <- score_markers(tmdata_annotated[[name]], markers_list)
  score_all <- scored$score
  detected_all <- scored$detected

  high_markers <- score_all > doublet_score_cut & detected_all >= doublet_detect_cut
  high_marker_count <- rowSums(high_markers, na.rm = TRUE)

  special_case <- apply(high_markers, 1, function(cell_row) {
    high_ct <- names(which(cell_row))
    if (length(high_ct) == 0) {
      return(FALSE)
    }
    any(vapply(non_exclusive_groups, function(group) {
      all(high_ct %in% group) && length(high_ct) <= length(group)
    }, logical(1)))
  })

  tmdata_annotated[[name]]$coexpression <- ifelse(
    high_marker_count > 1 & !special_case,
    "doublet",
    "singlet"
  )

  gap_input <- score_all
  gap_input[detected_all < doublet_detect_cut] <- 0
  top_two_diff <- apply(gap_input, 1, compute_top_two_diff)

  tmdata_annotated[[name]]$coexpression_loose <- ifelse(
    (high_marker_count > 1) & (top_two_diff < 1),
    "doublet",
    "singlet"
  )

  filter_status <- classify_filter_status(
    tmdata_annotated[[name]]$celltype_update,
    score_all,
    detected_all
  )

  tmdata_annotated[[name]]$marker_expression <- filter_status$marker_expression
  tmdata_annotated[[name]]$marker_match_score <- filter_status$marker_match_score
  tmdata_annotated[[name]]$marker_match_detected <- filter_status$marker_match_detected
  tmdata_annotated[[name]]$max_filter_score <- filter_status$max_filter_score
  tmdata_annotated[[name]]$max_filter_celltype <- filter_status$max_filter_celltype
  tmdata_annotated[[name]]$max_filter_detected <- filter_status$max_filter_detected
  tmdata_annotated[[name]]$filter_reason <- build_filter_reason(
    tmdata_annotated[[name]]$coexpression_loose,
    tmdata_annotated[[name]]$marker_expression
  )
  tmdata_annotated[[name]]$filter_keep <- tmdata_annotated[[name]]$filter_reason == "keep"
}

for (name in names(tmdata_annotated)) {
  message("checking suspicious clusters in sample: ", name)

  mask <- tmdata_annotated[[name]]$coexpression == "singlet" &
    !(tmdata_annotated[[name]]$marker_expression %in% c("poor", "poor_unresolved"))

  marker_expression_table <- table(
    tmdata_annotated[[name]]$seurat_clusters[mask],
    tmdata_annotated[[name]]$marker_expression[mask]
  )

  cluster_tbc <- character(0)

  if ("good_inconsistent" %in% colnames(marker_expression_table)) {
    frac_inconsistent <- marker_expression_table[, "good_inconsistent"] / rowSums(marker_expression_table)
    cluster_tbc <- c(
      cluster_tbc,
      rownames(marker_expression_table)[
        frac_inconsistent > cluster_inconsistent_frac_cut &
          rowSums(marker_expression_table) >= cluster_min_cells_cut
      ]
    )
  }

  if ("good_unresolved" %in% colnames(marker_expression_table)) {
    cluster_tbc <- c(
      cluster_tbc,
      rownames(marker_expression_table)[
        marker_expression_table[, "good_unresolved"] > cluster_unresolved_count_cut
      ]
    )
  }

  cluster_tbc <- unique(cluster_tbc)

  if (length(cluster_tbc) == 0) {
    next
  }

  for (cluster in cluster_tbc) {
    tm_sub <- subset(
      tmdata_annotated[[name]],
      subset = seurat_clusters == cluster & coexpression == "singlet"
    )

    ncells <- ncol(tm_sub)
    if (ncells < 3) {
      next
    }

    tm_sub <- FindVariableFeatures(tm_sub, nfeatures = 3000, verbose = FALSE)
    tm_sub <- ScaleData(tm_sub, verbose = FALSE)
    max_pcs <- min(ncells, nrow(tm_sub), 50)
    tm_sub <- RunPCA(tm_sub, npcs = max_pcs, verbose = FALSE)
    pc_use <- 1:min(30, max_pcs)
    k_use <- min(30, max(10, round(sqrt(ncells) - 2)))
    tm_sub <- FindNeighbors(tm_sub, dims = pc_use, k.param = k_use, verbose = FALSE)
    tm_sub <- FindClusters(tm_sub, resolution = 0.8, algorithm = 4, verbose = FALSE)
    tm_sub <- RunUMAP(tm_sub, dims = pc_use, n.neighbors = 30, min.dist = 0.3, verbose = FALSE)

    rescored <- score_markers(tm_sub, markers_list)
    sub_scores_long <- build_cluster_scores_long(tm_sub$seurat_clusters, rescored$score)
    step2_calls <- call_cluster_from_scores(sub_scores_long)

    cl_map <- step2_calls$step2
    names(cl_map) <- step2_calls$cluster

    tm_sub$celltype_update <- as.character(cl_map[as.character(tm_sub$seurat_clusters)])
    tm_sub$celltype_update <- normalize_pair_label(tm_sub$celltype_update)

    tm_sub_status <- classify_filter_status(
      tm_sub$celltype_update,
      rescored$score,
      rescored$detected
    )

    tm_sub$marker_expression <- tm_sub_status$marker_expression
    tm_sub$marker_match_score <- tm_sub_status$marker_match_score
    tm_sub$marker_match_detected <- tm_sub_status$marker_match_detected
    tm_sub$max_filter_score <- tm_sub_status$max_filter_score
    tm_sub$max_filter_celltype <- tm_sub_status$max_filter_celltype
    tm_sub$max_filter_detected <- tm_sub_status$max_filter_detected
    tm_sub$filter_reason <- build_filter_reason(tm_sub$coexpression_loose, tm_sub$marker_expression)
    tm_sub$filter_keep <- tm_sub$filter_reason == "keep"
    tm_sub$expr_filter_reclustered <- TRUE

    update_cols <- c(
      "celltype_update", "marker_expression",
      "marker_match_score", "marker_match_detected",
      "max_filter_score", "max_filter_celltype", "max_filter_detected",
      "filter_reason", "filter_keep", "expr_filter_reclustered"
    )

    tmdata_annotated[[name]]@meta.data[colnames(tm_sub), update_cols] <-
      tm_sub@meta.data[, update_cols, drop = FALSE]
  }
}

mclapply(names(tmdata_annotated), function(nm) {
  sample_dir <- file.path("by_samples", nm)
  if (!dir.exists(sample_dir)) {
    dir.create(sample_dir, recursive = TRUE)
  }
  saveRDS(
    tmdata_annotated[[nm]],
    file.path(sample_dir, paste0(nm, "_anno.rds")),
    compress = FALSE
  )
  NULL
}, mc.cores = n_workers, mc.preschedule = FALSE)

filtered_objects <- list()
filtered_summary <- list()

for (name in names(tmdata_annotated)) {
  status_path <- file.path("by_samples", name, "expr_filter_status.txt")
  keep_cells <- colnames(tmdata_annotated[[name]])[tmdata_annotated[[name]]$filter_keep]

  filtered_obj <- if (length(keep_cells) > 1) {
    subset(tmdata_annotated[[name]], cells = keep_cells)
  } else {
    NULL
  }

  if (is.null(filtered_obj) || ncol(filtered_obj) <= 1) {
    if (length(keep_cells) > 0) {
      keep_idx <- colnames(tmdata_annotated[[name]]) %in% keep_cells
      tmdata_annotated[[name]]$filter_keep[keep_idx] <- FALSE
      tmdata_annotated[[name]]$filter_reason[keep_idx] <- "too_few_cells"
    }

    writeLines("no_cell", status_path)
    saveRDS("No cells", file.path("by_samples", name, "no_cell"))
    filtered_summary[[name]] <- data.frame(
      sample = name,
      reference_batch = unique(tmdata_annotated[[name]]$reference_batch)[1],
      technology = unique(tmdata_annotated[[name]]$technology)[1],
      input_source = unique(tmdata_annotated[[name]]$input_source)[1],
      cells_before = ncol(tmdata_annotated[[name]]),
      cells_after = 0,
      cells_removed = ncol(tmdata_annotated[[name]]),
      pct_retained = 0,
      epithelial_after = 0,
      status = "no_cell",
      stringsAsFactors = FALSE
    )
    next
  }

  writeLines("ok", status_path)
  saveRDS(
    filtered_obj,
    file.path("by_samples", name, paste0(name, "_filtered.rds")),
    compress = FALSE
  )

  filtered_objects[[name]] <- filtered_obj
  filtered_summary[[name]] <- data.frame(
    sample = name,
    reference_batch = unique(filtered_obj$reference_batch)[1],
    technology = unique(filtered_obj$technology)[1],
    input_source = unique(filtered_obj$input_source)[1],
    cells_before = ncol(tmdata_annotated[[name]]),
    cells_after = ncol(filtered_obj),
    cells_removed = ncol(tmdata_annotated[[name]]) - ncol(filtered_obj),
    pct_retained = (ncol(filtered_obj) / ncol(tmdata_annotated[[name]])) * 100,
    epithelial_after = sum(filtered_obj$celltype_update == "epithelial", na.rm = TRUE),
    status = "ok",
    stringsAsFactors = FALSE
  )
}

if (length(filtered_summary) == 0) {
  stop("No samples were processed in expression filtering.")
}

filtered_sample_summary <- bind_rows(filtered_summary)
write.csv(filtered_sample_summary, "filtered_sample_summary.csv", row.names = FALSE)

status_summary <- bind_rows(lapply(names(tmdata_annotated), function(name) {
  tmdata <- tmdata_annotated[[name]]
  sample_after <- filtered_sample_summary$cells_after[filtered_sample_summary$sample == name]
  if (length(sample_after) == 0) {
    sample_after <- 0L
  }

  data.frame(
    sample = name,
    reference_batch = unique(tmdata$reference_batch)[1],
    technology = unique(tmdata$technology)[1],
    input_source = unique(tmdata$input_source)[1],
    cells_total = ncol(tmdata),
    singlet_strict = sum(tmdata$coexpression == "singlet", na.rm = TRUE),
    doublet_strict = sum(tmdata$coexpression == "doublet", na.rm = TRUE),
    singlet_loose = sum(tmdata$coexpression_loose == "singlet", na.rm = TRUE),
    doublet_loose = sum(tmdata$coexpression_loose == "doublet", na.rm = TRUE),
    good = sum(tmdata$marker_expression == "good", na.rm = TRUE),
    good_inconsistent = sum(tmdata$marker_expression == "good_inconsistent", na.rm = TRUE),
    good_unresolved = sum(tmdata$marker_expression == "good_unresolved", na.rm = TRUE),
    poor = sum(tmdata$marker_expression == "poor", na.rm = TRUE),
    poor_unresolved = sum(tmdata$marker_expression == "poor_unresolved", na.rm = TRUE),
    kept_cells = as.integer(sample_after[1]),
    reclustered_cells = sum(tmdata$expr_filter_reclustered, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write.csv(status_summary, "expression_filter_status_by_sample.csv", row.names = FALSE)

reason_by_sample <- bind_rows(lapply(names(tmdata_annotated), function(name) {
  tmdata <- tmdata_annotated[[name]]
  counts <- as.data.frame(table(tmdata$filter_reason), stringsAsFactors = FALSE)
  colnames(counts) <- c("filter_reason", "n_cells")
  counts$sample <- name
  counts$reference_batch <- unique(tmdata$reference_batch)[1]
  counts$technology <- unique(tmdata$technology)[1]
  counts$input_source <- unique(tmdata$input_source)[1]
  counts
}))
write.csv(reason_by_sample, "expression_filter_reason_by_sample.csv", row.names = FALSE)
write.csv(
  reason_by_sample %>%
    group_by(filter_reason) %>%
    summarize(n_cells = sum(n_cells), .groups = "drop") %>%
    arrange(desc(n_cells)),
  "expression_filter_reason_overall.csv",
  row.names = FALSE
)

celltypes_ordered <- c("macrophage", "endothelial", "fibroblast", "t.cell")
ct_col <- "celltype_update"

reference_summary <- list()
ref_batches <- unique(filtered_sample_summary$reference_batch[filtered_sample_summary$status == "ok"])

for (sid in ref_batches) {
  sample_ids <- filtered_sample_summary$sample[
    filtered_sample_summary$status == "ok" & filtered_sample_summary$reference_batch == sid
  ]
  sample_ids <- sample_ids[sample_ids %in% names(filtered_objects)]

  if (length(sample_ids) == 0) {
    next
  }

  message("writing reference for ", sid)

  all_genes <- Reduce(union, lapply(sample_ids, function(id) {
    rownames(get_layer(filtered_objects[[id]], "counts"))
  }))

  cpm_list <- lapply(sample_ids, function(id) {
    mat <- get_layer(filtered_objects[[id]], "CPM")
    colnames(mat) <- paste(id, colnames(mat), sep = "_")
    pad_matrix(mat, all_genes)
  })
  names(cpm_list) <- sample_ids

  meta_list <- lapply(sample_ids, function(id) {
    md <- filtered_objects[[id]]@meta.data[
      ,
      c(
        "orig.ident", "celltype_update", "coexpression", "celltype_initial",
        "reference_batch", "technology", "input_source"
      ),
      drop = FALSE
    ]
    rownames(md) <- paste(id, rownames(md), sep = "_")
    md
  })
  names(meta_list) <- sample_ids

  mat <- do.call(cbind, cpm_list)
  expected_ref_cells <- unlist(lapply(cpm_list, colnames), use.names = FALSE)
  if (ncol(mat) != length(expected_ref_cells)) {
    stop("Reference CPM matrices are not aligned in column count for batch: ", sid)
  }
  colnames(mat) <- expected_ref_cells
  meta <- do.call(rbind, meta_list)
  expected_meta_cells <- unlist(lapply(meta_list, rownames), use.names = FALSE)
  if (nrow(meta) != length(expected_meta_cells)) {
    stop("Reference metadata frames are not aligned in row count for batch: ", sid)
  }
  rownames(meta) <- expected_meta_cells
  meta <- subset(meta, coexpression == "singlet" & celltype_initial == celltype_update)

  existing_ct <- intersect(celltypes_ordered, unique(meta[[ct_col]]))
  if (length(existing_ct) < 2) {
    for (sample in sample_ids) {
      saveRDS("Not enough reference cell types", file.path("by_samples", sample, "no_ref"))
    }

    reference_summary[[sid]] <- data.frame(
      reference_batch = sid,
      samples = length(sample_ids),
      ref_types = length(existing_ct),
      ref_ct = if (length(existing_ct) > 0) paste(existing_ct, collapse = "|") else "",
      ref_cells = 0,
      status = "skipped",
      stringsAsFactors = FALSE
    )
    next
  }

  ref_ct <- head(existing_ct, 2)
  cells <- intersect(
    rownames(meta)[meta[[ct_col]] %in% ref_ct],
    colnames(mat)
  )

  if (length(cells) == 0) {
    for (sample in sample_ids) {
      saveRDS("No shared reference cells after alignment", file.path("by_samples", sample, "no_ref"))
    }

    reference_summary[[sid]] <- data.frame(
      reference_batch = sid,
      samples = length(sample_ids),
      ref_types = length(ref_ct),
      ref_ct = paste(ref_ct, collapse = "|"),
      ref_cells = 0,
      status = "skipped_no_shared_cells",
      stringsAsFactors = FALSE
    )
    next
  }

  mat_ref <- mat[, cells, drop = FALSE]
  meta_ref <- meta[cells, c("orig.ident", "celltype_update"), drop = FALSE]
  ref <- setNames(lapply(ref_ct, function(ct) {
    rownames(meta_ref)[meta_ref[[ct_col]] == ct]
  }), ref_ct)

  saveRDS(
    list(
      matrix = mat_ref,
      ref_ct = ref_ct,
      ref = ref,
      meta = meta_ref,
      reference_batch = sid
    ),
    paste0(sid, "_reference.rds")
  )

  reference_summary[[sid]] <- data.frame(
    reference_batch = sid,
    samples = length(sample_ids),
    ref_types = length(ref_ct),
    ref_ct = paste(ref_ct, collapse = "|"),
    ref_cells = ncol(mat_ref),
    status = "written",
    stringsAsFactors = FALSE
  )
}

if (length(reference_summary) > 0) {
  write.csv(bind_rows(reference_summary), "reference_summary.csv", row.names = FALSE)
}

if (length(filtered_objects) == 0) {
  stop("No samples passed expression filtering.")
}

all_genes <- Reduce(union, lapply(filtered_objects, function(obj) {
  rownames(get_layer(obj, "counts"))
}))

counts_list <- lapply(names(filtered_objects), function(id) {
  mat <- get_layer(filtered_objects[[id]], "counts")
  colnames(mat) <- paste(id, colnames(mat), sep = "_")
  pad_matrix(mat, all_genes)
})
names(counts_list) <- names(filtered_objects)

cpm_list <- lapply(names(filtered_objects), function(id) {
  mat <- get_layer(filtered_objects[[id]], "CPM")
  colnames(mat) <- paste(id, colnames(mat), sep = "_")
  pad_matrix(mat, all_genes)
})
names(cpm_list) <- names(filtered_objects)

lognorm_list <- lapply(names(filtered_objects), function(id) {
  mat <- get_layer(filtered_objects[[id]], "data")
  colnames(mat) <- paste(id, colnames(mat), sep = "_")
  pad_matrix(mat, all_genes)
})
names(lognorm_list) <- names(filtered_objects)

meta_list <- lapply(names(filtered_objects), function(id) {
  meta <- filtered_objects[[id]]@meta.data
  rownames(meta) <- paste(id, rownames(meta), sep = "_")
  meta
})
names(meta_list) <- names(filtered_objects)

combined_counts <- do.call(cbind, counts_list)
combined_cpm <- do.call(cbind, cpm_list)
combined_lognorm <- do.call(cbind, lognorm_list)
combined_meta <- do.call(rbind, meta_list)
expected_meta_cells <- unlist(lapply(meta_list, rownames), use.names = FALSE)
if (nrow(combined_meta) != length(expected_meta_cells)) {
  stop("Filtered metadata frames are not aligned in row count.")
}
rownames(combined_meta) <- expected_meta_cells

expected_cells <- unlist(lapply(counts_list, colnames), use.names = FALSE)

if (ncol(combined_counts) != length(expected_cells) ||
    ncol(combined_cpm) != length(expected_cells) ||
    ncol(combined_lognorm) != length(expected_cells)) {
  stop("Filtered count/data/CPM matrices are not aligned in column count.")
}

colnames(combined_counts) <- expected_cells
colnames(combined_cpm) <- expected_cells
colnames(combined_lognorm) <- expected_cells

shared_cells <- intersect(expected_cells, rownames(combined_meta))

if (length(shared_cells) == 0) {
  stop("No shared cells remain after aligning filtered matrices and metadata.")
}

combined_counts <- combined_counts[, shared_cells, drop = FALSE]
combined_cpm <- combined_cpm[, shared_cells, drop = FALSE]
combined_lognorm <- combined_lognorm[, shared_cells, drop = FALSE]
combined_meta <- combined_meta[shared_cells, , drop = FALSE]

merged_obj <- CreateSeuratObject(counts = combined_counts, meta.data = combined_meta)
merged_obj@assays$RNA$CPM <- combined_cpm
merged_obj@assays$RNA$data <- combined_lognorm

if (!"technology" %in% colnames(merged_obj@meta.data)) {
  stop("Missing technology metadata in merged object. Rerun step 1.")
}
if (!"reference_batch" %in% colnames(merged_obj@meta.data)) {
  stop("Missing reference_batch metadata in merged object. Rerun step 1.")
}

merged_obj$study <- merged_obj$reference_batch
merged_obj$timepoint <- ifelse(
  grepl("pre", merged_obj$orig.ident, ignore.case = TRUE),
  "pre",
  ifelse(grepl("post", merged_obj$orig.ident, ignore.case = TRUE), "post", NA)
)
merged_obj$tissue_type <- ifelse(grepl("N1|N2", merged_obj$orig.ident), "LN", "Tumour")

write.csv(
  as.data.frame(table(merged_obj$celltype_update), stringsAsFactors = FALSE) %>%
    rename(celltype_update = Var1, n_cells = Freq) %>%
    arrange(desc(n_cells)),
  "filtered_celltype_counts.csv",
  row.names = FALSE
)

if (ncol(merged_obj) > 20) {
  merged_obj <- FindVariableFeatures(merged_obj, verbose = FALSE)
  merged_obj <- ScaleData(merged_obj, verbose = FALSE)
  merged_obj <- RunPCA(merged_obj, verbose = FALSE)
  merged_obj <- FindNeighbors(merged_obj, dims = 1:20, verbose = FALSE)
  merged_obj <- FindClusters(merged_obj, resolution = 0.8, algorithm = 1, verbose = FALSE)
  merged_obj$leiden_clusters <- Idents(merged_obj)
  merged_obj <- RunUMAP(merged_obj, dims = 1:20, verbose = FALSE)

  p1 <- DimPlot(merged_obj, group.by = "leiden_clusters", label = TRUE) + ggtitle("Louvain Clustering")
  p2 <- DimPlot(merged_obj, group.by = "study", label = FALSE) + ggtitle("Reference Batch")
  p3 <- DimPlot(merged_obj, group.by = "celltype_update", label = FALSE) + ggtitle("Celltype")
  p4 <- DimPlot(merged_obj, group.by = "technology", label = FALSE) + ggtitle("Technology")
  ggsave("snSeq_filtered.png", plot = (p1 | p2) / (p3 | p4), width = 12, height = 6, dpi = 300)
}

saveRDS(merged_obj, "snSeq_merged.rds")
