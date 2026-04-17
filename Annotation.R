suppressPackageStartupMessages({
  library("Seurat")
  library("dplyr")
  library("purrr")
  library("tidyr")
  library("ggplot2")
  library("readxl")
})

set.seed(1)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs")

manifest_path <- "sample_manifest.csv"
if (!file.exists(manifest_path)) {
  stop("Missing sample_manifest.csv. Run step 1 first.")
}

sample_manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
sample_dirs <- unique(sample_manifest$sample)
tmdata_list <- list()

for (sample in sample_dirs) {
  rds_path <- file.path("by_samples", sample, paste0(sample, ".rds"))
  if (!file.exists(rds_path)) {
    stop("Missing sample RDS: ", rds_path)
  }

  tmdata_list[[sample]] <- readRDS(rds_path)
  if (!"celltype_initial" %in% colnames(tmdata_list[[sample]]@meta.data)) {
    stop("Missing celltype_initial in sample: ", sample)
  }
  if (ncol(tmdata_list[[sample]]) < 50) {
    tmdata_list[[sample]]$celltype_update <- rep("unresolved", ncol(tmdata_list[[sample]]))
  }
}

markers <- read_excel("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Marker_Genes.xlsx", sheet = 1)
celltype_map <- c(
  "Fibroblast" = "fibroblast",
  "Macrophage" = "macrophage",
  "Mast" = "mast",
  "B cell" = "b.cell",
  "T cell" = "t.cell",
  "Dendritic" = "dendritic",
  "Endothelial" = "endothelial",
  "Epithelial" = "epithelial",
  "NK cell" = "nk.cell",
  "Plasma" = "plasma"
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

combine_marker_scores <- function(df, w_specificity = 0.2, w_sensitivity = 0.8) {
  pr <- function(x) {
    r <- rank(x, ties.method = "average", na.last = "keep")
    r / (sum(!is.na(x)) + 1)
  }

  combined <- (w_specificity * pr(df$specificity) + w_sensitivity * pr(df$sensitivity)) /
    (w_specificity + w_sensitivity)

  df %>% mutate(Combined = combined) %>% arrange(desc(Combined))
}

build_initial_seed <- function(obj, sample_id, valid_celltypes, expr_cut = 1, min_gap = 0.1) {
  initial <- as.character(obj$celltype_initial)
  seed <- rep("unresolved", length(initial))
  seed_top <- rep(NA_character_, length(initial))
  seed_top_prop <- rep(NA_real_, length(initial))
  seed_gap <- rep(NA_real_, length(initial))

  names(seed) <- colnames(obj)
  names(seed_top) <- colnames(obj)
  names(seed_top_prop) <- colnames(obj)
  names(seed_gap) <- colnames(obj)

  marker_sets <- canonical_markers[intersect(valid_celltypes, names(canonical_markers))]
  marker_sets <- marker_sets[vapply(marker_sets, function(gs) {
    length(intersect(gs, rownames(obj))) > 0
  }, logical(1))]

  if (length(marker_sets) > 0) {
    marker_sets <- lapply(marker_sets, intersect, rownames(obj))
    seed_genes <- unique(unlist(marker_sets))
    expr <- GetAssayData(obj, slot = "data")[seed_genes, , drop = FALSE]
    prop_mat <- do.call(cbind, lapply(marker_sets, function(gs) {
      Matrix::colMeans(expr[gs, , drop = FALSE] > expr_cut)
    }))
    if (is.null(dim(prop_mat))) {
      prop_mat <- matrix(prop_mat, ncol = 1)
    }
    rownames(prop_mat) <- colnames(obj)
    colnames(prop_mat) <- names(marker_sets)

    valid_cells <- which(initial %in% colnames(prop_mat))
    for (idx in valid_cells) {
      vals <- prop_mat[colnames(obj)[idx], , drop = TRUE]
      max_prop <- max(vals, na.rm = TRUE)
      if (!is.finite(max_prop) || max_prop <= 0) {
        next
      }

      top_ct <- names(vals)[vals == max_prop]
      sorted_vals <- sort(vals, decreasing = TRUE)
      gap <- if (length(sorted_vals) >= 2) sorted_vals[1] - sorted_vals[2] else sorted_vals[1]

      seed_top[idx] <- paste(top_ct, collapse = "|")
      seed_top_prop[idx] <- max_prop
      seed_gap[idx] <- gap

      if (length(top_ct) == 1 && identical(top_ct, initial[idx]) && is.finite(gap) && gap >= min_gap) {
        seed[idx] <- initial[idx]
      }
    }
  }

  obj$celltype_initial_seed <- seed
  obj$canonical_seed_top <- seed_top
  obj$canonical_seed_prop <- seed_top_prop
  obj$canonical_seed_gap <- seed_gap

  summary <- tibble(
    sample = sample_id,
    celltype_initial = initial,
    celltype_initial_seed = seed
  ) %>%
    filter(celltype_initial %in% valid_celltypes) %>%
    count(sample, celltype_initial, celltype_initial_seed, name = "n_cells")

  list(obj = obj, summary = summary)
}

markers <- markers[markers$specificity > 0.2 & markers$cell_type != "Malignant", ]
markers_list <- markers %>%
  mutate(cell_type = recode(cell_type, !!!celltype_map)) %>%
  split(.$cell_type)

markers_ranked <- lapply(markers_list, function(df) {
  combine_marker_scores(df, w_specificity = 0.2, w_sensitivity = 0.8)
})

N <- 100
setsN <- markers_ranked |>
  imap(~ .x %>% arrange(desc(Combined)) %>% slice_head(n = N) %>% pull(gene))

celltypes <- names(setsN)
all_genes <- unique(unlist(setsN))
valid_seed_celltypes <- intersect(celltypes, names(canonical_markers))
seed_summary_list <- list()

for (sample in sample_dirs) {
  seed_res <- build_initial_seed(
    obj = tmdata_list[[sample]],
    sample_id = sample,
    valid_celltypes = valid_seed_celltypes,
    expr_cut = 1,
    min_gap = 0.1
  )
  tmdata_list[[sample]] <- seed_res$obj
  seed_summary_list[[sample]] <- seed_res$summary
}

seed_summary <- bind_rows(seed_summary_list)
write.csv(seed_summary, "celltype_initial_seed_summary.csv", row.names = FALSE)

seed_retention <- seed_summary %>%
  group_by(sample, celltype_initial) %>%
  summarise(
    n_total_exact = sum(n_cells),
    n_seed = sum(n_cells[celltype_initial_seed == celltype_initial], na.rm = TRUE),
    pct_seed = ifelse(n_total_exact > 0, 100 * n_seed / n_total_exact, 0),
    .groups = "drop"
  )
write.csv(seed_retention, "celltype_initial_seed_retention.csv", row.names = FALSE)

method <- "celltype_initial_seed"

get_excluded_off_cts <- function(home_ct) {
  excluded <- character(0)

  if (grepl("^b.cell", home_ct, ignore.case = TRUE)) {
    excluded <- c(excluded, "plasma", "dendritic")
  }
  if (grepl("^plasma", home_ct, ignore.case = TRUE)) {
    excluded <- c(excluded, "b.cell", "dendritic")
  }
  if (grepl("^t.cell", home_ct, ignore.case = TRUE)) {
    excluded <- c(excluded, "nk.cell", "dendritic")
  }
  if (grepl("^nk.cell", home_ct, ignore.case = TRUE)) {
    excluded <- c(excluded, "t.cell", "dendritic")
  }
  if (grepl("^macrophage", home_ct, ignore.case = TRUE)) {
    excluded <- c(excluded, "dendritic")
  }

  excluded
}

pct_tbl <- imap_dfr(tmdata_list, function(obj, sid) {
  meta <- obj@meta.data
  genes <- intersect(all_genes, rownames(obj))
  expr <- GetAssayData(obj, slot = "data")[genes, , drop = FALSE]
  meta_ct <- meta %>% pull(!!sym(method))
  split_cells <- split(colnames(obj), meta_ct)

  map_dfr(names(split_cells), function(ct) {
    cells <- split_cells[[ct]]
    if (length(cells) == 0) {
      return(NULL)
    }
    pct <- Matrix::rowMeans(expr[, cells, drop = FALSE] > 1) * 100
    tibble(sample = sid, celltype = ct, gene = genes, pct_expr = pct)
  })
})

pct_tbl <- pct_tbl %>% filter(celltype %in% celltypes)

markers_long <- tibble(celltype = names(setsN), gene = setsN) %>% unnest(gene)
marker_rank_tbl <- imap_dfr(markers_ranked, function(df, ct) {
  tibble(
    celltype = ct,
    gene = df$gene,
    combined = df$Combined
  )
})

expr_off_threshold <- 15
off_max_pct_samps <- 30

derive_exclusive_tbl <- function(expr_home_threshold, home_min_pct_samps, tier) {
  pct_with_homeflag <- pct_tbl %>%
    left_join(markers_long %>% mutate(is_home = TRUE), by = c("gene", "celltype")) %>%
    mutate(
      is_home = coalesce(is_home, FALSE),
      pass = ifelse(is_home, pct_expr > expr_home_threshold, pct_expr > expr_off_threshold)
    )

  gene_ct_summary <- pct_with_homeflag %>%
    group_by(gene, celltype, is_home) %>%
    summarise(
      n_total = sum(!is.na(pct_expr)),
      n_pass = sum(pass, na.rm = TRUE),
      pct_samples = ifelse(n_total > 0, 100 * n_pass / n_total, 0),
      .groups = "drop"
    )

  gene_ct_long <- gene_ct_summary %>%
    select(gene, celltype, pct_samples, is_home)

  home_tbl <- markers_long %>%
    left_join(gene_ct_long, by = c("gene", "celltype")) %>%
    mutate(home_pct = coalesce(pct_samples, 0)) %>%
    select(gene, celltype, home_pct)

  off_tbl_raw <- gene_ct_long %>%
    rename(off_celltype = celltype, off_pct = pct_samples) %>%
    inner_join(markers_long, by = "gene")

  off_tbl_filtered <- off_tbl_raw %>%
    rowwise() %>%
    mutate(
      drop_these = list(get_excluded_off_cts(celltype)),
      keep_row = !(off_celltype %in% drop_these)
    ) %>%
    ungroup() %>%
    filter(off_celltype != celltype, keep_row) %>%
    select(-drop_these, -keep_row)

  off_tbl <- off_tbl_filtered %>%
    group_by(gene, celltype) %>%
    summarise(off_pct_max = max(off_pct, na.rm = TRUE), .groups = "drop") %>%
    mutate(off_pct_max = ifelse(is.finite(off_pct_max), off_pct_max, 0))

  home_tbl %>%
    left_join(marker_rank_tbl, by = c("celltype", "gene")) %>%
    left_join(off_tbl, by = c("gene", "celltype")) %>%
    mutate(
      tier = tier,
      expr_home_threshold = expr_home_threshold,
      expr_off_threshold = expr_off_threshold,
      home_min_pct_samps = home_min_pct_samps,
      off_max_pct_samps = off_max_pct_samps,
      off_pct_max = coalesce(off_pct_max, 0),
      combined = coalesce(combined, -Inf),
      is_exclusive = (home_pct > home_min_pct_samps) & (off_pct_max < off_max_pct_samps)
    ) %>%
    arrange(desc(is_exclusive), desc(home_pct), off_pct_max, desc(combined))
}

relax_grid <- tibble(
  tier = c("strict", "home_expr25", "home_expr20_samp10", "home_expr15_samp5", "home_expr10_samp1"),
  expr_home_threshold = c(30, 25, 20, 15, 10),
  home_min_pct_samps = c(15, 15, 10, 5, 1)
)

exclusive_tbls <- pmap(
  relax_grid,
  \(tier, expr_home_threshold, home_min_pct_samps) {
    derive_exclusive_tbl(
      expr_home_threshold = expr_home_threshold,
      home_min_pct_samps = home_min_pct_samps,
      tier = tier
    )
  }
)
names(exclusive_tbls) <- relax_grid$tier

exclusive_tbls_long <- bind_rows(exclusive_tbls, .id = "tier_id")
saveRDS(exclusive_tbls, "exclusive_tbls.rds")
saveRDS(exclusive_tbls_long, "exclusive_tbls_long.rds")
write.csv(exclusive_tbls_long, "exclusive_tbls.csv", row.names = FALSE)

save <- exclusive_tbls[[1]] %>%
  filter(is_exclusive) %>%
  select(gene, home_celltype = celltype, home_pct, off_pct_max)

required_marker_celltypes <- setdiff(celltypes, "dendritic")
target_marker_count <- 4L

for (ct in required_marker_celltypes) {
  current_genes <- unique(save$gene[save$home_celltype == ct])
  if (length(current_genes) >= target_marker_count) {
    next
  }

  for (tier_name in names(exclusive_tbls)[-1]) {
    needed <- target_marker_count - length(current_genes)
    if (needed <= 0) {
      break
    }

    additions <- exclusive_tbls[[tier_name]] %>%
      filter(is_exclusive, celltype == ct, !gene %in% current_genes) %>%
      arrange(desc(home_pct), off_pct_max, desc(combined)) %>%
      transmute(gene, home_celltype = celltype, home_pct, off_pct_max) %>%
      slice_head(n = needed)

    if (nrow(additions) > 0) {
      save <- bind_rows(save, additions)
      current_genes <- unique(c(current_genes, additions$gene))
    }
  }
}

save <- save %>%
  distinct(home_celltype, gene, .keep_all = TRUE) %>%
  arrange(home_celltype, desc(home_pct), off_pct_max)

missing_required_celltypes <- setdiff(required_marker_celltypes, unique(save$home_celltype))
if (length(missing_required_celltypes) > 0) {
  stop(
    "Required annotation marker celltypes still missing after home-side relaxation with strict exclusivity: ",
    paste(missing_required_celltypes, collapse = ", ")
  )
}

if (nrow(save) == 0) {
  stop("No exclusive annotation markers were derived.")
}

saveRDS(save, "anno_markers.rds")
markers_list <- split(save$gene, save$home_celltype)
score_weight_floor <- 0.1
score_marker_tbl <- save %>%
  mutate(score_weight = pmax(home_pct - off_pct_max, score_weight_floor))
write.csv(score_marker_tbl, "annotation_score_markers.csv", row.names = FALSE)
score_marker_tbl_top4 <- score_marker_tbl %>%
  group_by(home_celltype) %>%
  arrange(desc(home_pct), off_pct_max, .by_group = TRUE) %>%
  slice_head(n = 4) %>%
  ungroup()
write.csv(score_marker_tbl_top4, "annotation_score_markers_top4.csv", row.names = FALSE)

score_markers_list <- split(score_marker_tbl$gene, score_marker_tbl$home_celltype)
score_weights_list <- split(score_marker_tbl$score_weight, score_marker_tbl$home_celltype)

gap_cut <- 0.8
allowed_pairs <- list(
  c("t.cell", "nk.cell"),
  c("nk.cell", "t.cell"),
  c("b.cell", "plasma"),
  c("plasma", "b.cell")
)

normalize_annotation_label <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- NA_character_
  x <- trimws(x)
  x[grepl("^unknown(\\s*<>\\s*unknown)?$", x)] <- "unknown"
  x
}

is_resolved_annotation <- function(x) {
  !is.na(x) &&
    nzchar(x) &&
    !(x %in% c("unknown", "unresolved"))
}

is_single_annotation <- function(x) {
  is_resolved_annotation(x) &&
    !grepl("\\|", x) &&
    !grepl("<>", x, fixed = TRUE)
}

split_annotation_tokens <- function(x) {
  if (is.na(x) || !nzchar(x)) {
    return(character(0))
  }
  raw <- gsub("\\s*<>\\s*", "|", x)
  toks <- trimws(unlist(strsplit(raw, "\\|", fixed = FALSE), use.names = FALSE))
  toks[nzchar(toks)]
}

label_contains_single <- function(label, single_label) {
  if (!is_single_annotation(single_label)) {
    return(FALSE)
  }
  single_label %in% split_annotation_tokens(label)
}

combine_annotation_calls <- function(initial_label, current_label) {
  s1 <- normalize_annotation_label(initial_label)
  s2 <- normalize_annotation_label(current_label)

  if (!is.na(s1) && !is.na(s2) && identical(s1, s2)) {
    return(s1)
  }

  if (is_single_annotation(s1) && (is.na(s2) || s2 == "unresolved")) {
    return(s1)
  }

  if (is_single_annotation(s2) && (is.na(s1) || s1 %in% c("unknown", "unresolved"))) {
    return(s2)
  }

  if (is_single_annotation(s1) && label_contains_single(s2, s1)) {
    return(s1)
  }

  if (is_single_annotation(s2) && label_contains_single(s1, s2)) {
    return(s2)
  }

  if (!is.na(s1) && s1 == "unknown" && (is.na(s2) || s2 == "unresolved")) {
    return("unresolved")
  }

  if (!is.na(s2) && s2 == "unresolved" && !is.na(s1) && !is_single_annotation(s1)) {
    return("unresolved")
  }

  if (is.na(s1) || !nzchar(s1)) {
    s1 <- "NA"
  }
  if (is.na(s2) || !nzchar(s2)) {
    s2 <- "NA"
  }

  paste0(s1, " <> ", s2)
}

call_clusters_from_scores <- function(scores_long) {
  scores_long %>%
    group_by(cluster) %>%
    summarize(
      cell_types = list(cell_type),
      zs = list(z),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      call = {
        az <- unlist(zs)
        act <- unlist(cell_types)

        if (length(az) == 0 || all(is.na(az))) {
          "unresolved"
        } else {
          az_safe <- ifelse(is.na(az), -Inf, az)
          max_idx <- which.max(az_safe)
          top_ct <- act[max_idx]
          top_z <- az[max_idx]
          other_idx <- setdiff(seq_along(az), max_idx)

          if (length(other_idx) == 0) {
            top_ct
          } else {
            margins <- top_z - az[other_idx]
            close_idx <- other_idx[margins < gap_cut]

            if (length(close_idx) == 0) {
              top_ct
            } else if (length(close_idx) == 1) {
              other_ct <- act[close_idx]
              pair_vec <- c(top_ct, other_ct)
              is_allowed <- any(vapply(allowed_pairs, function(p) identical(p, pair_vec), logical(1)))
              if (is_allowed) {
                paste0(top_ct, "|", other_ct)
              } else {
                "unresolved"
              }
            } else {
              "unresolved"
            }
          }
        }
      }
    ) %>%
    ungroup() %>%
    select(cluster, call)
}

build_weighted_score_long <- function(tmdata, score_markers_in_data, score_weights_in_data, agg_fun = c("median", "mean")) {
  agg_fun <- match.arg(agg_fun)

  score_gene_union <- unique(unlist(score_markers_in_data))
  expr_mat <- GetAssayData(tmdata, slot = "data")[score_gene_union, , drop = FALSE]
  expr_mat <- as.matrix(expr_mat)
  z_mat <- t(scale(t(expr_mat)))
  z_mat[!is.finite(z_mat)] <- 0

  ct_names <- names(score_markers_in_data)
  score_mat <- matrix(
    NA_real_,
    nrow = ncol(tmdata),
    ncol = length(ct_names),
    dimnames = list(colnames(tmdata), ct_names)
  )

  for (ct in ct_names) {
    genes <- score_markers_in_data[[ct]]
    weights <- as.numeric(score_weights_in_data[[ct]])
    weights <- weights / sum(weights)
    score_mat[, ct] <- as.numeric(crossprod(weights, z_mat[genes, , drop = FALSE]))
  }

  for (ct in ct_names) {
    tmdata@meta.data[[paste0("mod_", ct)]] <- score_mat[colnames(tmdata), ct]
  }

  mod_cols <- setNames(paste0("mod_", ct_names), ct_names)
  cluster_score_tbl <- tmdata@meta.data %>%
    mutate(cluster = as.character(tmdata$seurat_clusters)) %>%
    group_by(cluster) %>%
    summarize(
      across(
        all_of(mod_cols),
        \(x) if (agg_fun == "median") median(x, na.rm = TRUE) else mean(x, na.rm = TRUE)
      ),
      .groups = "drop"
    )

  cluster_score_tbl %>%
    pivot_longer(-cluster, names_to = "mod", values_to = "score") %>%
    mutate(cell_type = names(mod_cols)[match(mod, names(mod_cols))]) %>%
    group_by(cluster) %>%
    mutate(z = as.numeric(scale(score))) %>%
    ungroup() %>%
    select(cluster, cell_type, score, z)
}

build_canonical_score_long <- function(tmdata, expr_cut = 1) {
  expr <- GetAssayData(tmdata, slot = "data")
  clusters <- sort(unique(as.character(tmdata$seurat_clusters)))

  map_dfr(clusters, function(cl) {
    cells <- colnames(tmdata)[tmdata$seurat_clusters == cl]

    map_dfr(names(canonical_markers), function(ct) {
      genes <- intersect(canonical_markers[[ct]], rownames(tmdata))
      if (length(genes) == 0) {
        return(NULL)
      }

      gene_pct <- Matrix::rowMeans(expr[genes, cells, drop = FALSE] > expr_cut)
      tibble(
        cluster = cl,
        cell_type = ct,
        score = mean(gene_pct),
        ngenes_gt10 = sum(gene_pct > 0.10)
      )
    })
  })
}

call_canonical_rescue <- function(canonical_scores, cluster_sizes, min_cluster_size = 50L) {
  canonical_scores %>%
    arrange(cluster, desc(score), cell_type) %>%
    group_by(cluster) %>%
    summarize(
      cell_types = list(cell_type),
      scores = list(score),
      ngenes = list(ngenes_gt10),
      n_cells = cluster_sizes[match(cluster[1], names(cluster_sizes))],
      call = {
        act <- unlist(cell_types)
        az <- unlist(scores)
        ng <- unlist(ngenes)

        if (length(az) == 0 || is.na(n_cells) || n_cells < min_cluster_size) {
          "unresolved"
        } else {
          ord <- order(az, decreasing = TRUE)
          top_idx <- ord[1]
          second_idx <- if (length(ord) >= 2) ord[2] else NA_integer_

          top_ct <- act[top_idx]
          top_score <- az[top_idx]
          second_score <- if (is.na(second_idx)) -Inf else az[second_idx]
          top_ngenes <- ng[top_idx]

          is_immune <- top_ct %in% c("b.cell", "t.cell", "nk.cell", "plasma")
          top_thresh <- if (is_immune) 0.20 else 0.12
          second_thresh <- if (is_immune) 0.05 else 0.08
          ngene_thresh <- if (is_immune) 3L else 2L

          if (
            is.finite(top_score) &&
            top_score >= top_thresh &&
            second_score <= second_thresh &&
            top_ngenes >= ngene_thresh
          ) {
            top_ct
          } else {
            "unresolved"
          }
        }
      },
      .groups = "drop"
    ) %>%
    select(cluster, call)
}

choose_refined_call <- function(initial_label, median_call, mean_call, canonical_call, cluster_n, min_cluster_size = 50L) {
  initial_label <- normalize_annotation_label(initial_label)
  median_call <- normalize_annotation_label(median_call)
  mean_call <- normalize_annotation_label(mean_call)
  canonical_call <- normalize_annotation_label(canonical_call)

  if (is_resolved_annotation(median_call) && is_resolved_annotation(mean_call)) {
    if (identical(median_call, mean_call)) {
      return(median_call)
    }
    if (is_single_annotation(initial_label) && identical(median_call, initial_label)) {
      return(median_call)
    }
    if (is_single_annotation(initial_label) && identical(mean_call, initial_label)) {
      return(mean_call)
    }
    if (is_single_annotation(median_call) && grepl("\\|", mean_call)) {
      return(median_call)
    }
    if (grepl("\\|", median_call) && is_single_annotation(mean_call)) {
      return(mean_call)
    }
    return(median_call)
  }

  if (is_resolved_annotation(median_call)) {
    return(median_call)
  }

  if (is_resolved_annotation(mean_call)) {
    if (!is.na(cluster_n) && cluster_n >= min_cluster_size) {
      if (is_single_annotation(initial_label) && identical(mean_call, initial_label)) {
        return(mean_call)
      }
      if (is_resolved_annotation(canonical_call) && identical(mean_call, canonical_call)) {
        return(mean_call)
      }
    }
    return("unresolved")
  }

  "unresolved"
}

choose_final_call <- function(initial_label, refined_call, canonical_call, cluster_n, min_cluster_size = 50L) {
  initial_label <- normalize_annotation_label(initial_label)
  refined_call <- normalize_annotation_label(refined_call)
  canonical_call <- normalize_annotation_label(canonical_call)

  current_call <- refined_call
  if (
    !is_resolved_annotation(current_call) &&
    !is.na(cluster_n) &&
    cluster_n >= min_cluster_size &&
    is_resolved_annotation(canonical_call)
  ) {
    current_call <- canonical_call
  }

  if (is_resolved_annotation(current_call)) {
    return(current_call)
  }

  if (is_single_annotation(initial_label)) {
    return(initial_label)
  }

  "unresolved"
}

pct_list <- list()
count_list <- list()

for (sample in sample_dirs) {
  tmdata <- tmdata_list[[sample]]
  Idents(tmdata) <- tmdata$seurat_clusters
  available_genes <- rownames(tmdata)

  score_markers_in_data <- imap(score_markers_list, function(gene_set, ct) {
    intersect(gene_set, available_genes)
  })
  score_markers_in_data <- Filter(function(v) length(v) > 0, score_markers_in_data)
  ct_names <- names(score_markers_in_data)

  if (length(ct_names) == 0) {
    stop("No robust annotation markers found in sample: ", sample)
  }

  weight_lookup <- imap(score_weights_list, function(w, ct) {
    names(w) <- score_markers_list[[ct]]
    w
  })
  score_weights_in_data <- imap(score_markers_in_data, function(genes, ct) {
    weight_lookup[[ct]][genes]
  })

  scores_long_median <- build_weighted_score_long(
    tmdata = tmdata,
    score_markers_in_data = score_markers_in_data,
    score_weights_in_data = score_weights_in_data,
    agg_fun = "median"
  )
  scores_long_mean <- build_weighted_score_long(
    tmdata = tmdata,
    score_markers_in_data = score_markers_in_data,
    score_weights_in_data = score_weights_in_data,
    agg_fun = "mean"
  )
  canonical_scores <- build_canonical_score_long(tmdata, expr_cut = 1)

  write.csv(
    scores_long_median,
    file.path("by_samples", sample, paste0(sample, "_cluster_scores.csv")),
    row.names = FALSE
  )
  write.csv(
    scores_long_median,
    file.path("by_samples", sample, paste0(sample, "_cluster_scores_median.csv")),
    row.names = FALSE
  )
  write.csv(
    scores_long_mean,
    file.path("by_samples", sample, paste0(sample, "_cluster_scores_mean.csv")),
    row.names = FALSE
  )
  write.csv(
    canonical_scores,
    file.path("by_samples", sample, paste0(sample, "_cluster_canonical_scores.csv")),
    row.names = FALSE
  )

  median_calls <- call_clusters_from_scores(scores_long_median) %>%
    rename(celltype_refined_median = call)
  mean_calls <- call_clusters_from_scores(scores_long_mean) %>%
    rename(celltype_refined_mean = call)
  cluster_sizes <- table(as.character(tmdata$seurat_clusters))
  canonical_calls <- call_canonical_rescue(canonical_scores, cluster_sizes = cluster_sizes) %>%
    rename(celltype_canonical = call)

  cluster_initial <- tmdata@meta.data %>%
    mutate(cluster = as.character(tmdata$seurat_clusters)) %>%
    count(cluster, celltype_initial, name = "n_cells") %>%
    group_by(cluster) %>%
    slice_max(order_by = n_cells, n = 1, with_ties = FALSE) %>%
    ungroup()

  cluster_summary <- cluster_initial %>%
    left_join(median_calls, by = "cluster") %>%
    left_join(mean_calls, by = "cluster") %>%
    left_join(canonical_calls, by = "cluster") %>%
    rowwise() %>%
    mutate(
      celltype_refined = choose_refined_call(
        initial_label = celltype_initial,
        median_call = celltype_refined_median,
        mean_call = celltype_refined_mean,
        canonical_call = celltype_canonical,
        cluster_n = n_cells
      ),
      celltype_update = choose_final_call(
        initial_label = celltype_initial,
        refined_call = celltype_refined,
        canonical_call = celltype_canonical,
        cluster_n = n_cells
      ),
      celltype_update_detail = {
        detail_source <- if (is_resolved_annotation(celltype_refined)) {
          celltype_refined
        } else if (n_cells >= 50 && is_resolved_annotation(celltype_canonical)) {
          celltype_canonical
        } else {
          "unresolved"
        }
        combine_annotation_calls(celltype_initial, detail_source)
      }
    ) %>%
    ungroup()

  write.csv(
    cluster_summary,
    file.path("by_samples", sample, paste0(sample, "_cluster_annotation_summary.csv")),
    row.names = FALSE
  )

  cluster_maps <- list(
    celltype_refined_median = setNames(cluster_summary$celltype_refined_median, cluster_summary$cluster),
    celltype_refined_mean = setNames(cluster_summary$celltype_refined_mean, cluster_summary$cluster),
    celltype_canonical = setNames(cluster_summary$celltype_canonical, cluster_summary$cluster),
    celltype_refined = setNames(cluster_summary$celltype_refined, cluster_summary$cluster),
    celltype_update = setNames(cluster_summary$celltype_update, cluster_summary$cluster),
    celltype_update_detail = setNames(cluster_summary$celltype_update_detail, cluster_summary$cluster)
  )

  for (meta_col in names(cluster_maps)) {
    tmdata[[meta_col]] <- as.character(cluster_maps[[meta_col]][as.character(tmdata$seurat_clusters)])
  }

  saveRDS(tmdata, file.path("by_samples", sample, paste0(sample, "_anno.rds")), compress = FALSE)

  if ("umap" %in% names(tmdata@reductions)) {
    p1 <- DimPlot(tmdata, group.by = "seurat_clusters", label = TRUE) + ggtitle("Clusters")
    p2 <- DimPlot(tmdata, group.by = "celltype_initial", label = FALSE) + ggtitle("Celltype Initial")
    p3 <- DimPlot(tmdata, group.by = "celltype_update", label = FALSE) + ggtitle("Celltype Update")
    ggsave(
      filename = file.path("by_samples", sample, paste0(sample, "_cluster_overview.png")),
      plot = p1 + p2 + p3,
      width = 15,
      height = 4,
      dpi = 300
    )
  }

  tab <- table(tmdata$celltype_update)
  pct <- 100 * tab / sum(tab)
  pct_list[[sample]] <- data.frame(
    study = sample,
    celltype = names(pct),
    pct = as.numeric(pct),
    stringsAsFactors = FALSE
  )
  count_list[[sample]] <- data.frame(
    study = sample,
    celltype = names(tab),
    count = as.integer(tab),
    stringsAsFactors = FALSE
  )

  message("annotated ", sample)
}

pct_all <- bind_rows(pct_list)
count_all <- bind_rows(count_list)
pct_mat <- pct_all %>%
  pivot_wider(names_from = celltype, values_from = pct, values_fill = 0)
count_mat <- count_all %>%
  pivot_wider(names_from = celltype, values_from = count, values_fill = 0)

saveRDS(pct_mat, file = "pct_mat.rds")
saveRDS(count_mat, file = "count_mat.rds")
write.csv(pct_mat, file = "pct_mat.csv", row.names = FALSE)
write.csv(count_mat, file = "count_mat.csv", row.names = FALSE)
