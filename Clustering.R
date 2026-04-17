suppressPackageStartupMessages({
  library("Seurat")
  library("dplyr")
  library("purrr")
  library("tidyr")
  library("ggplot2")
  library("readxl")
})

reticulate::use_condaenv(
  "dmtcp",
  conda = "/rds/general/user/sg3723/home/anaconda3/bin/conda",
  required = TRUE
)

set.seed(1)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs")

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
message(sample)

sample_path <- file.path("by_samples", sample, paste0(sample, ".rds"))
if (!file.exists(sample_path)) {
  stop("Missing sample object: ", sample_path)
}

tmdata <- readRDS(sample_path)
ncells <- ncol(tmdata)

run_pca_safe <- function(obj, npcs) {
  tryCatch(
    RunPCA(obj, npcs = npcs, verbose = FALSE),
    error = function(e) {
      if (grepl("starting vector near the null space", conditionMessage(e), fixed = TRUE)) {
        RunPCA(obj, npcs = npcs, approx = FALSE, verbose = FALSE)
      } else {
        stop(e)
      }
    }
  )
}

if (ncells > 50) {
  tmdata <- FindVariableFeatures(tmdata, nfeatures = 3000)
  tmdata <- ScaleData(tmdata, verbose = FALSE)
  pca_npcs <- min(50, max(2, ncells - 1))
  tmdata <- run_pca_safe(tmdata, npcs = pca_npcs)

  if (ncells > 20000) {
    k_use <- 10
    pc_use <- 1:min(20, pca_npcs)
  } else {
    k_use <- min(30, max(10, round(sqrt(ncells) - 2)))
    pc_use <- 1:min(30, pca_npcs)
  }
  tmdata <- FindNeighbors(tmdata, dims = pc_use, k.param = k_use, verbose = FALSE)
  tmdata <- FindClusters(
    tmdata,
    resolution = 0.8,
    algorithm = 4,
    method = "igraph",
    verbose = FALSE
  )
  umap_neighbors <- if (ncells > 20000) 15 else 30
  tmdata <- RunUMAP(tmdata, dims = pc_use, n.neighbors = umap_neighbors, min.dist = 0.3, verbose = FALSE)
} else {
  tmdata <- FindVariableFeatures(tmdata, nfeatures = 3000, verbose = FALSE)
  tmdata <- ScaleData(tmdata, features = VariableFeatures(tmdata), verbose = FALSE)

  max_pcs <- min(length(VariableFeatures(tmdata)), ncells - 1)
  npcs <- max(2, max_pcs)
  tmdata <- run_pca_safe(tmdata, npcs = npcs)
  pc_use <- 1:npcs
  k_use <- max(5, min(15, ncells - 1))
  tmdata <- FindNeighbors(tmdata, dims = pc_use, k.param = k_use, verbose = FALSE)
  tmdata <- FindClusters(
    tmdata,
    resolution = 0.8,
    algorithm = 4,
    method = "igraph",
    verbose = FALSE
  )

  umap_neighbors <- max(5, min(15, ncells - 1))
  tmdata <- RunUMAP(tmdata, dims = pc_use, n.neighbors = umap_neighbors, min.dist = 0.3, verbose = FALSE)
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

combine_marker_scores <- function(df, w_specificity = 0.2, w_sensitivity = 0.8) {
  pr <- function(x) {
    r <- rank(x, ties.method = "average", na.last = "keep")
    r / (sum(!is.na(x)) + 1)
  }

  combined <- (w_specificity * pr(df$specificity) + w_sensitivity * pr(df$sensitivity)) /
    (w_specificity + w_sensitivity)

  df %>% mutate(Combined = combined) %>% arrange(desc(Combined))
}

markers <- markers[markers$specificity > 0.2 & markers$cell_type != "Malignant", ]
markers_list <- markers %>%
  mutate(cell_type = recode(cell_type, !!!celltype_map)) %>%
  split(.$cell_type)

markers_ranked <- lapply(markers_list, function(df) {
  combine_marker_scores(df, w_specificity = 0.2, w_sensitivity = 0.8)
})

N <- 50
lfc_th <- 0.5
z_cut <- 1
z_gap <- 0.5
diffone <- 3

Idents(tmdata) <- tmdata$seurat_clusters

setsN <- markers_ranked |>
  imap(~ .x %>% arrange(desc(Combined)) %>% slice_head(n = N) %>% pull(gene) %>% intersect(rownames(tmdata)))

gap_signature_call <- function(cell_types, zs, z_cut, z_gap) {
  az <- as.numeric(unlist(zs))
  act <- as.character(unlist(cell_types))

  if (length(az) == 0 || all(is.na(az))) {
    return("unknown")
  }

  valid_idx <- which(!is.na(az))
  if (length(valid_idx) == 0) {
    return("unknown")
  }

  ord <- valid_idx[order(az[valid_idx], decreasing = TRUE)]
  max_idx <- ord[1]
  top_ct <- act[max_idx]
  top_z <- az[max_idx]

  if (!is.finite(top_z) || top_z < z_cut) {
    return("unknown")
  }

  if (length(ord) == 1) {
    return(top_ct)
  }

  next_z <- az[ord[2]]
  if (!is.finite(next_z)) {
    return(top_ct)
  }

  if ((top_z - next_z) >= z_gap) {
    return(top_ct)
  }

  "unknown"
}

if (ncells > 50) {
  de <- FindAllMarkers(tmdata, only.pos = TRUE, min.pct = 0.1, logfc.threshold = lfc_th)
  if (!"pct.1" %in% colnames(de)) {
    stop("FindAllMarkers output is missing pct.1; cannot apply cluster-detection filter.")
  }
  de <- de %>% filter(pct.1 >= 0.5)
  universe <- rownames(tmdata)
  clusters_all <- sort(unique(as.character(tmdata$seurat_clusters)))

  enrich_one <- function(cl) {
    de_genes <- de %>% filter(cluster == cl) %>% pull(gene) %>% unique()
    if (length(de_genes) == 0) {
      return(tibble(cluster = cl, step1 = "unknown"))
    }

    res <- imap_dfr(setsN, ~ {
      a <- length(intersect(.x, de_genes))
      b <- length(setdiff(de_genes, .x))
      c <- length(setdiff(.x, de_genes))
      d <- length(universe) - a - b - c
      p <- fisher.test(matrix(c(a, b, c, d), 2, 2), alternative = "greater")$p.value
      tibble(cell_type = .y, overlap = a, pval = p)
    }) %>%
      mutate(padj = p.adjust(pval, "BH"), cluster = cl)

    sig <- res %>%
      filter(padj <= 0.05) %>%
      arrange(desc(overlap), padj)

    if (nrow(sig) == 0) {
      return(tibble(cluster = cl, step1 = "unknown"))
    }

    top_ov <- sig$overlap[1]
    if (nrow(sig) == 1) {
      return(tibble(cluster = cl, step1 = sig$cell_type[1]))
    }

    second_ov <- sig$overlap[2]
    if ((top_ov - second_ov) >= diffone) {
      return(tibble(cluster = cl, step1 = sig$cell_type[1]))
    }

    keep <- sig %>%
      filter((top_ov - overlap) <= diffone, overlap > diffone)

    if (nrow(keep) == 0) {
      return(tibble(cluster = cl, step1 = "unknown"))
    }

    calls <- keep %>%
      arrange(desc(overlap), padj) %>%
      pull(cell_type) %>%
      unique()

    tibble(cluster = cl, step1 = paste(calls, collapse = "|"))
  }

  if (nrow(de) == 0 || !"cluster" %in% colnames(de)) {
    step1_calls <- tibble(cluster = clusters_all, step1 = "unknown")
  } else {
    step1_calls <- clusters_all %>% map_dfr(enrich_one)
  }

  ct_names <- names(setsN)
  tmdata <- AddModuleScore(tmdata, features = unname(setsN), name = "mod_")
  mod_cols <- paste0("mod_", seq_along(setsN))
  names(mod_cols) <- ct_names

  scores_long <- tmdata@meta.data %>%
    mutate(cluster = tmdata$seurat_clusters) %>%
    group_by(cluster) %>%
    summarize(across(all_of(mod_cols), \(x) median(x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(-cluster, names_to = "mod", values_to = "score") %>%
    mutate(cell_type = names(mod_cols)[match(mod, names(mod_cols))]) %>%
    group_by(cluster) %>%
    mutate(z = as.numeric(scale(score))) %>%
    ungroup() %>%
    select(cluster, cell_type, score, z)

  step2_calls <- scores_long %>%
    group_by(cluster) %>%
    summarize(
      cell_types = list(cell_type),
      zs = list(z),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      step2 = gap_signature_call(cell_types, zs, z_cut = z_cut, z_gap = z_gap)
    ) %>%
    ungroup() %>%
    select(cluster, step2)

  final_calls <- step1_calls %>%
    left_join(step2_calls, by = "cluster") %>%
    rowwise() %>%
    mutate(
      final = {
        s1 <- step1
        s2 <- step2
        is_single <- function(x) !is.na(x) && x != "unknown" && !grepl("\\|", x)

        if (is_single(s1) && is_single(s2) && s1 == s2) {
          s1
        } else if (is_single(s1) && (is.na(s2) || s2 == "unknown")) {
          s1
        } else if (is_single(s2) && (is.na(s1) || s1 == "unknown")) {
          s2
        } else {
          paste0(s1, " <> ", s2)
        }
      }
    ) %>%
    ungroup()

  label_map <- final_calls$final
  names(label_map) <- final_calls$cluster
  tmdata$celltype_initial <- as.character(label_map[as.character(tmdata$seurat_clusters)])
  tmdata$celltype_initial[is.na(tmdata$celltype_initial)] <- "unknown"
} else {
  tmdata$celltype_initial <- rep("unknown", ncol(tmdata))
}

saveRDS(tmdata, sample_path, compress = FALSE)

if ("umap" %in% names(tmdata@reductions)) {
  p1 <- DimPlot(tmdata, group.by = "seurat_clusters", label = TRUE) + ggtitle("Clusters")
  p2 <- DimPlot(tmdata, group.by = "celltype_initial", label = FALSE) + ggtitle("Celltype Initial")
  ggsave(
    filename = file.path("by_samples", sample, paste0(sample, "_cluster_overview.png")),
    plot = p1 + p2,
    width = 10,
    height = 4,
    dpi = 300
  )
}
