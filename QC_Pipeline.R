suppressPackageStartupMessages({
  library("Seurat")
  library("dplyr")
  library("ggplot2")
  library("grid")
  library("gridExtra")
  library("Matrix")
  library("tibble")
  library("tidyr")
})

set.seed(1)

output_root <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs"
source_root <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples"

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
setwd(output_root)

batch <- "_sn"
names_tmdata <- readLines(file.path(source_root, "names_tmdata_sn.txt"))
max_mt <- 15
min_ngenes <- 200
min_hk_expr <- 0.5

fibroblast <- c("COL3A1", "COL1A1", "COL1A2", "LUM")
macrophage <- c("CSF1R", "TYROBP", "CD14", "CD163", "AIF1")
mast <- c("MS4A2", "CPA3", "TPSB2")
lymph <- c("CCL21")
epithelial <- c("KRT7", "MUC1", "KRT19", "EPCAM")
t.cell <- c("CD3E", "CD3D", "CD2")
b.cell <- c("MS4A1", "CD79B", "CD79A")
dendritic <- c("CLEC10A", "CCR7", "CD86")
endothelial <- c("ENG", "CLEC14A", "CLDN5", "VWF", "CDH5")
erythrocyte <- c("HBA1", "HBA2", "HBB")
keratinocyte <- c("FLG", "IVL")
neutrophil <- c("CTSG", "ELANE", "MPO", "AZU1")
housekeeping <- c(
  "ACTB", "GAPDH", "RPS11", "RPS13", "RPS14", "RPS15", "RPS16", "RPS18",
  "RPS19", "RPS20", "RPL10", "RPL13", "RPL15", "RPL18"
)

markers_list <- list(
  b.cell = b.cell,
  dendritic = dendritic,
  endothelial = endothelial,
  epithelial = epithelial,
  erythrocyte = erythrocyte,
  fibroblast = fibroblast,
  keratinocyte = keratinocyte,
  lymph = lymph,
  macrophage = macrophage,
  mast = mast,
  neutrophil = neutrophil,
  t.cell = t.cell
)

classify_technology <- function(source_batch) {
  dplyr::case_when(
    source_batch == "gemx" ~ "GEMX",
    source_batch == "cynthia_sn" ~ "snSeq",
    TRUE ~ "Multiome"
  )
}

read_sample <- function(path) {
  sample_id <- basename(path)
  source_batch <- dirname(path)
  raw_counts <- Read10X(file.path(source_root, path))
  if (is.list(raw_counts)) {
    if ("Gene Expression" %in% names(raw_counts)) {
      raw_counts <- raw_counts[["Gene Expression"]]
    } else {
      raw_counts <- raw_counts[[1]]
    }
  }

  tmdata <- CreateSeuratObject(counts = raw_counts)
  tmdata$orig.ident <- sample_id
  tmdata$input_source <- source_batch
  tmdata$reference_batch <- source_batch
  tmdata$technology <- classify_technology(source_batch)
  tmdata$sample_id <- sample_id
  tmdata
}

####################
# scRef-style QC plotting helpers
####################
build_inspection_plots <- function(tmdata, sample_id) {
  plot_data <- tmdata@meta.data
  mean_nfeature <- mean(tmdata$nFeature_RNA, na.rm = TRUE)
  median_nfeature <- median(tmdata$nFeature_RNA, na.rm = TRUE)

  mean_ncount <- mean(tmdata$nCount_RNA, na.rm = TRUE)
  median_ncount <- median(tmdata$nCount_RNA, na.rm = TRUE)

  mean_percent_mt <- mean(tmdata$percent.mt, na.rm = TRUE)
  median_percent_mt <- median(tmdata$percent.mt, na.rm = TRUE)

  base_theme <- theme(
    text = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    legend.position = "none"
  )

  build_single_violin <- function(column_name) {
    ggplot(plot_data, aes(x = orig.ident, y = .data[[column_name]])) +
      geom_violin(fill = "grey90", color = "grey40", scale = "width") +
      base_theme
  }

  feature_plot <- build_single_violin("nFeature_RNA") +
    geom_hline(yintercept = mean_nfeature, linetype = "dashed", color = "blue", linewidth = 0.5) +
    geom_hline(yintercept = median_nfeature, linetype = "solid", color = "red", linewidth = 0.5) +
    annotate(
      "text",
      x = 1.5,
      y = mean_nfeature,
      label = paste("Mean:", round(mean_nfeature, 1)),
      hjust = 0.5,
      vjust = -1,
      size = 3,
      color = "blue"
    ) +
    annotate(
      "text",
      x = 1.5,
      y = median_nfeature,
      label = paste("Median:", round(median_nfeature, 1)),
      hjust = 0.5,
      vjust = 1.5,
      size = 3,
      color = "red"
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste("NCells:", ncol(tmdata)),
      hjust = 1.1,
      vjust = 1.1,
      size = 3,
      color = "black"
    ) +
    ggtitle(sample_id)

  count_plot <- build_single_violin("nCount_RNA") +
    geom_hline(yintercept = mean_ncount, linetype = "dashed", color = "blue", linewidth = 0.5) +
    geom_hline(yintercept = median_ncount, linetype = "solid", color = "red", linewidth = 0.5) +
    annotate(
      "text",
      x = 1.5,
      y = mean_ncount,
      label = paste("Mean:", round(mean_ncount, 1)),
      hjust = 0.5,
      vjust = -1,
      size = 3,
      color = "blue"
    ) +
    annotate(
      "text",
      x = 1.5,
      y = median_ncount,
      label = paste("Median:", round(median_ncount, 1)),
      hjust = 0.5,
      vjust = 1.5,
      size = 3,
      color = "red"
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste("NCells:", ncol(tmdata)),
      hjust = 1.1,
      vjust = 1.1,
      size = 3,
      color = "black"
    ) +
    ggtitle(sample_id)

  mito_plot <- build_single_violin("percent.mt") +
    geom_hline(yintercept = mean_percent_mt, linetype = "dashed", color = "blue", linewidth = 0.5) +
    geom_hline(yintercept = median_percent_mt, linetype = "solid", color = "red", linewidth = 0.5) +
    annotate(
      "text",
      x = 1.5,
      y = mean_percent_mt,
      label = paste("Mean:", round(mean_percent_mt, 1)),
      hjust = 0.5,
      vjust = -1,
      size = 3,
      color = "blue"
    ) +
    annotate(
      "text",
      x = 1.5,
      y = median_percent_mt,
      label = paste("Median:", round(median_percent_mt, 1)),
      hjust = 0.5,
      vjust = 1.5,
      size = 3,
      color = "red"
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste("NCells:", ncol(tmdata)),
      hjust = 1.1,
      vjust = 1.1,
      size = 3,
      color = "black"
    ) +
    ggtitle(sample_id)

  list(
    feature = feature_plot,
    count = count_plot,
    mito = mito_plot
  )
}

build_cells_filtering_plot <- function(sample_id, n_genes, hk_mean, sl_cells_hk, ngenes, hkmean) {
  plot_data <- data.frame(
    hk_mean = hk_mean,
    n_genes = n_genes,
    valid = sl_cells_hk
  )

  ggplot(plot_data, aes(x = n_genes, y = hk_mean, color = valid)) +
    geom_point() +
    scale_x_continuous(trans = "log10", labels = scales::comma) +
    scale_y_continuous(trans = "log10", labels = scales::comma) +
    scale_color_manual(values = c("lightgrey", "black")) +
    labs(x = "Number of Genes", y = "HK Mean", color = "Valid Cells") +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(size = 8)) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste0("NCells passed: ", sum(sl_cells_hk)),
      hjust = 1.1,
      vjust = 1.1,
      size = 3,
      color = "black"
    ) +
    ggtitle(sample_id) +
    geom_vline(xintercept = ngenes, linetype = "dashed", color = "red") +
    geom_hline(yintercept = hkmean, linetype = "dashed", color = "red") +
    annotate(
      "text",
      x = ngenes,
      y = min(plot_data$hk_mean),
      label = paste0("Number of genes > ", ngenes),
      hjust = -0.1,
      vjust = 0,
      size = 3,
      color = "red"
    ) +
    annotate(
      "text",
      x = min(plot_data$n_genes),
      y = hkmean,
      label = paste0("HK Mean > ", hkmean),
      hjust = 0,
      vjust = -0.5,
      size = 3,
      color = "red"
    )
}

plot_chunks <- function(plot_list) {
  split(plot_list, ceiling(seq_along(plot_list) / 6))
}

get_grid_dims <- function(n) {
  if (n == 1) return(c(1, 1))
  if (n == 2) return(c(1, 2))
  if (n == 3) return(c(2, 2))
  if (n == 4) return(c(2, 2))
  c(2, 3)
}

write_chunked_pdf <- function(pdf_path, plot_lists, width = 8, height = 11) {
  grDevices::pdf(pdf_path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (plot_list in plot_lists) {
    if (length(plot_list) == 0) {
      next
    }
    for (chunk in plot_chunks(plot_list)) {
      dims <- get_grid_dims(length(chunk))
      gridExtra::grid.arrange(grobs = chunk, ncol = dims[1], nrow = dims[2])
    }
  }
}

normalise_object <- function(tmdata) {
  counts <- tmdata@assays$RNA$counts
  lib_size <- Matrix::colSums(counts)
  lib_size[lib_size == 0] <- 1
  cpm <- t(t(counts) / lib_size) * 1e6
  cpm <- as(cpm, "CsparseMatrix")
  expr <- log2((cpm / 10) + 1)
  expr <- as(expr, "CsparseMatrix")
  tmdata@assays$RNA$CPM <- cpm
  tmdata@assays$RNA$data <- expr
  tmdata
}

filter_cells <- function(tmdata, sample_id, ngenes, hkmean) {
  expr <- tmdata@assays$RNA$data
  n_genes <- Matrix::colSums(expr > 0)

  hk_list <- housekeeping[housekeeping %in% rownames(expr)]
  hk_mean <- if (length(hk_list) > 0) {
    Matrix::colMeans(expr[hk_list, , drop = FALSE])
  } else {
    rep(0, ncol(expr))
  }

  sl_cells_g <- n_genes >= ngenes
  sl_cells_hk <- hk_mean >= hkmean & sl_cells_g
  p <- build_cells_filtering_plot(
    sample_id = sample_id,
    n_genes = n_genes,
    hk_mean = hk_mean,
    sl_cells_hk = sl_cells_hk,
    ngenes = ngenes,
    hkmean = hkmean
  )

  filtered <- if (sum(sl_cells_hk) > 1) {
    subset(tmdata, cells = names(sl_cells_hk)[sl_cells_hk])
  } else {
    NULL
  }

  list(
    obj = filtered,
    g_filter = sum(sl_cells_g),
    hk_filter = sum(sl_cells_hk),
    plot = p
  )
}

manual_celltyping <- function(tmdata) {
  expr <- tmdata@assays$RNA$data
  tp_markers_list <- lapply(markers_list, function(marker_set) {
    marker_set[marker_set %in% rownames(expr)]
  })

  score_matrix <- matrix(
    0,
    nrow = ncol(expr),
    ncol = length(tp_markers_list),
    dimnames = list(colnames(expr), names(tp_markers_list))
  )

  for (celltype in names(tp_markers_list)) {
    marker_genes <- tp_markers_list[[celltype]]
    if (length(marker_genes) == 0) {
      next
    }
    marker_expr <- expr[marker_genes, , drop = FALSE]
    score_matrix[, celltype] <- Matrix::colMeans(marker_expr)
  }

  tmdata@meta.data <- cbind(tmdata@meta.data, score_matrix)
  tmdata$celltype_manual <- apply(score_matrix, 1, function(score_row) {
    valid_row <- score_row[is.finite(score_row)]
    if (length(valid_row) == 0) {
      return("unclassified")
    }
    max_score <- max(valid_row, na.rm = TRUE)
    max_name <- names(valid_row)[valid_row == max_score]
    if (length(max_name) != 1 || max_score < 1) {
      "unclassified"
    } else {
      max_name
    }
  })

  tmdata
}

x_filter <- setNames(numeric(length(names_tmdata)), basename(names_tmdata))
mt_filter <- setNames(numeric(length(names_tmdata)), basename(names_tmdata))
g_filter <- setNames(numeric(length(names_tmdata)), basename(names_tmdata))
hk_filter <- setNames(numeric(length(names_tmdata)), basename(names_tmdata))
inspection_feature_plots <- list()
inspection_count_plots <- list()
inspection_mito_plots <- list()
cells_filtering_plots <- list()

sample_manifest <- list()
dir.create("by_samples", recursive = TRUE, showWarnings = FALSE)
writeLines(names_tmdata, "names_tmdata_sn.txt")

for (path in names_tmdata) {
  sample_id <- basename(path)
  tmdata <- read_sample(path)
  message("finished reading ", path)

  x_filter[[sample_id]] <- ncol(tmdata)
  tmdata[["percent.mt"]] <- PercentageFeatureSet(tmdata, pattern = "^MT-")
  inspection_plots <- build_inspection_plots(tmdata, sample_id)
  inspection_feature_plots[[sample_id]] <- inspection_plots$feature
  inspection_count_plots[[sample_id]] <- inspection_plots$count
  inspection_mito_plots[[sample_id]] <- inspection_plots$mito

  if (sum(tmdata$percent.mt < max_mt) > 1) {
    tmdata <- subset(tmdata, percent.mt < max_mt)
    mt_filter[[sample_id]] <- ncol(tmdata)
  } else {
    mt_filter[[sample_id]] <- 0
    rm(tmdata)
    gc(verbose = FALSE)
    next
  }

  tmdata <- normalise_object(tmdata)
  message("finished normalising ", sample_id)

  filtered_out <- filter_cells(tmdata, sample_id, min_ngenes, min_hk_expr)
  cells_filtering_plots[[sample_id]] <- filtered_out$plot
  g_filter[[sample_id]] <- filtered_out$g_filter
  hk_filter[[sample_id]] <- filtered_out$hk_filter

  tmdata <- filtered_out$obj
  if (is.null(tmdata) || ncol(tmdata) <= 1) {
    rm(tmdata, filtered_out)
    gc(verbose = FALSE)
    next
  }

  tmdata <- manual_celltyping(tmdata)
  sample_dir <- file.path("by_samples", sample_id)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(tmdata, file.path(sample_dir, paste0(sample_id, ".rds")), compress = FALSE)

  sample_manifest[[sample_id]] <- data.frame(
    sample = sample_id,
    source_path = path,
    input_source = unique(tmdata$input_source)[1],
    technology = unique(tmdata$technology)[1],
    reference_batch = unique(tmdata$reference_batch)[1],
    cells = ncol(tmdata),
    stringsAsFactors = FALSE
  )

  rm(tmdata, filtered_out)
  gc(verbose = FALSE)
}

####################
# write scRef-style multi-panel QC PDFs
####################
write_chunked_pdf(
  pdf_path = paste0("Inspections", batch, ".pdf"),
  plot_lists = list(inspection_feature_plots, inspection_count_plots, inspection_mito_plots)
)
write_chunked_pdf(
  pdf_path = paste0("cells_filtering", batch, ".pdf"),
  plot_lists = list(cells_filtering_plots)
)

summary_tables <- lapply(
  list(x_filter, mt_filter, g_filter, hk_filter),
  function(x) tibble::enframe(x, name = "sample", value = "value")
)
sm_table <- Reduce(function(left, right) full_join(left, right, by = "sample"), summary_tables)
sm_table <- mutate(sm_table, across(-sample, ~ tidyr::replace_na(.x, 0)))

colnames(sm_table) <- c(
  "sample",
  "raw",
  paste0("mito_DNA\npercentage < ", max_mt),
  paste0("number of\ngenes > ", min_ngenes),
  paste0("housekeeping\nexpression > ", min_hk_expr)
)
write.csv(sm_table, paste0("filtering_summary", batch, ".csv"), row.names = FALSE)

sample_manifest_df <- if (length(sample_manifest) > 0) {
  do.call(rbind, sample_manifest)
} else {
  data.frame(
    sample = character(0),
    source_path = character(0),
    input_source = character(0),
    technology = character(0),
    reference_batch = character(0),
    cells = numeric(0),
    stringsAsFactors = FALSE
  )
}

write.csv(sample_manifest_df, "sample_manifest.csv", row.names = FALSE)
saveRDS(sample_manifest_df, paste0("snSeq_list", batch, ".rds"))
