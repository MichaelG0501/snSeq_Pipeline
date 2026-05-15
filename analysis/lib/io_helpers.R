####################
# Analysis registry:
#   Status: active shared library; not a standalone workflow script
#   Script: analysis/lib/io_helpers.R
#   Methodology: analysis/methodology/lib/shared_analysis_library_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Input: file paths supplied by calling scripts
#   Output: shared I/O helper functions in the calling R session
####################

####################
# Shared lightweight file and status helpers.
####################

read_status <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  trimws(readLines(path, warn = FALSE, n = 1))
}

read_optional_csv <- function(path, ...) {
  if (!file.exists(path)) {
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE, ...)
}

require_files <- function(paths, label = "required input") {
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths) > 0) {
    stop(
      "Missing ", label, if (length(missing_paths) > 1) "s" else "", ": ",
      paste(missing_paths, collapse = ", ")
    )
  }
  invisible(paths)
}

write_csv_checked <- function(x, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, quote = TRUE, ...)
  invisible(path)
}

save_rds_checked <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path)
  invisible(path)
}

extract_layer_matrix_subset <- function(seurat_obj, genes, assay = "RNA", layer_prefix = "data") {
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  if (length(genes) == 0) {
    stop("No genes provided for layer extraction.")
  }

  assay_obj <- seurat_obj[[assay]]
  layer_names <- Layers(assay_obj)
  target_layers <- layer_names[grepl(paste0("^", layer_prefix, "(\\.|$)"), layer_names)]

  if (length(target_layers) == 0) {
    stop("No layers found with prefix '", layer_prefix, "' in assay ", assay)
  }

  layer_mats <- lapply(target_layers, function(layer_name) {
    layer_mat <- LayerData(seurat_obj, assay = assay, layer = layer_name)
    common_genes <- intersect(genes, rownames(layer_mat))
    out_mat <- Matrix::Matrix(
      0,
      nrow = length(genes),
      ncol = ncol(layer_mat),
      sparse = TRUE,
      dimnames = list(genes, colnames(layer_mat))
    )
    if (length(common_genes) > 0) {
      out_mat[match(common_genes, genes), ] <- layer_mat[common_genes, , drop = FALSE]
    }
    out_mat
  })

  merged_mat <- do.call(cbind, layer_mats)
  merged_mat[, Cells(seurat_obj), drop = FALSE]
}
