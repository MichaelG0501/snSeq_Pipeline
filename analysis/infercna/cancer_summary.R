suppressPackageStartupMessages({
  library("dplyr")
  library("tidyr")
  library("ggplot2")
  library("scales")
  library("patchwork")
  library("Seurat")
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]))
} else {
  normalizePath(getwd())
}

analysis_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(analysis_dir, "..", ".."))
out_dir <- file.path(project_dir, "sn_outs")

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

manifest_path <- file.path(out_dir, "sample_manifest.csv")
filtered_summary_path <- file.path(out_dir, "filtered_sample_summary.csv")
infercna_status_path <- file.path(out_dir, "infercna_status_by_sample.csv")

required_files <- c(manifest_path, filtered_summary_path, infercna_status_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required inputs: ", paste(basename(missing_files), collapse = ", "))
}

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
filtered_summary <- read.csv(filtered_summary_path, stringsAsFactors = FALSE)
infercna_status <- read.csv(infercna_status_path, stringsAsFactors = FALSE)

read_status <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  trimws(readLines(path, warn = FALSE, n = 1))
}

read_optional_csv <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE)
}

tech_colors <- c(
  "snSeq" = "steelblue",
  "GEMX" = "darkorange",
  "Multiome" = "forestgreen"
)

malignancy_colors <- c(
  malignant_level_1 = "#CB181D",
  malignant_level_2 = "#FB6A4A",
  malignant_unresolved = "#FCAE91",
  non_malignant_level_1 = "#2171B5",
  non_malignant_level_2 = "#6BAED6",
  non_malignant_unresolved = "#BDD7E7",
  unresolved = "grey55"
)
malignancy_levels <- names(malignancy_colors)

status_rows <- lapply(manifest$sample, function(sample) {
  expr_status <- read_status(file.path(out_dir, "by_samples", sample, "expr_filter_status.txt"))
  infer_status <- read_status(file.path(out_dir, "by_samples", sample, "infercna_status.txt"))
  malign_status <- read_status(file.path(out_dir, "by_samples", sample, "malignancy_status.txt"))
  epi_f_path <- file.path(out_dir, "by_samples", sample, paste0(sample, "_epi_f.rds"))
  sample_summary_path <- file.path(out_dir, "by_samples", sample, paste0(sample, "_malignancy_summary.csv"))
  cluster_summary_path <- file.path(out_dir, "by_samples", sample, paste0(sample, "_malignancy_cluster_summary.csv"))

  effective_malignancy_status <- malign_status
  if (identical(expr_status, "ok") && identical(infer_status, "ok") && identical(malign_status, "ok") &&
      !(file.exists(epi_f_path) && file.exists(sample_summary_path) && file.exists(cluster_summary_path))) {
    effective_malignancy_status <- "ok_missing_outputs"
  }

  data.frame(
    sample = sample,
    expr_filter_status = expr_status,
    infercna_status = infer_status,
    malignancy_status = malign_status,
    effective_malignancy_status = effective_malignancy_status,
    epi_f_present = file.exists(epi_f_path),
    malignancy_summary_present = file.exists(sample_summary_path),
    malignancy_cluster_summary_present = file.exists(cluster_summary_path),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows() %>%
  left_join(
    filtered_summary %>%
      select(sample, cells_before, cells_after, pct_retained, technology, reference_batch, filtered_status = status),
    by = "sample"
  )

write.csv(
  status_rows,
  file.path(out_dir, "malignancy_status_by_sample.csv"),
  row.names = FALSE,
  quote = TRUE
)

eligible_samples <- status_rows$sample[
  status_rows$expr_filter_status == "ok" &
    status_rows$infercna_status == "ok" &
    status_rows$effective_malignancy_status == "ok"
]

sample_summaries <- lapply(eligible_samples, function(sample) {
  read_optional_csv(file.path(out_dir, "by_samples", sample, paste0(sample, "_malignancy_summary.csv")))
}) %>% bind_rows()

cluster_summaries <- lapply(eligible_samples, function(sample) {
  out <- read_optional_csv(file.path(out_dir, "by_samples", sample, paste0(sample, "_malignancy_cluster_summary.csv")))
  if (is.null(out)) {
    return(NULL)
  }
  out$seurat_clusters <- as.character(out$seurat_clusters)
  out$malignant_clus <- as.character(out$malignant_clus)
  out
}) %>% bind_rows()

meta_rows <- list()
for (sample in eligible_samples) {
  epi_f_path <- file.path(out_dir, "by_samples", sample, paste0(sample, "_epi_f.rds"))
  epi <- readRDS(epi_f_path)
  keep_cols <- intersect(
    c(
      "orig.ident", "reference_batch", "technology", "classification",
      "malignant_clus", "cc_status", "cs_status", "cc_score", "cs_score",
      "malignancy", "seurat_clusters"
    ),
    colnames(epi@meta.data)
  )
  meta <- epi@meta.data[, keep_cols, drop = FALSE]
  meta$sample <- sample
  meta_rows[[sample]] <- meta
}
meta_all <- bind_rows(meta_rows)

if (nrow(meta_all) == 0) {
  stop("No completed malignancy outputs found.")
}

write.csv(
  sample_summaries,
  file.path(out_dir, "malignancy_sample_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  cluster_summaries,
  file.path(out_dir, "malignancy_cluster_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)

meta_all$malignancy <- factor(as.character(meta_all$malignancy), levels = malignancy_levels)
meta_all$malignancy_simple <- ifelse(
  meta_all$malignancy %in% c("malignant_level_1", "malignant_level_2"),
  "malignant",
  ifelse(
    meta_all$malignancy %in% c("non_malignant_level_1", "non_malignant_level_2"),
    "non_malignant",
    "other"
  )
)

overall_counts <- meta_all %>%
  count(malignancy, name = "cells") %>%
  complete(malignancy = factor(malignancy_levels, levels = malignancy_levels), fill = list(cells = 0))

overall_total <- sum(overall_counts$cells, na.rm = TRUE)
overall_counts$pct_of_epithelial <- if (overall_total > 0) {
  overall_counts$cells / overall_total
} else {
  NA_real_
}

write.csv(
  overall_counts,
  file.path(out_dir, "malignancy_classification_overall.csv"),
  row.names = FALSE,
  quote = TRUE
)

overall_summary <- data.frame(
  expr_filter_ok_samples = sum(status_rows$expr_filter_status == "ok", na.rm = TRUE),
  infercna_ok_samples = sum(status_rows$expr_filter_status == "ok" & status_rows$infercna_status == "ok", na.rm = TRUE),
  malignancy_ok_samples = length(eligible_samples),
  no_epi_after_infercna = sum(status_rows$expr_filter_status == "ok" & status_rows$infercna_status == "no_epi", na.rm = TRUE),
  epithelial_cells = nrow(meta_all),
  malignant_level_1_cells = sum(meta_all$malignancy == "malignant_level_1", na.rm = TRUE),
  malignant_level_2_cells = sum(meta_all$malignancy == "malignant_level_2", na.rm = TRUE),
  malignant_unresolved_cells = sum(meta_all$malignancy == "malignant_unresolved", na.rm = TRUE),
  malignant_total_cells = sum(meta_all$malignancy %in% c("malignant_level_1", "malignant_level_2"), na.rm = TRUE),
  non_malignant_level_1_cells = sum(meta_all$malignancy == "non_malignant_level_1", na.rm = TRUE),
  non_malignant_level_2_cells = sum(meta_all$malignancy == "non_malignant_level_2", na.rm = TRUE),
  non_malignant_unresolved_cells = sum(meta_all$malignancy == "non_malignant_unresolved", na.rm = TRUE),
  non_malignant_total_cells = sum(meta_all$malignancy %in% c("non_malignant_level_1", "non_malignant_level_2"), na.rm = TRUE),
  unresolved_cells = sum(meta_all$malignancy == "unresolved", na.rm = TRUE),
  cs_malignant_cells = sum(meta_all$cs_status == "cs_malignant", na.rm = TRUE),
  cc_malignant_cells = sum(meta_all$cc_status == "cc_malignant", na.rm = TRUE),
  stringsAsFactors = FALSE
)
overall_summary$pct_malignant_total <- overall_summary$malignant_total_cells / overall_summary$epithelial_cells
overall_summary$pct_non_malignant_total <- overall_summary$non_malignant_total_cells / overall_summary$epithelial_cells
overall_summary$pct_unresolved <- overall_summary$unresolved_cells / overall_summary$epithelial_cells

write.csv(
  overall_summary,
  file.path(out_dir, "malignancy_overall_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)

sample_level <- meta_all %>%
  count(sample, technology, reference_batch, malignancy, name = "cells") %>%
  complete(sample, malignancy = malignancy_levels, fill = list(cells = 0)) %>%
  left_join(
    filtered_summary %>% select(sample, cells_after, technology, reference_batch),
    by = "sample",
    suffix = c("", ".filtered")
  ) %>%
  mutate(
    technology = coalesce(technology, technology.filtered),
    reference_batch = coalesce(reference_batch, reference_batch.filtered)
  ) %>%
  select(-technology.filtered, -reference_batch.filtered) %>%
  group_by(sample) %>%
  mutate(
    epithelial_cells = sum(cells),
    pct_of_sample_epithelial = ifelse(epithelial_cells > 0, cells / epithelial_cells, NA_real_)
  ) %>%
  ungroup()

write.csv(
  sample_level,
  file.path(out_dir, "malignancy_by_sample.csv"),
  row.names = FALSE,
  quote = TRUE
)

sample_simple <- meta_all %>%
  group_by(sample) %>%
  summarise(
    epithelial_cells = n(),
    malignant_cells = sum(malignancy %in% c("malignant_level_1", "malignant_level_2"), na.rm = TRUE),
    non_malignant_cells = sum(malignancy %in% c("non_malignant_level_1", "non_malignant_level_2"), na.rm = TRUE),
    unresolved_cells = sum(malignancy == "unresolved", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_malignant = malignant_cells / epithelial_cells,
    pct_non_malignant = non_malignant_cells / epithelial_cells,
    pct_unresolved = unresolved_cells / epithelial_cells
  ) %>%
  left_join(
    filtered_summary %>% select(sample, cells_after, technology, reference_batch),
    by = "sample"
  ) %>%
  arrange(desc(pct_malignant), desc(epithelial_cells))

write.csv(
  sample_simple,
  file.path(out_dir, "malignancy_simple_by_sample.csv"),
  row.names = FALSE,
  quote = TRUE
)

overall_plot <- ggplot(overall_counts, aes(x = malignancy, y = cells, fill = malignancy)) +
  geom_col(width = 0.8, color = "gray20", linewidth = 0.25) +
  geom_text(aes(label = comma(cells)), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = malignancy_colors, drop = FALSE, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = comma) +
  labs(title = "Overall Malignancy Classification", x = NULL, y = "Epithelial cells") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

sample_simple$sample <- factor(sample_simple$sample, levels = sample_simple$sample)
malignant_fraction_plot <- ggplot(sample_simple, aes(x = sample, y = pct_malignant, fill = technology)) +
  geom_col(width = 0.85, color = "gray20", linewidth = 0.2) +
  scale_fill_manual(values = tech_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.04))) +
  labs(title = "Malignant Fraction Per Sample", x = NULL, y = "% malignant epithelial cells", fill = "Technology") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

sample_level$sample <- factor(sample_level$sample, levels = sample_simple$sample)
composition_plot <- ggplot(sample_level, aes(x = sample, y = pct_of_sample_epithelial, fill = malignancy)) +
  geom_col(width = 0.85, color = "gray20", linewidth = 0.1) +
  scale_fill_manual(values = malignancy_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Malignancy Composition Per Sample", x = NULL, y = "Fraction of epithelial cells", fill = "Malignancy") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

batch_level <- sample_level %>%
  group_by(reference_batch, malignancy) %>%
  summarise(cells = sum(cells), .groups = "drop") %>%
  group_by(reference_batch) %>%
  mutate(pct = cells / sum(cells)) %>%
  ungroup()

batch_plot <- ggplot(batch_level, aes(x = reference_batch, y = pct, fill = malignancy)) +
  geom_col(width = 0.8, color = "gray20", linewidth = 0.1) +
  scale_fill_manual(values = malignancy_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Malignancy Composition By Reference Batch", x = NULL, y = "Fraction of epithelial cells", fill = "Malignancy") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "right"
  )

ggsave(file.path(analysis_dir, "malignancy_overall_counts.png"), overall_plot, width = 9, height = 5.5, dpi = 300)
ggsave(file.path(analysis_dir, "malignant_fraction_by_sample.png"), malignant_fraction_plot, width = 12, height = 4.8, dpi = 300)
ggsave(file.path(analysis_dir, "malignancy_composition_by_sample.png"), composition_plot, width = 13, height = 5.4, dpi = 300)
ggsave(file.path(analysis_dir, "malignancy_composition_by_reference_batch.png"), batch_plot, width = 9, height = 5.5, dpi = 300)

overview <- (overall_plot + batch_plot) / (malignant_fraction_plot + composition_plot) +
  plot_layout(heights = c(1, 1.15))
ggsave(file.path(analysis_dir, "cancer_summary_overview.pdf"), overview, width = 14, height = 10.5)
