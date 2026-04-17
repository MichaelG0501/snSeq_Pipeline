suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("patchwork")
  library("Seurat")
  library("infercna")
  library("scales")
})

# WD
setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline")

out_dir <- "sn_outs"
output_pdf <- file.path(out_dir, "Auto_malignancy_scatter_comparison.png")

# Samples: Good, Medium, Bad
samples <- c("C_post_N2_biopsy", "D_pre_T0_biopsy", "N_post_T1_biopsy")
names(samples) <- c("Good (Mostly CNA Positive)", 
                    "Medium (Mixed Profile)", 
                    "Bad (Mostly CNA Negative)")

plot_list <- list()

# Define color palette matching professional infercna style
cna_colors <- c(
  "CNA Positive" = "#D73027",
  "CNA Unresolved" = "#FDAE61",
  "CNA Negative" = "#4575B4",
  "Reference" = "#D9D9D9"
)

for (i in seq_along(samples)) {
  sample_id <- samples[i]
  sample_label <- names(samples)[i]
  
  outs_path <- file.path(out_dir, "by_samples", sample_id, paste0(sample_id, "_outs.rds"))
  epi_path <- file.path(out_dir, "by_samples", sample_id, paste0(sample_id, "_epi.rds"))
  summary_path <- file.path(out_dir, "by_samples", sample_id, paste0(sample_id, "_infercna_summary.csv"))
  
  if (!file.exists(outs_path) || !file.exists(summary_path)) {
    message("Skipping ", sample_id, " - files not found")
    next
  }
  
  message("Processing ", sample_id, " (", sample_label, ")")
  outs <- readRDS(outs_path)
  summary_df <- read.csv(summary_path, stringsAsFactors = FALSE)
  
  coord <- cnaScatterPlot(outs)
  plot_df <- coord %>% as.data.frame()
  plot_df$cell_id <- rownames(plot_df)
  
  thr_cor <- summary_df$thr_cor[1]
  thr_sig <- summary_df$thr_sig[1]
  
  if (file.exists(epi_path)) {
    epi <- readRDS(epi_path)
    epi_meta <- epi@meta.data %>% as.data.frame()
    epi_meta$cell_id <- rownames(epi_meta)
    plot_df <- plot_df %>% left_join(epi_meta[, c("cell_id", "classification")], by = "cell_id")
  } else {
    plot_df$classification <- NA_character_
  }
  
  plot_df <- plot_df %>%
    mutate(
      Status = case_when(
        classification == "cna_malignant" ~ "CNA Positive",
        classification == "cna_unresolved" ~ "CNA Unresolved",
        classification == "cna_non_malignant" ~ "CNA Negative",
        TRUE ~ "Reference"
      )
    ) %>%
    mutate(
      Status = factor(Status, levels = names(cna_colors))
    )
  
  # Plot
  p <- ggplot(plot_df, aes(x = cna.cor, y = cna.signal)) +
    geom_point(data = subset(plot_df, Status == "Reference"), 
               color = "#E0E0E0", size = 0.4, alpha = 0.3) +
    geom_point(data = subset(plot_df, Status != "Reference"), 
               aes(color = Status), size = 0.75, alpha = 0.75) +
    geom_vline(xintercept = thr_cor, linetype = "dashed", color = "gray30", linewidth = 0.4) +
    geom_hline(yintercept = thr_sig, linetype = "dashed", color = "gray30", linewidth = 0.4) +
    scale_color_manual(values = cna_colors, drop = FALSE) +
    labs(
      title = sample_label,
      x = "CNA Correlation",
      y = "CNA Signal"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      legend.position = "none", # Individual legends hidden
      axis.line = element_line(linewidth = 0.6),
      axis.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "#F8F8F8", linewidth = 0.4)
    )
  
  plot_list[[i]] <- p
}

if (length(plot_list) == 3) {
  # Layout: 1x3 with collected legend
  combined_plot <- (plot_list[[1]] | plot_list[[2]] | plot_list[[3]]) + 
    plot_layout(guides = "collect") +
    plot_annotation(
      theme = theme(
        legend.position = "bottom",
        legend.text = element_text(size = 14, face = "bold"),
        legend.margin = margin(t = 10)
      )
    )
  
  ggsave(output_pdf, combined_plot, width = 18, height = 7, dpi = 300)
  message("Success: Created 3-sample CNA comparison plot at ", output_pdf)
}
