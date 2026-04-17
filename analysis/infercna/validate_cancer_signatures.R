#!/usr/bin/env Rscript
################################################################################
# Validate Malignant Signatures
# Computes the cancer signature score from scratch on Step 5 outputs (_epi.rds)
# Checks the consistency against InferCNA calls before running Malignancy.R
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# Move to the snSeq outputs folder
setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs")

# Helper function to compute the signature score
score_signature <- function(obj, genes) {
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(obj)))
  }
  colMeans(as.matrix(obj@assays$RNA$data[genes, , drop = FALSE]), na.rm = TRUE)
}

sig_file <- "cancer_signatures.txt"
global_cs_genes <- character(0)
if (file.exists(sig_file) && file.info(sig_file)$size > 0) {
  global_cs_genes <- read.table(sig_file, header = FALSE, stringsAsFactors = FALSE)$V1
} else {
  stop("cancer_signatures.txt not found. Cannot validate signature expression.")
}

# Identify samples that successfully derived signatures
# These are the samples listed in cancer_signatures_summary.csv
sum_path <- "cancer_signatures_summary.csv"
if (!file.exists(sum_path)) {
  stop("cancer_signatures_summary.csv not found. Step 5 (InferCNA.R) must complete for target samples first.")
}
sig_summary <- read.csv(sum_path, stringsAsFactors = FALSE)
valid_samples <- unique(unlist(strsplit(sig_summary$samples, "\\|")))

# Load cell cycle genes
cc_ref_path <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv"
if (!file.exists(cc_ref_path)) {
  stop("Cell cycle reference file not found.")
}
cell_cycle_ref <- read.csv(cc_ref_path, header = TRUE, stringsAsFactors = FALSE)[, 1:3]
consensus_cc_genes <- cell_cycle_ref$Gene[cell_cycle_ref$Consensus == 1]

all_meta <- list()

message(sprintf("Evaluating cancer and cell cycle signatures on %d qualified samples...", length(valid_samples)))

for (sample in valid_samples) {
  epi_file <- file.path("by_samples", sample, paste0(sample, "_epi.rds"))
  
  if (!file.exists(epi_file)) {
    message("  Skipping ", sample, " - RDS not found.")
    next
  }
  
  epi <- readRDS(epi_file)
  
  # Ensure the object has the requisite columns from InferCNA.R (Step 5)
  if (all(c("classification", "malignant_clus") %in% colnames(epi@meta.data))) {
    
    # Signature 1: Cancer Signatures (CS)
    cs_genes <- intersect(global_cs_genes, rownames(epi))
    if (length(cs_genes) > 0) {
      gene_means_cs <- rowMeans(epi@assays$RNA$data[cs_genes, , drop = FALSE], na.rm = TRUE)
      cs_genes_top <- names(sort(gene_means_cs, decreasing = TRUE))[seq_len(min(50, length(gene_means_cs)))]
    } else {
      cs_genes_top <- character(0)
    }
    epi$cs_score <- score_signature(epi, cs_genes_top)

    # Signature 2: Cell Cycle (CC)
    cc_genes <- intersect(consensus_cc_genes, rownames(epi))
    if (length(cc_genes) > 0) {
      gene_means_cc <- rowMeans(epi@assays$RNA$data[cc_genes, , drop = FALSE], na.rm = TRUE)
      cc_genes_top <- names(sort(gene_means_cc, decreasing = TRUE))[seq_len(min(10, length(gene_means_cc)))]
    } else {
      cc_genes_top <- character(0)
    }
    epi$cc_score <- score_signature(epi, cc_genes_top)
    
    # Store metadata for plotting
    meta <- epi@meta.data[, c("orig.ident", "classification", "malignant_clus", "cs_score", "cc_score"), drop = FALSE]
    all_meta[[sample]] <- meta
    message("  Processed ", sample)
  }
}

if (length(all_meta) == 0) {
  stop("No valid _epi.rds objects found to analyze.")
}

df <- bind_rows(all_meta)

message(sprintf("Computed scores for %d single cells across samples.", nrow(df)))

# ----------------------------------------------------------------------------------
# Plot 1: Compare the Sample Cancer Signatures at the Cell-Level InferCNA Classification
# ----------------------------------------------------------------------------------
df_cna <- df %>% filter(classification %in% c("cna_malignant", "cna_non_malignant", "cna_unresolved"))

df_cna$classification <- factor(df_cna$classification, 
                                 levels = c("cna_non_malignant", "cna_unresolved", "cna_malignant"),
                                 labels = c("Non-Malignant\n(CNA-)", "Unresolved\n(CNA?)", "Malignant\n(CNA+)"))

class_colors <- c("Non-Malignant\n(CNA-)" = "#377EB8", "Unresolved\n(CNA?)" = "#BDBDBD", "Malignant\n(CNA+)" = "#E41A1C")

p1 <- ggplot(df_cna, aes(x = classification, y = cs_score, fill = classification)) +
  geom_violin(alpha = 0.5, scale = "width", trim = TRUE, color = NA) +
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values = class_colors) +
  theme_classic(base_size = 14) +
  labs(
    title = "Cancer Signatures vs Cell-Level InferCNA",
    subtitle = "Cell-level continuous CNA baseline",
    x = "InferCNA Signal Status",
    y = "Cancer Signature Score (cs_score)"
  ) +
  theme(
    legend.position = "none", 
    plot.title = element_text(face = "bold", hjust = 0.5), 
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    axis.text.x = element_text(face = "bold")
  )

# ----------------------------------------------------------------------------------
# Plot 2: Compare the Sample Cancer Signatures at the Cluster-Level Consensus
# ----------------------------------------------------------------------------------
df_clus <- df %>% filter(malignant_clus %in% c("malignant_clus", "non_malignant_clus", "unresolved_clus"))

df_clus$malignant_clus <- factor(df_clus$malignant_clus, 
                                 levels = c("non_malignant_clus", "unresolved_clus", "malignant_clus"),
                                 labels = c("Non-Malignant\nCluster", "Unresolved\nCluster", "Malignant\nCluster"))

clus_colors <- c("Non-Malignant\nCluster" = "#377EB8", "Unresolved\nCluster" = "#BDBDBD", "Malignant\nCluster" = "#E41A1C")

p3 <- ggplot(df_clus, aes(x = malignant_clus, y = cs_score, fill = malignant_clus)) +
  geom_violin(alpha = 0.5, scale = "width", trim = TRUE, color = NA) +
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values = clus_colors) +
  theme_classic(base_size = 14) +
  labs(
    title = "Cancer Signatures vs Cluster Consensus",
    subtitle = "Cluster-level (K-NN) CNA voting",
    x = "Seurat Cluster Identity",
    y = "Cancer Signature Score (cs_score)"
  ) +
  theme(
    legend.position = "none", 
    plot.title = element_text(face = "bold", hjust = 0.5), 
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    axis.text.x = element_text(face = "bold")
  )
  
# ----------------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------------

out_dir <- "analysis/infercna"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
pdf_path <- file.path(out_dir, "validation_cancer_signatures_vs_cna.pdf")

message("Done! Processing output...")

# ----------------------------------------------------------------------------------
# Page 2: Cell Cycle Scores
# ----------------------------------------------------------------------------------

p_cc1 <- ggplot(df_cna, aes(x = classification, y = cc_score, fill = classification)) +
  geom_violin(alpha = 0.5, scale = "width", trim = TRUE, color = NA) +
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values = class_colors) +
  theme_classic(base_size = 14) +
  labs(
    title = "Cell Cycle vs Cell-Level InferCNA",
    subtitle = "Cell-level continuous CNA baseline",
    x = "InferCNA Signal Status",
    y = "Cell Cycle Score (cc_score)"
  ) +
  theme(
    legend.position = "none", 
    plot.title = element_text(face = "bold", hjust = 0.5), 
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    axis.text.x = element_text(face = "bold")
  )

p_cc2 <- ggplot(df_clus, aes(x = malignant_clus, y = cc_score, fill = malignant_clus)) +
  geom_violin(alpha = 0.5, scale = "width", trim = TRUE, color = NA) +
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values = clus_colors) +
  theme_classic(base_size = 14) +
  labs(
    title = "Cell Cycle vs Cluster Consensus",
    subtitle = "Cluster-level (K-NN) CNA voting",
    x = "Seurat Cluster Identity",
    y = "Cell Cycle Score (cc_score)"
  ) +
  theme(
    legend.position = "none", 
    plot.title = element_text(face = "bold", hjust = 0.5), 
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    axis.text.x = element_text(face = "bold")
  )

pdf(pdf_path, width = 11, height = 12)
print((p1 + p3) / (p_cc1 + p_cc2))
dev.off()

message("Full report saved to: sn_outs/", pdf_path)
