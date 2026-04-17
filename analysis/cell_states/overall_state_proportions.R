####################
# Auto_overall_state_proportions.R
# Overall proportion barplot for noreg Approach B cell states.
####################

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs")

message("Loading data ...")
tmdata_all <- readRDS("snSeq_malignant_epi.rds")
final_states_path <- "Auto_final_states.rds"
if (file.exists(final_states_path)) {
  state_B <- readRDS(final_states_path)
} else {
  state_B <- readRDS("Auto_topmp_v2_noreg_states_B.rds")
}

common_cells <- intersect(names(state_B), Cells(tmdata_all))
state_B <- state_B[common_cells]

group_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "Immune Infiltrating" = "#377EB8",
  "3CA_EMT_and_Protein_maturation" = "#666666",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)

base_order <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "Immune Infiltrating"
)
extra_states <- setdiff(unique(as.character(state_B)), c(base_order, "Unresolved", "Hybrid"))
extra_states <- c("3CA_EMT_and_Protein_maturation", setdiff(extra_states, "3CA_EMT_and_Protein_maturation"))
extra_states <- extra_states[extra_states %in% unique(as.character(state_B))]
if (length(setdiff(extra_states, names(group_cols))) > 0) {
  new_cols <- setNames(scales::hue_pal()(length(setdiff(extra_states, names(group_cols)))), setdiff(extra_states, names(group_cols)))
  group_cols <- c(group_cols, new_cols)
}
state_order <- unique(c(base_order, extra_states, "Unresolved", "Hybrid"))

prop_df <- data.frame(
  state = factor(as.character(state_B), levels = state_order),
  stringsAsFactors = FALSE
)

overall <- prop_df %>%
  count(state, .drop = FALSE) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  mutate(label = sprintf("%.1f%%", pct))

message("Generating overall barplot ...")
p <- ggplot(overall, aes(x = "Overall", y = pct, fill = state)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.6) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 6,
    fontface = "bold",
    color = ifelse(overall$state %in% c("Hybrid"), "white", "black")
  ) +
  scale_fill_manual(values = group_cols, drop = FALSE) +
  labs(
    title = "scATLAS State Proportions",
    subtitle = paste0("Total cells: ", sum(overall$n)),
    x = NULL,
    y = "% of cells",
    fill = "State"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  coord_cartesian(expand = FALSE)

pdf("Auto_overall_state_proportions.pdf", width = 7, height = 9)
print(p)
dev.off()

message("Finished. Plot saved to sn_outs/Auto_overall_state_proportions.pdf")
