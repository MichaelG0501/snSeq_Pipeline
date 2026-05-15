####################
# Analysis registry:
#   Status: active terminal figure-generation
#   Script: analysis/cell_states/plot_state_overall_proportions.R
#   Methodology: analysis/methodology/cell_states/plot_state_overall_proportions_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Output tiers: sn_outs/cell_states/overall_proportions/{intermediate,tables,figures,logs,reports}
####################

####################
# plot_state_overall_proportions.R
#
# Generate the overall malignant epithelial state composition barplot using
# the selected Approach-B noreg final states.
#
# Input:
#   sn_outs/snSeq_malignant_epi.rds
#   sn_outs/Auto_final_states.rds or sn_outs/Auto_topmp_v2_noreg_states_B.rds
#
# Output:
#   sn_outs/Auto_overall_state_proportions.pdf
#   sn_outs/cell_states/overall_proportions/tables/overall_state_proportions.csv
#   sn_outs/cell_states/overall_proportions/logs/plot_state_overall_proportions.log
#
# Usage:
#   Rscript analysis/cell_states/plot_state_overall_proportions.R
####################

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

source("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/analysis/lib/config.R")
source(file.path(ANALYSIS_DIR, "lib", "state_helpers.R"))
source(file.path(ANALYSIS_DIR, "lib", "logging.R"))

setwd(SN_OUTS_DIR)
output_dirs <- ensure_output_dirs("cell_states/overall_proportions")

run_summary <- start_run_summary(
  script = "analysis/cell_states/plot_state_overall_proportions.R",
  inputs = c(
    file.path(SN_OUTS_DIR, "snSeq_malignant_epi.rds"),
    PREFERRED_STATE_DEFINITION$final_state_file,
    PREFERRED_STATE_DEFINITION$primary_state_file
  ),
  outputs = c(
    file.path(SN_OUTS_DIR, "Auto_overall_state_proportions.pdf"),
    file.path(output_dirs["tables"], "overall_state_proportions.csv")
  ),
  parameters = list(state_method = PREFERRED_STATE_DEFINITION$label)
)

message("Loading data ...")
tmdata_all <- readRDS("snSeq_malignant_epi.rds")
state_B <- read_preferred_states()

common_cells <- intersect(names(state_B), Cells(tmdata_all))
state_B <- state_B[common_cells]

group_cols <- STATE_COLORS

base_order <- names(STATE_GROUPS)
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

write.csv(overall, file.path(output_dirs["tables"], "overall_state_proportions.csv"), row.names = FALSE)

run_summary <- finish_run_summary(run_summary, status = "ok")
write_run_summary(
  run_summary,
  file.path(output_dirs["logs"], "plot_state_overall_proportions.log")
)

message("Finished. Plot saved to sn_outs/Auto_overall_state_proportions.pdf")
