####################
# Auto_states_hybrid_pairwise_nodeplot.R
# Node plot for real states and pairwise hybrid proportions.
# Multi-class hybrids (>2 groups) are excluded from edge construction.
####################

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs")
task_prefix <- "task6"
out_dir <- paste0(task_prefix, "_hybrid_pairwise")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

real_states <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "Immune Infiltrating"
)

state_groups <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive" = c("MP13", "MP12"),
  "SMG-like Metaplasia" = c("MP18", "MP16"),
  "Immune Infiltrating" = c("MP15")
)

group_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "Immune Infiltrating" = "#377EB8",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)

state_B <- readRDS("Auto_topmp_v2_noreg_states_B.rds")
mp_adj <- readRDS("Auto_topmp_v2_noreg_mp_adj.rds")

group_max <- sapply(state_groups, function(mps) {
  mps_avail <- intersect(mps, colnames(mp_adj))
  if (length(mps_avail) == 0) return(rep(NA_real_, nrow(mp_adj)))
  if (length(mps_avail) == 1) return(as.numeric(mp_adj[, mps_avail]))
  apply(mp_adj[, mps_avail, drop = FALSE], 1, max)
})
group_max <- as.matrix(group_max)
rownames(group_max) <- rownames(mp_adj)

hybrid_cells <- names(state_B)[state_B == "Hybrid"]

assign_pair <- function(x, names_vec) {
  ord <- names(sort(x, decreasing = TRUE))[1:2]
  ord <- names_vec[names_vec %in% ord]
  paste(ord, collapse = "__")
}

pair_labels <- vapply(hybrid_cells, function(cl) assign_pair(group_max[cl, real_states], real_states), character(1))

pair_df <- data.frame(pair = pair_labels, stringsAsFactors = FALSE) %>%
  count(pair, name = "hybrid_cells") %>%
  separate(pair, into = c("from", "to"), sep = "__", remove = FALSE)

state_df <- data.frame(state = state_B, stringsAsFactors = FALSE) %>%
  filter(state %in% real_states) %>%
  count(state, name = "cells")

tot_cells <- length(state_B)
state_df <- state_df %>% mutate(pct = 100 * cells / tot_cells)
pair_df <- pair_df %>% mutate(pct = 100 * hybrid_cells / tot_cells)

n <- length(real_states)
theta <- seq(0, 2 * pi, length.out = n + 1)[1:n]
layout_df <- data.frame(
  state = real_states,
  x = cos(theta),
  y = sin(theta),
  stringsAsFactors = FALSE
)

node_df <- left_join(layout_df, state_df, by = c("state" = "state"))
node_df$cells[is.na(node_df$cells)] <- 0
node_df$pct[is.na(node_df$pct)] <- 0
node_df$label_x <- node_df$x * 1.25
node_df$label_y <- node_df$y * 1.25

edge_df <- pair_df %>%
  left_join(layout_df, by = c("from" = "state")) %>%
  rename(x = x, y = y) %>%
  left_join(layout_df, by = c("to" = "state"), suffix = c("", "_to")) %>%
  rename(xend = x_to, yend = y_to)

p <- ggplot() +
  geom_segment(
    data = edge_df,
    aes(x = x, y = y, xend = xend, yend = yend, linewidth = pct),
    color = "grey35",
    alpha = 0.8
  ) +
  geom_point(
    data = node_df,
    aes(x = x, y = y, size = pct, color = state)
  ) +
  geom_text(
    data = node_df,
    aes(x = label_x, y = label_y, label = paste0(state, "\n", sprintf("%.1f%%", pct))),
    size = 3.5,
    fontface = "bold"
  ) +
  geom_label(
    data = edge_df,
    aes(
      x = (x + xend) / 2,
      y = (y + yend) / 2,
      label = sprintf("%.1f%%", pct)
    ),
    size = 2.6,
    fill = "white",
    label.size = 0,
    fontface = "bold"
  ) +
  scale_color_manual(values = group_cols) +
  scale_size(range = c(8, 22), guide = "none") +
  scale_linewidth(range = c(0.6, 6), guide = "none") +
  coord_equal() +
  expand_limits(x = c(-1.4, 1.4), y = c(-1.4, 1.4)) +
  theme_void(base_size = 14) +
  labs(title = "Pairwise hybrid network - noreg") +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(t = 10, b = 10))
  )

pdf(file.path(out_dir, "Auto_task6_hybrid_pairwise_nodeplot_noreg.pdf"), width = 6, height = 6)
print(p)
dev.off()

mat <- matrix(0, nrow = length(real_states), ncol = length(real_states), dimnames = list(real_states, real_states))
if (nrow(pair_df) > 0) {
  for (k in seq_len(nrow(pair_df))) {
    a <- pair_df$from[k]
    b <- pair_df$to[k]
    v <- pair_df$pct[k]
    if (a %in% real_states && b %in% real_states) {
      mat[a, b] <- v
      mat[b, a] <- v
    }
  }
}

hm_df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
colnames(hm_df) <- c("StateA", "StateB", "Pct")
p_hm <- ggplot(hm_df, aes(StateB, StateA, fill = Pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", Pct)), size = 3) +
  scale_fill_gradient(low = "white", high = "firebrick3") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(title = "Pairwise hybrid heatmap - noreg", x = NULL, y = NULL, fill = "% of all cells")

pdf(file.path(out_dir, "Auto_task6_hybrid_pairwise_heatmap_noreg.pdf"), width = 8, height = 6)
print(p_hm)
dev.off()

write.csv(
  bind_rows(
    state_df %>% transmute(type = "state", label = state, cells = cells, pct = pct),
    pair_df %>% transmute(type = "pairwise_hybrid", label = pair, cells = hybrid_cells, pct = pct)
  ),
  file.path(out_dir, "Auto_task6_hybrid_pairwise_nodeplot_summary.csv"),
  row.names = FALSE
)
