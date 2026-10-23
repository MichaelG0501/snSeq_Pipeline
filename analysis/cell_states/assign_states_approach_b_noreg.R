####################
# Analysis registry:
#   Status: active upstream/terminal
#   Script: analysis/cell_states/assign_states_approach_b_noreg.R
#   Methodology: analysis/methodology/cell_states/assign_states_approach_b_noreg_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Project current centred scRef MP groups into malignant snRNA cells.
#   Inputs (live sn_outs): snSeq_malignant_epi_meta.rds;
#     cell_states/centred/intermediate/merged_refined_ucell_scores.rds;
#     cell_states/centred/tables/centred_refined_mp_state_grouping.csv.
#   Outputs under live sn_outs/cell_states/centred/:
#     intermediate/centred_refined_noreg_{states,mp_adj,group_max}.rds;
#     intermediate/state_assignment_cache.rds;
#     tables/centred_refined_noreg_{sample_state_abundance,cell_assignments}.csv;
#     figures/Auto_final/01_centred_state_summary.pdf; logs/state_assignment.log.
#   Downstream: PDO legacy_Auto_pdo_sn_matched_pair_comparison.R (current route).
#   Cache: default rebuild; plot_only=true requires matching input fingerprints.
#   Run: qsub cell_states_master.sh (4 CPUs, 64 GB, 08:00:00; dmtcp).
####################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
})
source("/rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/analysis/lib/config.R")
source(file.path(ANALYSIS_DIR, "lib", "io_helpers.R"))
source(file.path(ANALYSIS_DIR, "lib", "state_helpers.R"))
source(file.path(ANALYSIS_DIR, "lib", "logging.R"))
setwd(SN_OUTS_DIR)

####################
args_kv <- parse_key_value_args()
plot_only <- arg_flag(args_kv, "plot_only", FALSE)
output_dirs <- ensure_output_dirs("cell_states/centred")
meta_path <- file.path(SN_OUTS_DIR, "snSeq_malignant_epi_meta.rds")
group_path <- file.path(output_dirs["tables"], "centred_refined_mp_state_grouping.csv")
input_files <- c(meta_path, CENTRED_UCELL, group_path)
require_files(input_files)
fingerprints <- tools::md5sum(input_files)
cache_path <- file.path(output_dirs["intermediate"], "state_assignment_cache.rds")
run_summary <- start_run_summary(
  "analysis/cell_states/assign_states_approach_b_noreg.R", inputs = input_files,
  outputs = c(CENTRED_STATES, CENTRED_MP_ADJ, CENTRED_GROUP_MAX),
  parameters = c(STATE_THRESHOLDS, list(batch = "reference_batch", plot_only = plot_only)),
  cached_reused = plot_only
)
meta <- readRDS(meta_path)
grouping <- read.csv(group_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("mp", "state", "description") %in% names(grouping)),
          !anyDuplicated(grouping$mp), !anyNA(grouping),
          all(c("orig.ident", "reference_batch", "malignancy") %in% names(meta)),
          !anyDuplicated(rownames(meta)),
          all(meta$malignancy %in% c("malignant_level_1", "malignant_level_2")),
          !anyNA(meta[, c("orig.ident", "reference_batch")]))
biological_states <- unique(grouping$state[grouping$state != "Cell cycle"])
state_levels <- c(biological_states, "Unresolved", "Hybrid")
stopifnot(setequal(state_levels, names(CENTRED_STATE_COLORS)))
state_groups <- setNames(lapply(biological_states, function(x) grouping$mp[grouping$state == x]), biological_states)
non_cc_mps <- unlist(state_groups, use.names = FALSE)
if (plot_only) {
  require_files(cache_path)
  cached <- readRDS(cache_path)
  if (!identical(cached$fingerprints, fingerprints)) stop("Input files changed; rebuild state assignment before plot_only.")
  states <- cached$states
  mp_adj <- cached$mp_adj
  group_max <- cached$group_max
} else {
  scores <- as.matrix(readRDS(CENTRED_UCELL))
  stopifnot(!anyDuplicated(rownames(scores)), !anyDuplicated(colnames(scores)),
            setequal(rownames(scores), rownames(meta)),
            setequal(colnames(scores), grouping$mp), all(is.finite(scores)))
  scores <- scores[rownames(meta), grouping$mp, drop = FALSE]
  samples <- setNames(as.character(meta$orig.ident), rownames(meta))
  batches <- setNames(as.character(meta$reference_batch), rownames(meta))
  mp_adj <- z_normalise(scores, samples, batches)
  group_max <- build_group_max(mp_adj, state_groups)
  stopifnot(identical(rownames(group_max), rownames(meta)), all(is.finite(group_max)))
  states <- call_approach_b_states(group_max)
  saveRDS(states, CENTRED_STATES)
  saveRDS(mp_adj, CENTRED_MP_ADJ)
  saveRDS(group_max, CENTRED_GROUP_MAX)
  saveRDS(list(states = states, mp_adj = mp_adj, group_max = group_max,
               fingerprints = fingerprints, thresholds = STATE_THRESHOLDS), cache_path)
}
stopifnot(identical(names(states), rownames(meta)), identical(rownames(mp_adj), names(states)),
          all(states %in% state_levels))
top_mp <- colnames(mp_adj)[max.col(mp_adj, ties.method = "first")]
top_state_mp <- non_cc_mps[max.col(mp_adj[, non_cc_mps, drop = FALSE], ties.method = "first")]
assignments <- data.frame(
  cell = names(states), sample = as.character(meta$orig.ident),
  reference_batch = as.character(meta$reference_batch), state = unname(states),
  top_mp = top_mp, top_state_mp = top_state_mp,
  top_state_mp_group = grouping$state[match(top_state_mp, grouping$mp)]
)
write.csv(assignments, file.path(output_dirs["tables"], "centred_refined_noreg_cell_assignments.csv"), row.names = FALSE)
abundance <- assignments %>% count(sample, reference_batch, state, name = "n") %>%
  complete(nesting(sample, reference_batch), state = state_levels, fill = list(n = 0L)) %>%
  group_by(sample, reference_batch) %>% mutate(total_n = sum(n), pct = 100 * n / total_n) %>% ungroup()
stopifnot(sum(abundance$n) == nrow(meta))
write.csv(abundance, file.path(output_dirs["tables"], "centred_refined_noreg_sample_state_abundance.csv"), row.names = FALSE)
write.csv(data.frame(state = state_levels, n = as.integer(table(factor(states, levels = state_levels)))),
          file.path(output_dirs["tables"], "centred_refined_noreg_state_counts.csv"), row.names = FALSE)

final_dir <- file.path(output_dirs["figures"], "Auto_final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
pdf_path <- file.path(final_dir, "01_centred_state_summary.pdf")
sample_order <- unique(abundance$sample[order(abundance$reference_batch, abundance$sample)])
abundance$sample <- factor(abundance$sample, levels = sample_order)
abundance$state <- factor(abundance$state, levels = state_levels)
pdf(pdf_path, width = 18, height = 9, onefile = TRUE)
print(ggplot(abundance, aes(sample, pct, fill = state)) + geom_col() +
        facet_grid(~ reference_batch, scales = "free_x", space = "free_x") +
        scale_fill_manual(values = CENTRED_STATE_COLORS, drop = FALSE) +
        labs(x = NULL, y = "% malignant epithelial cells", fill = "State",
             title = "Current centred scRef states projected onto snRNA-seq") +
        theme_bw(base_size = 15) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10)))
set.seed(1)
plot_cells <- unlist(lapply(state_levels, function(s) {
  cells <- names(states)[states == s]
  if (length(cells) > 500L) sample(cells, 500L) else cells
}), use.names = FALSE)
plot_matrix <- t(mp_adj[plot_cells, grouping$mp, drop = FALSE])
rownames(plot_matrix) <- paste(grouping$mp, grouping$description)
lim <- max(as.numeric(quantile(abs(plot_matrix), 0.98)), 0.01)
ComplexHeatmap::draw(ComplexHeatmap::Heatmap(
  plot_matrix, name = "Adjusted score", col = circlize::colorRamp2(c(-lim, 0, lim), c("navy", "white", "firebrick3")),
  column_split = factor(states[plot_cells], levels = state_levels),
  cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
  row_names_gp = grid::gpar(fontsize = 13), column_title_gp = grid::gpar(fontsize = 13),
  column_title_rot = 30, use_raster = TRUE,
  column_title = "At most 500 cells per state; display sample only"
))
dev.off()
write.csv(data.frame(file = pdf_path, role = "Final cohort abundance and sampled MP heatmap"),
          file.path(output_dirs["reports"], "Auto_figure_index.csv"), row.names = FALSE)
run_summary <- finish_run_summary(run_summary, "ok")
write_run_summary(run_summary, file.path(output_dirs["logs"], "state_assignment.log"))
message("CENTRED_STATE_ASSIGNMENT_OK: ", length(states), " cells; ", length(biological_states), " biological states")
####################
