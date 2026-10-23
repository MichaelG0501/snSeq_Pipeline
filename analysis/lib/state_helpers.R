####################
# Analysis registry:
#   Status: active shared library; not a standalone workflow script
#   Script: analysis/lib/state_helpers.R
#   Methodology: analysis/methodology/lib/shared_analysis_library_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Input: matrices/vectors supplied by calling state scripts
#   Output: shared state helper functions in the calling R session
####################

####################
# Shared state-analysis helper functions.
####################

z_normalise <- function(mat, sample_var, study_var) {
  clust_df <- as.data.frame(mat)
  clust_df$.cell <- rownames(mat)
  clust_df$.sample <- sample_var[rownames(mat)]
  clust_df$.study <- study_var[rownames(mat)]

  study_sd <- clust_df %>%
    group_by(.study) %>%
    summarise(across(all_of(colnames(mat)), ~ sd(.x, na.rm = TRUE)), .groups = "drop") %>%
    tibble::column_to_rownames(".study") %>%
    as.matrix()
  study_sd[is.na(study_sd) | study_sd == 0] <- 1

  clust_centered <- clust_df %>%
    group_by(.sample) %>%
    mutate(across(all_of(colnames(mat)), ~ .x - mean(.x, na.rm = TRUE))) %>%
    ungroup()

  mp_adj <- as.matrix(clust_centered[, colnames(mat), drop = FALSE])
  rownames(mp_adj) <- clust_centered$.cell
  for (mp in colnames(mp_adj)) {
    mp_adj[, mp] <- mp_adj[, mp] / study_sd[clust_centered$.study, mp]
  }
  mp_adj[!is.finite(mp_adj)] <- 0
  mp_adj
}

clean_3ca_name <- function(x) {
  x <- gsub("^X3CA_", "3CA_", x)
  gsub("\\.", " ", x)
}

format_3ca_label <- function(x) {
  x <- clean_3ca_name(x)
  sub("^3CA_mp_([0-9]+) ", "3CA MP\\1 ", x)
}

label_mp <- function(mp_vec, mp_descriptions = MP_DESCRIPTIONS) {
  desc <- mp_descriptions[mp_vec]
  desc[is.na(desc)] <- mp_vec[is.na(desc)]
  out <- paste0(mp_vec, "_", desc)
  names(out) <- names(mp_vec)
  out
}

state_order <- function(states = NULL) {
  base_order <- c(names(STATE_GROUPS), "Unresolved", "Hybrid")
  if (is.null(states)) {
    return(base_order)
  }
  unique(c(base_order[base_order %in% states], setdiff(states, base_order)))
}

state_colors_for <- function(states) {
  cols <- STATE_COLORS
  missing_states <- setdiff(states, names(cols))
  if (length(missing_states) > 0) {
    cols <- c(cols, setNames(scales::hue_pal()(length(missing_states)), missing_states))
  }
  cols[states]
}

build_group_max <- function(mp_adj, state_groups = STATE_GROUPS) {
  group_max <- sapply(state_groups, function(mps) {
    mps_avail <- intersect(mps, colnames(mp_adj))
    if (length(mps_avail) == 0) {
      return(rep(NA_real_, nrow(mp_adj)))
    }
    if (length(mps_avail) == 1) {
      return(as.numeric(mp_adj[, mps_avail]))
    }
    apply(mp_adj[, mps_avail, drop = FALSE], 1, max)
  })
  group_max <- as.matrix(group_max)
  rownames(group_max) <- rownames(mp_adj)
  group_max
}

read_preferred_states <- function() {
  path <- PREFERRED_STATE_DEFINITION$final_state_file
  if (!file.exists(path)) stop("Missing current centred state file: ", path)
  readRDS(path)
}

####################
# Current centred Approach-B thresholds; preserve exact barcode names.
call_approach_b_states <- function(group_max, thresholds = STATE_THRESHOLDS) {
  stopifnot(is.matrix(group_max), ncol(group_max) >= 2L,
            !is.null(rownames(group_max)), !anyDuplicated(rownames(group_max)),
            !is.null(colnames(group_max)), all(is.finite(group_max)))
  maximum <- apply(group_max, 1, max)
  ordered <- t(apply(group_max, 1, sort, decreasing = TRUE))
  gap <- ordered[, 1] - ordered[, 2]
  states <- colnames(group_max)[max.col(group_max, ties.method = "first")]
  states[maximum < thresholds$approach_b_min_best_score] <- "Unresolved"
  states[maximum >= thresholds$approach_b_min_best_score & gap < thresholds$approach_b_hybrid_gap] <- "Hybrid"
  names(states) <- rownames(group_max)
  states
}
####################
