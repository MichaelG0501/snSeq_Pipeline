####################
# Analysis registry:
#   Status: active shared library; not a standalone workflow script
#   Script: analysis/lib/config.R
#   Methodology: analysis/methodology/lib/shared_analysis_library_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Input: none
#   Output: shared constants/functions in the calling R session
####################

####################
# Shared analysis configuration for snSeq_Pipeline.
#
# This file centralizes constants used by scripts under analysis/.
# Scripts should source this file before defining workflow-specific paths.
####################

PROJECT_DIR <- normalizePath(
  Sys.getenv(
    "SNSEQ_PIPELINE_DIR",
    "/rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline"
  ),
  mustWork = FALSE
)
SN_OUTS_DIR <- file.path(PROJECT_DIR, "sn_outs")
ANALYSIS_DIR <- file.path(PROJECT_DIR, "analysis")
ANALYSIS_MAP_PATH <- file.path(ANALYSIS_DIR, "ANALYSIS_MAP.md")
SCREF_OUTS_DIR <- normalizePath(
  Sys.getenv(
    "SCREF_OUTS_DIR",
    "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs"
  ),
  mustWork = FALSE
)

PREFERRED_STATE_DEFINITION <- list(
  approach = "B",
  regression_mode = "noreg",
  label = "Centred refined scRef projection, noreg",
  final_state_file = file.path(SN_OUTS_DIR, "cell_states/centred/intermediate/centred_refined_noreg_states.rds"),
  primary_state_file = file.path(SN_OUTS_DIR, "cell_states/centred/intermediate/centred_refined_noreg_states.rds")
)

OUTPUT_TIERS <- c("intermediate", "tables", "figures", "logs", "reports")

METADATA_COLUMNS <- list(
  sample = "orig.ident",
  reference_batch = "reference_batch",
  technology = "technology",
  infercna_classification = "classification",
  malignancy = "malignancy",
  initial_annotation = "celltype_initial",
  updated_annotation = "celltype_update"
)

QC_THRESHOLDS <- list(
  max_percent_mt = 15,
  min_genes = 200,
  min_housekeeping_expression = 0.5
)

STATE_THRESHOLDS <- list(
  approach_b_min_best_score = 0.5,
  approach_b_hybrid_gap = 0.3
)

PLOT_DEFAULTS <- list(
  dpi = 300,
  png_dpi = 400,
  pptx_base_size = 16,
  pptx_axis_size = 14,
  pptx_legend_size = 14,
  pptx_title_size = 18,
  pptx_point_size = 1.2,
  heatmap_row_font = 11,
  heatmap_column_font = 9,
  heatmap_legend_font = 11
)

####################
# Current fixed centred scRef display definitions; source grouping is validated
# and copied by the scoring producer before current assignment.
MP_DESCRIPTIONS <- c(
  "MP1" = "G2/M cell cycle",
  "MP5" = "G1/S cell cycle",
  "MP13+" = "Single-nucleus cell cycle",
  "MP2+" = "MYC driven biosynthesis",
  "MP14" = "Squamoid/basal transition",
  "MP3+" = "Basal-columnar invasive epithelium",
  "MP6+" = "Inflammatory-reactive columnar epithelium",
  "MP11+" = "Epithelial type I interferon response",
  "MP9+" = "Metabolic columnar epithelium",
  "MP10+" = "Intestinal metaplasia",
  "MP18b" = "Mucous-secretory differentiation",
  "MP16" = "Mucous-secretory glandular epithelium",
  "MP17" = "MHC-II glandular progenitor",
  "MP8b" = "Metabolic intestinal metaplasia",
  "MP8+" = "Glandular intestinal metaplasia",
  "MP12" = "Hypoxic-inflammatory adaptive plasticity",
  "MP15" = "Cancer-cell immune mimicry"
)

CC_MPS <- c("MP1", "MP5", "MP13+")
STATE_GROUPS <- list(
  "Classic proliferation" = c("MP2+"),
  "Squamous-to-intestinal" = c("MP14", "MP3+", "MP6+", "MP11+", "MP9+", "MP10+"),
  "Glandular-to-intestinal" = c("MP18b", "MP16", "MP17", "MP8b", "MP8+"),
  "Stress-adaptive" = c("MP12"),
  "Cancer-cell immune mimicry" = c("MP15")
)
STATE_COLORS <- c(
  "Classic proliferation" = "#E41A1C", "Squamous-to-intestinal" = "#4DAF4A",
  "Glandular-to-intestinal" = "#FF7F00", "Stress-adaptive" = "#984EA3",
  "Cancer-cell immune mimicry" = "#377EB8", "Unresolved" = "grey80", "Hybrid" = "black"
)
####################

# Historical retained 3CA mappings are diagnostic/legacy provenance only.
RETAINED_3CA_STATE_MAP <- c(
  "X3CA_mp_30.Respiration.1" = "Classic Proliferative",
  "X3CA_mp_12.Protein.maturation" = "3CA_EMT_and_Protein_maturation",
  "X3CA_mp_17.EMT.III" = "3CA_EMT_and_Protein_maturation"
)

THREE_CA_CELL_CYCLE_MPS <- c(
  "X3CA_mp_1.Cell.Cycle...G2.M",
  "X3CA_mp_2.Cell.Cycle...G1.S",
  "X3CA_mp_3.Cell.Cylce.HMG.rich",
  "X3CA_mp_4.Chromatin",
  "X3CA_mp_5.Cell.cycle.single.nucleus"
)

analysis_output_dirs <- function(stage, root = SN_OUTS_DIR) {
  base_dir <- file.path(root, stage)
  dirs <- setNames(file.path(base_dir, OUTPUT_TIERS), OUTPUT_TIERS)
  c(base = base_dir, dirs)
}

ensure_output_dirs <- function(stage, root = SN_OUTS_DIR) {
  dirs <- analysis_output_dirs(stage, root = root)
  for (path in dirs) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  dirs
}

parse_key_value_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  if (length(args) == 0) {
    return(out)
  }
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) {
      next
    }
    parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
    key <- parts[1]
    value <- paste(parts[-1], collapse = "=")
    out[[key]] <- value
  }
  out
}

arg_value <- function(args, key, default = NULL) {
  value <- args[[key]]
  if (is.null(value) || length(value) == 0 || is.na(value) || !nzchar(value)) {
    return(default)
  }
  value
}

arg_flag <- function(args, key, default = FALSE) {
  value <- arg_value(args, key, default = NA_character_)
  if (is.na(value)) {
    return(default)
  }
  tolower(value) %in% c("true", "1", "yes", "y")
}

sn_out_path <- function(...) {
  file.path(SN_OUTS_DIR, ...)
}

analysis_path <- function(...) {
  file.path(ANALYSIS_DIR, ...)
}

####################
# Current centred reference paths. Historical 3CA maps above are diagnostic
# provenance only and do not modify current centred state labels.
CENTRED_REFERENCE_DIR <- file.path(SCREF_OUTS_DIR, "Metaprogrammes_Results", "centred", "mp_refinement")
CENTRED_OUTPUT_DIR <- file.path(SN_OUTS_DIR, "cell_states", "centred")
CENTRED_MP_GENES_SOURCE <- file.path(CENTRED_REFERENCE_DIR, "intermediate", "merged_refined_mp_genes.rds")
CENTRED_MP_GROUPING_SOURCE <- file.path(CENTRED_REFERENCE_DIR, "tables", "centred_refined_mp_state_grouping.csv")
CENTRED_MP_GENES <- file.path(CENTRED_OUTPUT_DIR, "intermediate", "merged_refined_mp_genes.rds")
CENTRED_UCELL <- file.path(CENTRED_OUTPUT_DIR, "intermediate", "merged_refined_ucell_scores.rds")
CENTRED_STATES <- file.path(CENTRED_OUTPUT_DIR, "intermediate", "centred_refined_noreg_states.rds")
CENTRED_MP_ADJ <- file.path(CENTRED_OUTPUT_DIR, "intermediate", "centred_refined_noreg_mp_adj.rds")
CENTRED_GROUP_MAX <- file.path(CENTRED_OUTPUT_DIR, "intermediate", "centred_refined_noreg_group_max.rds")
CENTRED_STATE_COLORS <- STATE_COLORS
####################
