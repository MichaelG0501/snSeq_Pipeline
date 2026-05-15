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
    "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline"
  ),
  mustWork = FALSE
)
SN_OUTS_DIR <- file.path(PROJECT_DIR, "sn_outs")
ANALYSIS_DIR <- file.path(PROJECT_DIR, "analysis")
ANALYSIS_MAP_PATH <- file.path(ANALYSIS_DIR, "ANALYSIS_MAP.md")
SCREF_OUTS_DIR <- normalizePath(
  Sys.getenv(
    "SCREF_OUTS_DIR",
    "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs"
  ),
  mustWork = FALSE
)

PREFERRED_STATE_DEFINITION <- list(
  approach = "B",
  regression_mode = "noreg",
  label = "Approach B noreg",
  final_state_file = file.path(SN_OUTS_DIR, "Auto_final_states.rds"),
  primary_state_file = file.path(SN_OUTS_DIR, "Auto_topmp_v2_noreg_states_B.rds")
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

MP_DESCRIPTIONS <- c(
  "MP1"  = "G2M Cell Cycle",
  "MP9"  = "G1S Cell Cycle",
  "MP2"  = "MYC-related Proliferation",
  "MP17" = "Basal-like Transition",
  "MP14" = "Hypoxia Adapted Epi.",
  "MP5"  = "Epithelial IFN Resp.",
  "MP10" = "Columnar Diff.",
  "MP8"  = "Intestinal Diff.",
  "MP13" = "Hypoxic Inflam. Epi.",
  "MP7"  = "DNA Damage Repair",
  "MP18" = "Secretory Diff. (Intest.)",
  "MP16" = "Secretory Diff. (Gastric)",
  "MP15" = "Immune Infiltration",
  "MP12" = "Neuro-responsive Epi"
)

CC_MPS <- c("MP1", "MP7", "MP9")

STATE_GROUPS <- list(
  "Classic Proliferative" = c("MP2"),
  "Basal to Intestinal Metaplasia" = c("MP17", "MP14", "MP5", "MP10", "MP8"),
  "Stress-adaptive" = c("MP13", "MP12"),
  "SMG-like Metaplasia" = c("MP18", "MP16"),
  "Immune Infiltrating" = c("MP15")
)

STATE_COLORS <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "Immune Infiltrating" = "#377EB8",
  "3CA_EMT_and_Protein_maturation" = "#666666",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)

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
