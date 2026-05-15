# snSeq_Pipeline Analysis Map

Updated: 2026-05-14

This map describes the current `analysis/` layout, script responsibilities, dependencies, output tiers, and run order. New analysis scripts must update this file and add a matching methodology document under `analysis/methodology/<analysis-subfolder>/`.

## Current Conventions

- Preferred malignant epithelial state definition: **Approach B, noreg**.
- Final state labels: `sn_outs/Auto_final_states.rds`.
- State normalization batch field: `reference_batch`.
- Shared constants and reusable helpers: `analysis/lib/`.
- Canonical downstream compatibility outputs remain in `sn_outs/` where existing pipeline stages depend on them.
- New terminal and review outputs should also be mirrored into tiered folders:
  - `intermediate/`
  - `tables/`
  - `figures/`
  - `logs/`
  - `reports/`

## Shared Library

| Path | Purpose |
|---|---|
| `analysis/lib/config.R` | Project paths, state-method selection, shared thresholds, metadata column names, state groups/colors, retained 3CA map, plot defaults, output-tier helpers, argument parsing. |
| `analysis/lib/io_helpers.R` | Status-file readers, optional CSV readers, required-file checks, checked CSV/RDS writers, Seurat layer matrix extraction. |
| `analysis/lib/logging.R` | Lightweight run-summary logging with start/end time, inputs, outputs, parameters, cache reuse, and session info. |
| `analysis/lib/state_helpers.R` | Shared state logic: sample/reference_batch normalization, 3CA label cleaning, MP labelling, state ordering/colors, group-max calculation, final-state loading. |
| `analysis/lib/plot_helpers.R` | Presentation-oriented ggplot and heatmap text defaults. |

Methodology: `analysis/methodology/lib/shared_analysis_library_methodology.md`.

## Run Order

| Order | Script | Depends On | Produces | Downstream Consumers |
|---:|---|---|---|---|
| 1 | `analysis/summary/expression_filter_summary.R` | Step 4 outputs: `filtering_summary_sn.csv`, `filtered_sample_summary.csv`, `expression_filter_reason_overall.csv`, `snSeq_merged.rds` | Step-4 summary figures/report under `sn_outs/summary/expression_filter/` | Terminal only |
| 2 | `analysis/qc/qc_heatmap.R` | Step 1 and step 4 sample objects, `sample_manifest.csv`, `expression_filter_markers.csv` | QC heatmaps under root compatibility outputs and `sn_outs/qc/heatmaps/` | Terminal only |
| 3 | `analysis/annotation/plot_annotation_featureplots.R` | Step 2/3 annotated sample objects, 3CA marker workbook, annotation/filter marker CSVs | Annotation review PDF and status tables | Terminal only |
| 4 | `analysis/infercna/compile_infercna_summary.R` | Step 5 InferCNA per-sample summaries/statuses/signature RDS files | InferCNA cohort summary tables and `cancer_signatures.*` files | `validate_malignancy_signatures.R`, `plot_malignancy_summary.R` |
| 5 | `analysis/infercna/validate_malignancy_signatures.R` | `cancer_signatures.txt`, `cancer_signatures_summary.csv`, step-5 `_epi.rds` objects, cell-cycle reference CSV | Signature validation PDF | Terminal audit |
| 6 | `analysis/infercna/verify_malignancy_integrity.R` | Step 6 `_epi_f.rds`, malignancy summaries, status files | Malignancy integrity CSVs | `prepare_malignant_epithelial_ucell.R` gate/review |
| 7 | `analysis/infercna/plot_malignancy_summary.R` | Step 6 summaries plus compiled InferCNA status | Malignancy summary tables and figures | Terminal only |
| 8 | `analysis/infercna/plot_infercna_scatter_examples.R` | Selected step-5 `_outs.rds`, `_epi.rds`, InferCNA summaries | Example CNA scatter figure | Terminal only |
| 9 | `analysis/infercna/plot_merged_cnv_profile.R` | Final `_epi_f.rds` objects and `_outs.rds` files | Merged CNA profile PDF | Terminal only |
| 10 | `analysis/cell_states/prepare_malignant_epithelial_ucell.R` | Step 6 complete samples, scRef MP definitions, 3CA MP CSV, cell-cycle reference CSV | `snSeq_malignant_epi.rds`, UCell score RDS files, cell-cycle score RDS files | State assignment |
| 11 | `analysis/cell_states/assign_states_approach_b_noreg.R` | Merged malignant epithelial object, scRef UCell scores, MP definitions, CC score | `Auto_topmp_v2_noreg_states_B.rds`, adjusted MP matrix, group max matrix, Approach-B figures | Unresolved relabel, terminal state plots |
| 12 | `analysis/cell_states/relabel_unresolved_retained_3ca.R` | Approach-B noreg labels, adjusted MP matrix, 3CA UCell scores, retained scRef 3CA map | `Auto_final_states.rds`, unresolved relabel tables/figures | Final state plotting and downstream state analyses |
| 13 | `analysis/cell_states/plot_state_overall_proportions.R` | `Auto_final_states.rds` or primary Approach-B labels | Overall state proportion PDF/table | Terminal only |
| 14 | `analysis/cell_states/plot_state_sample_abundance.R` | Final state labels, adjusted MP matrix, UCell scores, MP definitions | Per-sample abundance PDF/table | Terminal only |
| 15 | `analysis/cell_states/plot_state_hybrid_pairwise.R` | Final/primary state labels and adjusted MP matrix | Pairwise hybrid nodeplot/heatmap/table | Terminal only |

`cell_states_master.sh` runs scripts 11-15 sequentially after script 10 has produced UCell inputs.

## Legacy And Superseded Scripts

| Script | Status | Why Kept | Downstream Consumers |
|---|---|---|---|
| `analysis/legacy/cell_states/legacy_prefix_unresolved_pan_cancer_noreg_comparison.R` | Legacy comparison only | Preserves the older all-pan-cancer-3CA unresolved-cell comparison plots. It must not write root-level downstream files. | None |

## Outdated References Fixed

- `8_geneNMF.sh` previously pointed to missing `geneNMF.R`; it now runs `analysis/cell_states/prepare_malignant_epithelial_ucell.R`.
- Step-4 summary execution now points to `analysis/summary/expression_filter_summary.R`.
- The old all-pan-cancer unresolved comparison was moved under `analysis/legacy/cell_states/` and no longer writes `Auto_final_states.rds` or other downstream files.

## Terminal Figure-Generation Scripts

- `analysis/summary/expression_filter_summary.R`
- `analysis/qc/qc_heatmap.R`
- `analysis/annotation/plot_annotation_featureplots.R`
- `analysis/infercna/validate_malignancy_signatures.R`
- `analysis/infercna/plot_malignancy_summary.R`
- `analysis/infercna/plot_infercna_scatter_examples.R`
- `analysis/infercna/plot_merged_cnv_profile.R`
- `analysis/cell_states/plot_state_overall_proportions.R`
- `analysis/cell_states/plot_state_sample_abundance.R`
- `analysis/cell_states/plot_state_hybrid_pairwise.R`
- `analysis/legacy/cell_states/legacy_prefix_unresolved_pan_cancer_noreg_comparison.R`

## Cached Visualization Reuse

Scripts with long-running computation expose intermediate/cache paths in their headers and logs. Current cache-enabled scripts:

- `analysis/cell_states/assign_states_approach_b_noreg.R plot_only=true`
- `analysis/cell_states/relabel_unresolved_retained_3ca.R plot_only=true`

Future long-running analysis scripts must save computational results under `intermediate/` so plots can be regenerated without rerunning the heavy compute step.
