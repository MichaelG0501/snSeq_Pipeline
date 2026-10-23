# snSeq_Pipeline analysis map

Updated: 2026-10-19. Canonical project and outputs are in live storage. This map
supersedes the historical uncentred nMP19/3CA-relabel route. Stages 1–6 regenerate
from raw counts and preserve their existing QC, annotation, filtering, reference,
and malignancy rules. Current state assignment projects fixed current centred
refined scRef MPs using Approach B, noreg, and reference_batch normalization.

## Execution and canonical outputs

`1_QC_Pipeline.sh` performs raw QC. `Auto_run_standard_pipeline.sh` is the tracked
PBS orchestration entry after successful QC. It submits standard stages 2–6,
compiles current InferCNA signatures, audits malignancy integrity, runs step 7
scoring and step 8 states, then submits the updated PDO matched-pair comparison.
It throttles against the user's total PBS jobs (maximum 46), waits for each
stage's scheduler results, and checks outputs/status markers before proceeding.
Only step-5 explicit no_epi/no_cell biological skips may have nonzero exits.
Missing/insufficient reference, unexpected errors, and missing matched samples
stop the run. Controller status and job ledger are in sn_outs/logs and reports.

Current downstream inputs:

| Path under live sn_outs/ | Producer | Consumers |
|---|---|---|
| snSeq_malignant_epi.rds; snSeq_malignant_epi_meta.rds | prepare_malignant_epithelial_ucell.R | Current state assignment; external PDO comparison uses metadata |
| cell_states/centred/intermediate/merged_refined_mp_genes.rds | step 7, copied from current scRef | External PDO comparison |
| cell_states/centred/intermediate/merged_refined_ucell_scores.rds | step 7, full-gene ranking | Current state assignment and PDO comparison |
| cell_states/centred/tables/centred_refined_mp_state_grouping.csv | step 7, copied from current scRef | Current state assignment and PDO comparison |
| cell_states/centred/intermediate/centred_refined_noreg_states.rds | step 8 | External PDO comparison |
| cell_states/centred/intermediate/centred_refined_noreg_mp_adj.rds; centred_refined_noreg_group_max.rds | step 8 | External PDO comparison; reproducible diagnostics |
| cell_states/centred/tables/centred_refined_noreg_cell_assignments.csv; centred_refined_noreg_sample_state_abundance.csv | step 8 | Terminal source tables |
| cell_states/centred/figures/Auto_final/01_centred_state_summary.pdf | step 8 | Terminal cohort figures |

Current reference is the live scRef centred mp_refinement/intermediate gene list
and tables/centred_refined_mp_state_grouping.csv. Step 7 copies both and records
provenance. The reference's 17 final MPs are not filtered again in snSeq. Five
biological states plus Unresolved/Hybrid use the reference names, not PDO names.
No current state is relabelled by 3CA. Diagnostic UCell_3CA_MPs and CC scores
remain available but never modify current centred labels.

## Maintained scripts and run order

All paths below are relative to analysis/ unless stated otherwise.

| Order | Script | Dependency | Status / downstream use |
|---|---|---|---|
| 1 | summary/expression_filter_summary.R | Step-4 merged object and filtering summaries | Active terminal; invoked by step-4 wrapper |
| 2 | qc/qc_heatmap.R | Step-1/4 objects, manifest, filter markers | Optional terminal QC |
| 3 | annotation/plot_annotation_featureplots.R | Step-2/3 objects, markers | Optional terminal annotation review |
| 4 | infercna/compile_infercna_summary.R | Step-5 status, outputs, signatures | Active upstream/terminal; cancer_signatures.txt feeds step 6 |
| 5 | infercna/validate_malignancy_signatures.R | Compiled signatures and epithelial objects | Optional terminal audit |
| 6 | infercna/verify_malignancy_integrity.R | Step-6 objects and summaries | Required audit before malignant merge |
| 7 | infercna/plot_malignancy_summary.R | Step-6 and compiled summaries | Optional terminal review |
| 8 | infercna/plot_infercna_scatter_examples.R | Selected step-5 objects | Optional terminal review |
| 9 | infercna/plot_merged_cnv_profile.R | Final epithelial/CNA objects | Optional terminal review |
| 10 | cell_states/prepare_malignant_epithelial_ucell.R | Manifest, three ok statuses, final epithelial objects, current scRef gene/group reference, diagnostic 3CA/CC genes | Active current step 7; merged malignant object and full-gene UCell |
| 11 | cell_states/assign_states_approach_b_noreg.R | Step-7 metadata, scores, grouping | Active current step 8; current centred states, matrices, tables, figures |
| 12 (external) | PDOs_Pipeline/analysis/cell_states/legacy_Auto_pdo_sn_matched_pair_comparison.R | Current PDO centred inputs and current snSeq outputs above | Active terminal despite historical filename; independent native state/MP proportions for two matched patients |

Current entry points are 8_geneNMF.sh and cell_states_master.sh. New analysis-only
launchers live under job_scripts/ with matching relative paths and basenames;
existing pipeline entry points remain in place. All use explicit live-root PBS
log destinations. Analysis source and complex methods remain tracked in Git.

## Shared library

lib/config.R contains live paths, current centred reference/output paths, current
colours, thresholds, and historical constants for explicitly legacy scripts.
lib/state_helpers.R supplies exact sample-centre/batch-SD normalization and current
Approach-B threshold calls. lib/io_helpers.R supplies checked I/O and layer
assembly; lib/logging.R supplies session/run logs; lib/plot_helpers.R supplies
presentation defaults. Method: methodology/lib/shared_analysis_library_methodology.md.

## Historical scripts retained without deletion or renaming

These scripts describe the old uncentred/3CA-relabel route and have no current
downstream consumers. They must not be run by the current state entry point:

- cell_states/relabel_unresolved_retained_3ca.R
- cell_states/plot_state_overall_proportions.R
- cell_states/plot_state_sample_abundance.R
- cell_states/plot_state_hybrid_pairwise.R
- legacy/cell_states/legacy_prefix_unresolved_pan_cancer_noreg_comparison.R

Historical Auto_final_states.rds, Auto_topmp_v2_noreg_* files, uncentred nMP19 gene
objects, and UCell_nMP19_filtered.rds are not inputs to the current workflow.
Historical methodology is retained as provenance and is not current run guidance.

## Caches and validation

Step-7 UCell always rebuilds. Step-8 plot_only=true requires unchanged metadata,
score and grouping MD5 fingerprints. PDO_REPLOT_ONLY=1 validates its recorded
input fingerprints. Smoke tests and probes remain in ignored Auto_session_checks/
and never become scientific inputs. Final figures, source tables, named states,
reference copies, and matrices are canonical live files. Historical root-level
compatibility files persist only where core pipeline consumers still need them.

## Final visualization additions (2026-10-20)

| Script | Status | Dependencies | Outputs / consumers |
|---|---|---|---|
| summary/plot_final_uncorrected_umaps.R | Active terminal | Stage-4 snSeq_merged.rds with its uncorrected full-cell embedding; final malignant epithelial merge and current states | summary/final_uncorrected_umaps/: full coordinates, PCA/UMAP cache and two-page figures/Auto_final/01_full_uncorrected_umaps.pdf; terminal/replot inputs |
| infercna/compare_matched_pdo_sn_cna.R | Active terminal | Current matched biopsy CNA/statuses/malignant metadata/states; validated native PDO CNA and subclone caches; maintained PDO calling method; hg38 gene order | infercna/matched_pdo_sn/: full cell/bin/clone tables, native calls and binned cache, figures/Auto_final/01_matched_cna_subclones.pdf; no active downstream consumer |
| lib/cna_helpers.R | Active shared support | Validated cached PDO method and maintained driver | Isolated exact calling environment and shared-bin/plot-cell helpers; consumed by matched CNA comparison |

Both scripts have mirrored methodology and matching local job_scripts launchers.
The full stage-4 cohort UMAP already exists without batch correction. The new
malignant map uses every final malignant epithelial cell and uncorrected PCA,
without regression/integration. CNA figures preserve native references; clone
labels are independently defined, and small unassigned groups stay in denominators.

## Matched PDO / biopsy state concordance

| Script | Status / run order | Dependencies | Outputs / consumers |
|---|---|---|---|
| cell_states/compare_matched_pdo_sn_state_concordance.R | Active terminal; after local stages 6–8 and current PDO state assignment | Matched RNA counts/states; verified current scRef MP genes/grouping; current PDO MPs; manifest/statuses; exact original PDO/scRef effect-gene pool | cell_states/matched_pdo_sn_concordance/: three 10-page PDFs (SUR680, SUR791, pooled), frozen reference/available genes, scores/pseudobulks/source tables and relative-linked reports; no active downstream consumer |
| lib/concordance_helpers.R | Active shared support | Native RNA count layers, grouped MP definitions and matched state pseudobulks | Validated counts, descriptive within-sample effects, paired summaries and conditional overlap/specificity functions; consumed by matched state concordance |

The reference evidence layers are preserved, with descriptive effect/marker
contrasts replacing replicate-based DGE inference for single-sample patients.
The frozen original 3,538-gene pool is restricted only to assayed genes shared
by the matched samples (currently 3,268), without new informative-gene selection.
Pooled effects use identical matched patients and equal patient weights; UCell
means and projected proportions pool cells as in the reference. The original
page order is restored. Page 4 projects individual current PDO MPs onto native
sn states; pages 5–10 each compare one PDO state against all five sn states.
Rare states remain in UCell means; the 20-cell gate applies only to pseudobulk
effects. All source cells remain audited; native uncertain labels are excluded
from receiving groups and projected uncertainty stays in denominators. Current
results are stable under sn_outs/cell_states/matched_pdo_sn_concordance/
{figures,tables,reports,logs,intermediate,reference}. Current scRef source MPs are
verified identical to canonical sn copies. The launcher supports plot_only=true;
canonical native labels remain unchanged.
