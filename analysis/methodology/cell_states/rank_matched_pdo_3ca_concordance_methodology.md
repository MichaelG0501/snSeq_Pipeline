# Pan-cancer 3CA relative matched-PDO concordance

Script: `analysis/cell_states/rank_matched_pdo_3ca_concordance.R`.
Status: active terminal analysis; no downstream pipeline consumer.

## Question and prespecified endpoints

For SUR680/H_post_T1_biopsy and SUR791/L_post_T1_biopsy, does the true matched
PDO rank highest among the entire current PDO cohort when comparing pan-cancer
3CA metaprogramme activity? True PDO identities are SUR680T3_PDO and SUR791T3_PDO.
Matched identities annotate results and never choose features or parameters.
Report the observed ranks even when the true pair is not the highest match.

The primary endpoint is Spearman correlation across all 67 fixed 3CA MPs between
the mean UCell profiles of each biopsy and each PDO. This follows the dataset-mean
MP-profile comparison in the historical external
`PDOs_Pipeline/analysis/metaprograms/PDO_mp_correlation_crossdata.R`, adapted
to individual samples. It does not reuse that script's historical MP/state
definitions, analytical fallbacks, or score-magnitude significance labels.

A prespecified sensitivity excludes MP1–MP5, the supplied cell-cycle/chromatin
programmes, and uses all remaining 62 MPs. A supporting comparison conditions on
the four native-state correspondences already specified for transcriptional
effects: Classic proliferation to itself; Squamous-to-intestinal to
Columnar-to-intestinal; Glandular-to-intestinal to Glandular differentiation;
Stress-adaptive to itself. Native labels are never reassigned. Unmapped immune
mimicry is omitted from comparison figures; all native counts are exported.

## Exact inputs

Local inputs under live `snSeq_Pipeline/sn_outs/`:

- `snSeq_malignant_epi.rds` and `sample_manifest.csv`.
- `cell_states/centred/intermediate/centred_refined_noreg_states.rds`.
- Both matched samples' `expr_filter_status.txt`, `infercna_status.txt`, and
  `malignancy_status.txt`, each required to be `ok`.

Read-only external inputs:

- Live `PDOs_Pipeline/PDOs_outs/PDOs_merged.rds`.
- Its `centred_mp_refinement/centred_refined_noreg_states.rds`.
- `PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_config.R`, including the current
  exclusion of SUR843T3_PDO and independent native PDO state order.
- `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv`:
  the fixed 67-column pan-cancer gene definitions, not current scRef/PDO MPs.

Tracked local configuration, I/O, logging and count-extraction helpers are also
fingerprinted. Samples are identified exclusively through `orig.ident`; cell
names are used only to align already identified cells with scores and states.

## Scoring and missing genes

Recompute both modalities using the same installed UCell version, `maxRank=1500`,
and the identical gene vectors for every MP. Read the full assayed-gene RNA count
matrices. Positive per-cell normalization and log1p preserve gene order, so these
within-cell count ranks equal log-normalized RNA ranks. UCell never receives a
signature-only expression matrix. Existing historical score caches are not
silently mixed across different gene coverage or package versions.

Restrict each original MP once to the gene intersection assayed in all included
sources. Each cell is still ranked against its complete source assay gene
universe. This preserves full-gene ranking and identical programme definitions
without inventing expression zeros for unassayed genes. Stop if any MP has no
remaining genes; do not discard or substitute it. Export original and available
gene counts, exact retained names and every skipped gene. Differences in overall
assay coverage and platform detection remain scientific limitations.

All retained malignant cells contribute to whole-sample means, including Hybrid
and Unresolved cells. State means use cells with the unchanged native state.
Require at least 20 cells in both compared groups, for whole-sample and
state-conditioned results. Missing and constant-profile comparisons remain NA.
No patient-state versus rest transformation is applied: this endpoint compares
MP activity levels directly, as requested.

## Ranking and examples

For each query and prespecified MP set, rank every finite candidate rho in
descending order. Ties within 1e-12 receive the same minimum rank. Export the
eligible denominator and matched-minus-highest-other rho. All row maxima,
including ties, receive a red outline (linewidth 1.4); the true matched PDO
receives a black outline (0.7). No stars or page-footer paragraphs are added.

Whole-sample scatter examples compare Matched PDO, Highest matching and Median
matching. Highest matching means the highest-scoring OTHER PDO. Median matching
is a distinct other PDO whose rank is nearest the cohort median rank. Resolve
example ties alphabetically; retain every tie in source tables. Each scatter
point is one of the same 67 MPs, with the same axis limits within each page.
Label the ten MPs with the highest activity across any panel on that page,
using the same ten labels in all three panels and 13.6-point annotation text.
This affects annotation only: all 67 points enter the correlation. Heatmap and
state-rank group labels show biopsy cell counts; the unsupported proliferation
row explicitly states fewer than 20 cells.

## Selected report and interpretation

One six-page, 18-by-10-inch PDF contains the primary whole-sample heatmap, ranks
with cell-cycle sensitivity, one three-panel mean-MP scatter page per patient,
the state-conditioned heatmap, and state-conditioned ranks. Keep comparisons on
their assigned page, with readable labels, and do not export duplicate standalone
figures. Source tables retain full precision.

High correlations can reflect shared programme baselines or epithelial biology.
Whole-sample scores also reflect state composition. Neither MPs nor individual
cells supply independent patient replicates. PDO candidates include repeated
donors and treatment conditions. Ranks describe this cohort and do not establish
absolute biological equivalence, ancestry, or population-level significance.
Do not tune MP subsets, thresholds, or normalization to improve true-match ranks.

This endpoint is separate from Celligner. The cited Nature study
([doi:10.1038/s41586-026-10806-y](https://doi.org/10.1038/s41586-026-10806-y))
uses bulk expression with a pan-cancer reference, contrastive PCA, MNN alignment,
and distances in 70-PC space. No simplified correlation or unaligned PCA is
labelled Celligner. A distinct feasibility assessment or explicitly authorized
adaptation is required before applying that method to this small malignant-only
single-cell/nucleus cohort.

## Execution and stable outputs

Run from the live project root:

```bash
qsub job_scripts/cell_states/rank_matched_pdo_3ca_concordance.sh
qsub -v plot_only=true job_scripts/cell_states/rank_matched_pdo_3ca_concordance.sh
```

The minimal ignored launcher activates `dmtcp`, requests 4 CPUs, 64 GB, four
hours, and writes live PBS stdout/stderr in the project root with `-koed`.
Analytical logic is entirely in the maintained R script.

Stable output root: `sn_outs/cell_states/matched_pdo_3ca_concordance/`.

- `intermediate/pan_cancer_mp_cache.rds`: per-cell scores, metadata, signatures,
  assay gene lists, annotation and original input/parameter fingerprints.
- `intermediate/pan_cancer_mp_rank_results.rds`: profiles, queries, scores,
  ranks and selected examples for independent downstream inspection.
- `tables/`: all candidate correlations, matched ranks, complete group counts,
  gene availability, all mean profiles, and selected scatter identities.
- `reference/`: exact scoring signatures and the original signature-gene union.
- `figures/Auto_final/01_matched_PDO_3CA_relative_concordance.pdf`: selected report.
- `reports/`: current result index, page/figure indices and validation results.
- `logs/pan_cancer_mp_concordance.log`: start/end, inputs, outputs, parameters,
  cache reuse and R session information.

`plot_only=true` requires unchanged scientific input paths, sizes, modification
times and compute parameters. It reuses per-cell scores and recomputes only
lightweight aggregation, correlations and plots. Finish source edits before
submission; do not edit scripts while a running R process reads them.
