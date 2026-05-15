# Malignant Epithelial Merge and UCell Scoring Methodology

## Goal

`analysis/cell_states/prepare_malignant_epithelial_ucell.R` prepares the malignant epithelial object used by all downstream state analyses. It merges only completed step-6 samples and computes UCell scores for the retained scRef metaprograms and the 3CA metaprograms.

## Inputs

- `sn_outs/sample_manifest.csv`
- Per-sample status files under `sn_outs/by_samples/<sample>/`
- Final step-6 epithelial objects: `<sample>_epi_f.rds`
- scRef metaprograms:
  `/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds`
- 3CA metaprograms:
  `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv`
- Cell-cycle reference:
  `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv`

## Sample Eligibility

A sample is included only if all of the following are true:

- `expr_filter_status.txt = ok`
- `infercna_status.txt = ok`
- `malignancy_status.txt = ok`
- `<sample>_epi_f.rds` exists

Only cells labelled `malignant_level_1` or `malignant_level_2` are retained.

## Method

1. Read the current manifest and status files.
2. Load eligible final epithelial objects.
3. Subset each object to malignant epithelial cells.
4. Merge all malignant epithelial objects into `snSeq_malignant_epi.rds`.
5. Copy/read the retained scRef metaprogram definitions.
6. Build an explicit sparse expression matrix from Seurat v5 `data.*` layers for the genes needed by downstream scoring.
7. Score retained scRef MPs with `ScoreSignatures_UCell()`.
8. Score 3CA metaprograms with `ScoreSignatures_UCell()`.
9. Derive the top 50 consensus cell-cycle genes by mean expression and save a per-cell cell-cycle score.

## Outputs

- `sn_outs/snSeq_malignant_epi.rds`
- `sn_outs/snSeq_malignant_epi_meta.rds`
- `sn_outs/snSeq_malignant_epi_sample_summary.csv`
- `sn_outs/snSeq_malignant_epi_status_input.csv`
- `sn_outs/snSeq_malignant_epi_overall_summary.csv`
- `sn_outs/snSeq_malignant_epi_cc_top50.rds`
- `sn_outs/snSeq_malignant_epi_cc_score.rds`
- `sn_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds`
- `sn_outs/UCell_3CA_MPs.rds`
- Run log under `sn_outs/cell_states/ucell_scoring/logs/`

## Downstream Dependencies

The output RDS files are required by `assign_states_approach_b_noreg.R` and `relabel_unresolved_retained_3ca.R`.
