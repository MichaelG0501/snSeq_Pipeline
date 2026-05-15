# Malignancy Summary Plotting Methodology

## Goal

`analysis/infercna/plot_malignancy_summary.R` summarizes and plots final step-6 malignancy calls across samples and reference batches.

## Inputs

- `sn_outs/sample_manifest.csv`
- `sn_outs/filtered_sample_summary.csv`
- `sn_outs/infercna_status_by_sample.csv`
- Per-sample malignancy status and summary files
- Per-sample final `_epi_f.rds` objects

## Method

1. Derive effective malignancy status for every manifest sample.
2. Bind sample-level and cluster-level malignancy summaries.
3. Count final malignancy labels from final epithelial objects.
4. Summarize malignancy composition overall, by sample, and by reference batch.
5. Render presentation-facing barplots and a combined PDF report.

## Outputs

- Root compatibility tables under `sn_outs/malignancy_*.csv`
- Figures under `sn_outs/infercna/malignancy_summary/figures/`
- Report under `sn_outs/infercna/malignancy_summary/reports/`
- Run log under `sn_outs/infercna/malignancy_summary/logs/`

## Downstream Dependencies

None. This is a terminal figure/reporting script.
