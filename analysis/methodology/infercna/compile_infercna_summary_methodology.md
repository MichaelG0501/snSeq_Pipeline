# InferCNA Summary Compilation Methodology

## Goal

`analysis/infercna/compile_infercna_summary.R` compiles per-sample step-5 InferCNA outputs into cohort-level status, sample, cluster, and signature summary tables.

## Inputs

- `sn_outs/sample_manifest.csv`
- `sn_outs/filtered_sample_summary.csv`
- `sn_outs/by_samples/<sample>/expr_filter_status.txt`
- `sn_outs/by_samples/<sample>/infercna_status.txt`
- `sn_outs/by_samples/<sample>/<sample>_infercna_summary.csv`
- `sn_outs/by_samples/<sample>/<sample>_infercna_cluster_summary.csv`
- `sn_outs/by_samples/<sample>/<sample>_epi.rds`
- `sn_outs/by_samples/<sample>/<sample>_signatures.rds`

## Method

1. Read the manifest and expression-filter summary.
2. Derive effective InferCNA status for every sample.
3. Flag `ok_missing_outputs` when status says `ok` but required output files are absent.
4. Bind per-sample InferCNA summaries for expression-filter `ok` samples.
5. Bind per-cluster InferCNA summaries.
6. Compile tumour signature gene membership from per-sample signature RDS files.
7. Write status tables, cohort summaries, and signature reference files.

## Outputs

- `sn_outs/infercna_status_by_sample.csv`
- `sn_outs/infercna_sample_summary.csv`
- `sn_outs/infercna_overall_summary.csv`
- `sn_outs/infercna_classification_overall.csv`
- `sn_outs/infercna_cluster_label_overall.csv`
- `sn_outs/infercna_cluster_summary.csv`
- `sn_outs/cancer_signatures_by_sample.csv`
- `sn_outs/cancer_signatures_summary.csv`
- `sn_outs/cancer_signatures.txt`
- Run log under `sn_outs/infercna/compile_summary/logs/`

## Downstream Dependencies

`validate_malignancy_signatures.R` consumes `cancer_signatures.txt` and `cancer_signatures_summary.csv`.
