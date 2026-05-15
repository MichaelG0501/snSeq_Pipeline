# Malignancy Integrity Verification Methodology

## Goal

`analysis/infercna/verify_malignancy_integrity.R` audits step-6 malignancy outputs before malignant epithelial merge and UCell scoring.

## Inputs

- `sn_outs/sample_manifest.csv`
- Per-sample `expr_filter_status.txt`
- Per-sample `infercna_status.txt`
- Per-sample `malignancy_status.txt`
- Per-sample `<sample>_epi_f.rds`
- Per-sample `<sample>_malignancy_summary.csv`
- Per-sample `<sample>_malignancy_cluster_summary.csv`

## Method

1. Select samples with all required statuses equal to `ok`.
2. Require final `_epi_f.rds`, sample summary, and cluster summary files.
3. Verify required metadata columns:
   - `orig.ident`
   - `reference_batch`
   - `classification`
   - `malignant_clus`
   - `malignancy`
4. Confirm object-level, summary-level, and cluster-level cell counts agree.
5. Confirm each final object contains a single sample and single reference batch.
6. Stop if any sample fails integrity checks.

## Outputs

- `sn_outs/Auto_malignancy_integrity_by_sample.csv`
- `sn_outs/Auto_malignancy_integrity_summary.csv`
- Tiered table/log copies under `sn_outs/infercna/malignancy_integrity/`

## Downstream Dependencies

The audit supports safe execution of `prepare_malignant_epithelial_ucell.R`.
