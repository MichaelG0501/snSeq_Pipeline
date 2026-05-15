# Expression-Filter Summary Methodology

## Goal

`analysis/summary/expression_filter_summary.R` validates and summarizes step-4 expression filtering and the final merged post-filter object.

## Inputs

- `sn_outs/filtering_summary_sn.csv`
- `sn_outs/filtered_sample_summary.csv`
- `sn_outs/expression_filter_reason_overall.csv`
- `sn_outs/snSeq_merged.rds`
- `sn_outs/names_tmdata_sn.txt`

## Integrity Checks

The script stops if:

- `snSeq_merged.rds` contains legacy broad cell types that should not survive final annotation/filtering.
- The sum of `filtered_sample_summary.csv$cells_after` does not equal the merged object cell count.
- Required inputs are missing.

## Method

1. Read pre-filter and post-filter sample summary tables.
2. Recover source batch and technology metadata from `names_tmdata_sn.txt`.
3. Build total-cell filtering summaries.
4. Summarize expression-filter failure reasons.
5. Summarize cell-type composition overall and by technology/reference batch.
6. Generate UMAP panels for final cell type and technology.
7. Render sample-level retention figures.
8. Combine major figures into a PDF report.

## Outputs

- `sn_outs/summary/expression_filter/figures/*.png`
- `sn_outs/summary/expression_filter/reports/summary_overview.pdf`
- `sn_outs/summary/expression_filter/logs/expression_filter_summary.log`

## Downstream Dependencies

None. This is a terminal summary/reporting script.
