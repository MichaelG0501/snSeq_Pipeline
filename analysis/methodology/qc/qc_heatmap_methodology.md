# QC Heatmap Methodology

## Goal

`analysis/qc/qc_heatmap.R` generates cohort-level QC heatmaps before and after expression filtering using the current snRNA-seq marker and QC rules.

## Inputs

- `sn_outs/sample_manifest.csv`
- `sn_outs/filtered_sample_summary.csv`
- `sn_outs/expression_filter_markers.csv`
- Per-sample post-QC objects: `sn_outs/by_samples/<sample>/<sample>.rds`
- Per-sample annotated/filtering objects: `sn_outs/by_samples/<sample>/<sample>_anno.rds`

## Method

1. Load the expression-filter marker panel and housekeeping genes.
2. Read post-QC objects for the pre-expression-filter heatmap.
3. Read annotated objects and retain final kept cells for the final heatmap.
4. Extract the required RNA layers while padding missing marker genes with zeros.
5. Compute marker, housekeeping, mitochondrial, and metadata annotations.
6. Render ComplexHeatmap reports with large presentation-facing dimensions.
7. Write stage/cell-type count summaries.

## Outputs

Compatibility outputs:

- `sn_outs/Auto_QC_snSeq_prefilter.png`
- `sn_outs/Auto_QC_snSeq_final.png`
- `sn_outs/Auto_QC_snSeq_heatmaps.pdf`
- `sn_outs/Auto_qc_heatmap_stage_summary.csv`

Tiered copies:

- `sn_outs/qc/heatmaps/figures/`
- `sn_outs/qc/heatmaps/reports/`
- `sn_outs/qc/heatmaps/tables/`
- `sn_outs/qc/heatmaps/logs/`

## Downstream Dependencies

None. This is a terminal QC review script.
