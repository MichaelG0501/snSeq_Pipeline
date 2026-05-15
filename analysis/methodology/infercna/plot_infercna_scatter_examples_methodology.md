# InferCNA Scatter Example Methodology

## Goal

`analysis/infercna/plot_infercna_scatter_examples.R` creates a three-panel InferCNA scatter comparison for representative high, mixed, and low CNA signal samples.

## Inputs

For each configured example sample:

- `<sample>_outs.rds`
- `<sample>_epi.rds`
- `<sample>_infercna_summary.csv`

## Method

1. Load InferCNA output matrices and thresholds.
2. Convert InferCNA coordinates with `cnaScatterPlot()`.
3. Join cell-level CNA classifications from the epithelial object.
4. Plot CNA correlation versus CNA signal.
5. Overlay threshold lines from the sample summary.
6. Combine the three examples into one presentation-sized PNG.

## Outputs

- `sn_outs/infercna/scatter_examples/figures/malignancy_scatter_comparison.png`
- `sn_outs/infercna/scatter_examples/logs/plot_infercna_scatter_examples.log`

## Downstream Dependencies

None. This is a terminal review figure.
