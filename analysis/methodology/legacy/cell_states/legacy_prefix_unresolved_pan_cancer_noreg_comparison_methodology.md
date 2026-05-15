# Legacy Unresolved Pan-Cancer Comparison Methodology

## Status

`analysis/legacy/cell_states/legacy_prefix_unresolved_pan_cancer_noreg_comparison.R` is **legacy comparison only**. It must not be used to define final states and must not write files consumed by downstream scripts.

## Superseded By

- `analysis/cell_states/relabel_unresolved_retained_3ca.R`

## Inputs

- `sn_outs/snSeq_malignant_epi.rds`
- `sn_outs/Auto_topmp_v2_noreg_states_B.rds`
- `sn_outs/UCell_3CA_MPs.rds`
- `sn_outs/snSeq_malignant_epi_cc_score.rds`

## Legacy Method

1. Select cells labelled `Unresolved` by the primary Approach-B noreg state call.
2. Exclude cell-cycle-like 3CA metaprograms.
3. Assign each unresolved cell to its highest raw 3CA UCell score.
4. Plot raw top-3CA assignments and an unsupervised heatmap for comparison.

## Difference From Current Method

The current method uses only the fixed scRef-retained 3CA map. This legacy method uses all non-cell-cycle 3CA programs and can therefore produce labels that are not part of the current finalized state convention.

## Outputs

All outputs are isolated under:

- `sn_outs/legacy/unresolved_pan_cancer_noreg/`

No output is written to the root of `sn_outs/`.
