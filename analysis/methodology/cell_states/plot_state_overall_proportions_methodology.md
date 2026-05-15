# Overall State Proportions Methodology

## Goal

`analysis/cell_states/plot_state_overall_proportions.R` renders an overall composition barplot for malignant epithelial states.

## Inputs

- `sn_outs/snSeq_malignant_epi.rds`
- Preferred: `sn_outs/Auto_final_states.rds`
- Fallback: `sn_outs/Auto_topmp_v2_noreg_states_B.rds`

## Method

1. Load preferred final state labels.
2. Align labels to cells in the merged malignant epithelial object.
3. Count cells per state.
4. Convert counts to percentages of all malignant epithelial cells.
5. Apply shared state order and shared state colors from `analysis/lib/config.R`.
6. Render a presentation-sized stacked barplot.

## Outputs

- `sn_outs/Auto_overall_state_proportions.pdf`
- `sn_outs/cell_states/overall_proportions/tables/overall_state_proportions.csv`
- `sn_outs/cell_states/overall_proportions/logs/plot_state_overall_proportions.log`

## Downstream Dependencies

None. This is a terminal figure-generation script.
