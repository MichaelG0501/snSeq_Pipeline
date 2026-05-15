# Approach-B noreg State Assignment Methodology

## Goal

`analysis/cell_states/assign_states_approach_b_noreg.R` implements the selected primary malignant epithelial state definition for this snRNA-seq workspace: **Approach B, noreg**.

## Inputs

- `sn_outs/snSeq_malignant_epi.rds`
- `sn_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds`
- `sn_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds`
- `sn_outs/snSeq_malignant_epi_cc_score.rds`

## Shared Constants

Constants are sourced from `analysis/lib/config.R`:

- `STATE_GROUPS`
- `STATE_COLORS`
- `CC_MPS`
- `STATE_THRESHOLDS`
- preferred method label: `Approach B noreg`

## Method

1. Align malignant epithelial cells across the merged object, UCell score matrix, and cell-cycle score vector.
2. Drop bad scRef MPs based on negative silhouette in the retained scRef metaprogram object.
3. Split retained MPs into cell-cycle MPs (`MP1`, `MP7`, `MP9`) and non-cell-cycle MPs.
4. Normalize non-cell-cycle UCell scores by:
   - centering within `orig.ident`
   - dividing by the per-`reference_batch` standard deviation
5. Calculate the maximum adjusted score within each state group.
6. Assign each cell to the state group with the highest score.
7. Set a cell to `Unresolved` when the best group score is below `0.5`.
8. Set a cell to `Hybrid` when the gap between the top two group scores is below `0.3`.
9. Save canonical downstream RDS outputs and a cache for fast plot-only reruns.
10. Render the state heatmap, state proportion figure, and cell-cycle score boxplot.

## Cache / Fast Plotting

Use:

```bash
Rscript analysis/cell_states/assign_states_approach_b_noreg.R plot_only=true
```

This reloads:

- `sn_outs/cell_states/state_assignment/intermediate/approach_b_noreg_cache.rds`

and regenerates figures without recomputing adjusted score matrices.

## Outputs

- `sn_outs/Auto_topmp_v2_noreg_states_B.rds`
- `sn_outs/Auto_topmp_v2_noreg_mp_adj.rds`
- `sn_outs/Auto_topmp_v2_noreg_group_max.rds`
- `sn_outs/Auto_topmp_v2_noreg_heatmap_B_cconly.pdf`
- `sn_outs/Auto_topmp_v2_noreg_proportion_B_withpie.pdf`
- `sn_outs/Auto_topmp_v2_noreg_ccscore_boxplot_B.pdf`
- `sn_outs/Auto_topmp_v2_noreg_summary.csv`
- Tiered cache/table/log outputs under `sn_outs/cell_states/state_assignment/`

## Downstream Dependencies

`relabel_unresolved_retained_3ca.R`, state abundance plots, and hybrid plots consume these outputs.
