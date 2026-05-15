# Pairwise Hybrid State Plot Methodology

## Goal

`analysis/cell_states/plot_state_hybrid_pairwise.R` summarizes hybrid cells as pairwise relationships among resolved state groups.

## Inputs

- `sn_outs/Auto_final_states.rds` or `sn_outs/Auto_topmp_v2_noreg_states_B.rds`
- `sn_outs/Auto_topmp_v2_noreg_mp_adj.rds`

## Method

1. Load final state labels and adjusted MP scores.
2. Calculate state-group maxima from adjusted MP scores.
3. Select cells labelled `Hybrid`.
4. For each hybrid cell, assign the two highest-scoring resolved state groups.
5. Count pairwise hybrid edges and resolved-state node sizes.
6. Render:
   - node/edge network plot
   - pairwise heatmap
7. Save the combined table of node and pairwise hybrid counts.

## Outputs

- `sn_outs/task6_hybrid_pairwise/Auto_task6_hybrid_pairwise_nodeplot_noreg.pdf`
- `sn_outs/task6_hybrid_pairwise/Auto_task6_hybrid_pairwise_heatmap_noreg.pdf`
- `sn_outs/task6_hybrid_pairwise/Auto_task6_hybrid_pairwise_nodeplot_summary.csv`
- `sn_outs/cell_states/hybrid_pairwise/tables/hybrid_pairwise_summary.csv`
- `sn_outs/cell_states/hybrid_pairwise/logs/plot_state_hybrid_pairwise.log`

## Downstream Dependencies

None. This is a terminal figure-generation script.
