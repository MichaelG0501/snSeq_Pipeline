# Unresolved Relabeling With Retained 3CA Metaprograms Methodology

## Goal

`analysis/cell_states/relabel_unresolved_retained_3ca.R` finalizes malignant epithelial state labels by relabelling only cells initially called `Unresolved` by Approach B noreg.

## Constraint

The snRNA-seq workflow must use the fixed scRef-retained 3CA metaprograms. It must not re-derive a retained 3CA set from the snRNA-seq cohort.

## Inputs

- `sn_outs/snSeq_malignant_epi.rds`
- `sn_outs/Auto_topmp_v2_noreg_states_B.rds`
- `sn_outs/Auto_topmp_v2_noreg_mp_adj.rds`
- `sn_outs/UCell_3CA_MPs.rds`
- `sn_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds`
- `sn_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds`
- `sn_outs/snSeq_malignant_epi_cc_score.rds`

## Retained 3CA Map

The fixed mapping is:

- `X3CA_mp_30.Respiration.1` -> `Classic Proliferative`
- `X3CA_mp_12.Protein.maturation` -> `3CA_EMT_and_Protein_maturation`
- `X3CA_mp_17.EMT.III` -> `3CA_EMT_and_Protein_maturation`

Cell-cycle-like 3CA programs are excluded before top-program assignment.

## Method

1. Align cells across primary states, adjusted scRef scores, 3CA scores, scRef UCell scores, and cell-cycle score.
2. Identify cells whose primary state is `Unresolved`.
3. For those cells, identify the top non-cell-cycle 3CA metaprogram by raw UCell score.
4. Relabel only cells whose top 3CA MP appears in the fixed retained map.
5. Leave all other unresolved cells as `Unresolved`.
6. Save `Auto_final_states.rds` as the final downstream state label vector.
7. Render updated state proportions, cell-cycle boxplots, and a combined scRef/3CA heatmap.

## Cache / Fast Plotting

Use:

```bash
Rscript analysis/cell_states/relabel_unresolved_retained_3ca.R plot_only=true
```

This reloads:

- `sn_outs/cell_states/unresolved_relabel/intermediate/unresolved_relabel_cache.rds`

and regenerates figures/tables without recomputing relabel assignments.

## Outputs

- `sn_outs/Auto_final_states.rds`
- `sn_outs/task4_unresolved_states/Auto_task4_unresolved_relabel_states.rds`
- `sn_outs/task4_unresolved_states/Auto_task4_unresolved_relabel_mp_coverage.csv`
- `sn_outs/task4_unresolved_states/Auto_task4_unresolved_relabel_proportion.pdf`
- `sn_outs/task4_unresolved_states/Auto_task4_unresolved_relabel_cc_boxplot.pdf`
- `sn_outs/task4_unresolved_states/Auto_task4_unresolved_relabel_heatmap.pdf`
- `sn_outs/task4_unresolved_states/Auto_task4_unresolved_relabel_summary.csv`
- Tiered cache/table/log outputs under `sn_outs/cell_states/unresolved_relabel/`

## Downstream Dependencies

Final state plotting scripts should prefer `Auto_final_states.rds` and fall back to `Auto_topmp_v2_noreg_states_B.rds` only when the final file has not yet been generated.
