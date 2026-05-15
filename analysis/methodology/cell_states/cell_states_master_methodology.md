# Cell-State Master PBS Wrapper Methodology

## Goal

`cell_states_master.sh` runs the active state assignment and terminal state plotting scripts in the correct order after malignant epithelial UCell inputs have been prepared.

## Inputs

- `sn_outs/snSeq_malignant_epi.rds`
- `sn_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds`
- `sn_outs/UCell_3CA_MPs.rds`
- `sn_outs/snSeq_malignant_epi_cc_score.rds`

## Run Order

1. `analysis/cell_states/assign_states_approach_b_noreg.R`
2. `analysis/cell_states/relabel_unresolved_retained_3ca.R`
3. `analysis/cell_states/plot_state_overall_proportions.R`
4. `analysis/cell_states/plot_state_sample_abundance.R`
5. `analysis/cell_states/plot_state_hybrid_pairwise.R`

## Outputs

The wrapper itself does not create independent outputs. It orchestrates the outputs described in each script methodology file.

## Notes

The wrapper uses PBS live logging (`#PBS -koed`) and the `dmtcp` conda environment.
