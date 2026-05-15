# Sample-Level State Abundance Methodology

## Goal

`analysis/cell_states/plot_state_sample_abundance.R` generates sample-level abundance plots for final states and metaprogram-level views.

## Inputs

- `sn_outs/snSeq_malignant_epi.rds`
- `sn_outs/Auto_final_states.rds` or primary Approach-B states
- `sn_outs/Auto_topmp_v2_noreg_mp_adj.rds`
- `sn_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds`
- `sn_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds`

## Method

1. Load final state labels when available.
2. Align state labels, adjusted scRef scores, raw UCell scores, and merged metadata.
3. Derive:
   - top non-cell-cycle MP per cell
   - top MP including cell-cycle MPs
   - final state per cell
   - Basal-to-Intestinal internal MP breakdown
4. Build per-sample proportions for each view.
5. Generate sample orders by diversity-style ranking and by `reference_batch`.
6. Render the multi-page abundance PDF.
7. Save summary tables for reproducibility.

## Outputs

- `sn_outs/task3_sample_abundance/Auto_task3_sample_abundance.pdf`
- `sn_outs/task3_sample_abundance/Auto_sample_abundance_summary.csv`
- `sn_outs/cell_states/sample_abundance/tables/sample_abundance_summary.csv`
- `sn_outs/cell_states/sample_abundance/logs/plot_state_sample_abundance.log`

## Downstream Dependencies

None. This is a terminal figure-generation script.
