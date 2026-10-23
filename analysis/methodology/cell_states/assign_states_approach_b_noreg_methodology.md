# Current centred refined scRef projection onto malignant snRNA-seq

## Scope and references

`analysis/cell_states/assign_states_approach_b_noreg.R` consumes the malignant-only
metadata and full-gene UCell scores produced by step 7 in this workspace. It uses
only the current scRef centred refined 17-MP gene collection and grouping table,
copied from the live scRef `ref_outs/Metaprogrammes_Results/centred/mp_refinement/`
`intermediate/merged_refined_mp_genes.rds` and
`tables/centred_refined_mp_state_grouping.csv`. The reference is fixed externally;
no snRNA-specific NMF, extra MP selection, or cohort-based 3CA retention occurs.

## Normalization and assignment

Every cell is malignant_level_1 or malignant_level_2 from completed local stages
1–6. Cell identity comes from exact barcodes; samples come from `orig.ident`.
Scores are centred by their within-sample mean, then divided by the standard
deviation of the original uncentred scores within `reference_batch`. Zero or
undefined standard deviations become one, matching the current scRef noreg
implementation. Cell-cycle regression is never performed.

The grouping table supplies the five biological groups: Classic proliferation,
Squamous-to-intestinal, Glandular-to-intestinal, Stress-adaptive, and Cancer-cell
immune mimicry. MP1, MP5, and MP13+ are cell-cycle/QC rows only. Each biological
group score is the maximum of its constituent adjusted MP scores. A maximum
below 0.5 yields Unresolved. Otherwise, a best-minus-second group-score gap below
0.3 yields Hybrid; other cells receive the highest group. Ties use table order.
No 3CA relabelling modifies these current centred states.

## Integrity and outputs

Require complete, unique barcodes and all expected MPs; do not intersect away
missing cells or silently drop MPs. Require finite scores and recognized states.
Persist named states, adjusted scores (including QC-only cell-cycle rows), group
maxima, cell assignments, top adjusted MP and non-CC MP, zero-completed sample
abundances, state counts, a reference fingerprinted cache, and a run/session log
under live `sn_outs/cell_states/centred/` output tiers. Abundance counts must sum to
the complete malignant cell count. The PDO comparison reads these current files.

## Execution and display

Run `qsub cell_states_master.sh`: dmtcp, 4 CPUs, 64 GB, 8 hours. The script rebuilds
assignment by default; `plot_only=true` requires the saved cache and unchanged
MD5 fingerprints of metadata, scores, and grouping. Plotting never reruns UCell.
One two-page PDF shows all-sample proportions and a deterministic heatmap sample
of at most 500 cells per state. Plot subsampling never changes assignment or
abundance denominators. These projected states are descriptive; their
sample-centred scores do not estimate absolute cross-modality expression shifts.
