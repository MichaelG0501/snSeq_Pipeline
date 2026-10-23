# Malignant epithelial merge and current centred scRef UCell scoring

## Input scope

Use only the current local `sn_outs/sample_manifest.csv`. Eligible samples require
expr_filter_status, infercna_status, and malignancy_status all equal to `ok`, plus
the final per-sample `_epi_f.rds`. Only malignant_level_1 and malignant_level_2
cells enter the merge. Preserve `orig.ident`, `reference_batch`, technology,
classification, and malignancy. No upstream annotation or derived cohort object
is reused. Stages 1–6 remain raw-only and their thresholds are unchanged.

## Current reference

Load only the live scRef centred refined gene collection at
`/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds`
and the corresponding `tables/centred_refined_mp_state_grouping.csv`. Copy genes
and grouping to live `sn_outs/cell_states/centred/`, record source paths and MD5s,
and require matching unique MP identities. Do not reapply historical silhouette
or coverage filters. Stop if any required MP has no assayed genes. Record full
signature size and detected-gene overlap for every MP.

3CA signatures come from
`/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv`.
Consensus cell-cycle genes come from
`/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv`.
These are diagnostic outputs, not additional current state-defining groups.

## Ranking and scoring

Assemble normalized RNA data for every assayed gene and every merged cell,
retaining exact cell names across Seurat layers. UCell must rank the full gene
universe. Never restrict expression rows to signature genes before ranking:
that changes ranks and makes scores dependent on the chosen signature panel.
Use the gnmf environment, four worker cores, and fixed seed 1. Save current MP
UCell scores to `cell_states/centred/intermediate/merged_refined_ucell_scores.rds`.
Score 3CA signatures on the same full expression matrix and retain
`UCell_3CA_MPs.rds` as a diagnostic compatibility output. Use the top 50 expressed
consensus cell-cycle genes to retain the existing diagnostic CC score outputs.

## Storage, execution, and validation

`qsub 8_geneNMF.sh` requests 4 CPUs, 128 GB, 8 hours, and gnmf. It always rebuilds
from current stage-6 objects and reference files; it never reuses historical
nMP19 scores. Canonical malignant Seurat/metadata objects and contribution/status
summaries remain in live `sn_outs/`; current MP inputs/scores, overlap tables, and
provenance use `cell_states/centred/` tiers. The scoring run log lives under
`cell_states/ucell_scoring/logs/`. The controller runs the existing malignancy
integrity audit before this merge. Assignment then requires complete barcode/MP
coverage and finite scores. Do not present partial-cohort outputs as complete.
