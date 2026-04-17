# AGENTS.md — snSeq_Pipeline

Single-nucleus / single-cell RNA-seq QC and analysis pipeline for the ITH all-samples dataset. Runs on Imperial College HPC with PBS Pro. The step layout follows the top-level `scRef_Pipeline` workflow, but this workspace regenerates stages 1 to 6 directly from raw inputs.

## Repository Structure

```text
*.R              — Core pipeline scripts for steps 1 to 6, executed with Rscript
*_master.sh      — PBS orchestrators for per-sample steps
N_<Step>.sh      — PBS job scripts for each pipeline step
Auto_*.sh        — Auxiliary PBS helpers for review / diagnostics
analysis/summary/    — Step-4 summary plotting scripts and outputs
analysis/plotting/   — QC heatmap and cohort-level plotting helpers
analysis/annotation/ — Annotation review plotting scripts and outputs
analysis/cell_states/ — Malignant epithelial state-mapping scripts and outputs
analysis/methodology/ — Workspace methodology write-up
sn_outs/         — All pipeline outputs for this workspace
temp/            — Optional PBS stdout/stderr staging area
AGENTS.md        — Workspace rules and execution guidance
```

## External Data Dependencies

The current pipeline is allowed to read these upstream inputs:

- `/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/pipeline.R`
- `/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/names_tmdata_sn.txt`
- Raw 10x directories referenced by `names_tmdata_sn.txt`

The current pipeline must **not** reuse precomputed annotation or downstream result files from the upstream workspace. In particular, do not depend on existing `celltypist/`, `all_samples_list_filtered.rds`, `all_filtered_outs.rds`, `all_epithelial.rds`, or related derived outputs. Stages 1 to 6 must be regenerated in this workspace from raw counts plus the rules encoded in `pipeline.R` and the first-level `scRef_Pipeline` scripts.

## Pipeline Execution Order

| Step | Shell Script | R Script | Scope | Conda Env |
|------|--------------|----------|-------|-----------|
| 1 | `1_QC_Pipeline.sh` | `QC_Pipeline.R` | All samples | `dmtcp` |
| 2 | `2_master.sh` → `2_Clustering.sh` | `Clustering.R` | Per-sample | `dmtcp` |
| 3 | `3_Annotation.sh` | `Annotation.R` | All samples | `dmtcp` |
| 4 | `4_Expr_filtering.sh` | `Expr_filtering.R` | All samples | `dmtcp` |
| 5 | `5_master.sh` → `5_InferCNA.sh` | `InferCNA.R` | Per-sample | `dmtcp` |
| 6 | `6_master.sh` → `6_Malignancy.sh` | `Malignancy.R` | Per-sample | `dmtcp` |
| 7 | `8_geneNMF.sh` | `geneNMF.R` | All malignant epithelial cells | `dmtcp` |
| 8 | `Auto_cell_states.sh` | `analysis/cell_states/*.R` | All malignant epithelial cells | `dmtcp` |

All persistent outputs for this workspace go under `sn_outs/`.

## Build / Run / Test Commands

There is no build system, linter, or unit-test suite. Execution is via PBS `qsub` or by running the wrappers directly.

```bash
# Step 1
qsub 1_QC_Pipeline.sh

# Step 2 per-sample clustering
qsub 2_master.sh

# Step 3 annotation
qsub 3_Annotation.sh

# Step 4 expression filtering and reference generation
qsub 4_Expr_filtering.sh

# Annotation review PDF
# Annotation review PDF (edit parameters directly in the R script)
qsub Auto_annotation_featureplots.sh

# Cohort QC heatmaps (post-QC/pre-expression-filter and final filtered)
qsub Auto_qc_heatmap.sh

# Step 5 per-sample inferCNA
qsub 5_master.sh

# Step 6 per-sample malignancy labelling
qsub 6_master.sh

# Step-6 integrity verification
Rscript analysis/cell_states/Auto_verify_malignancy_results.R

# Step 7 malignant epithelial merge + UCell scoring
qsub 8_geneNMF.sh

# Step 8 state mapping / visualisation
qsub Auto_cell_states.sh

# Manual per-sample submission
qsub -v sample="A_post_T1_biopsy" -N A_post_T1_biopsy 5_InferCNA.sh

# Interactive R debugging in the shared environment
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript QC_Pipeline.R
```

## HPC & File Safety Rules

These rules are inherited from the reference workspace and remain mandatory here:

1. **Working directory**: all outputs go to `sn_outs/`. Never write pipeline outputs outside project paths.
2. **Conda init**: always run `eval "$(~/miniforge3/bin/conda shell.bash hook)"` before activating environments.
3. **Interactive first**: tasks under 8 cores / 64 GB should default to `.R` script changes only when the user wants interactive execution.
4. **PBS required**: heavy tasks must have PBS `.sh` wrappers with `#PBS` resource headers.
5. **File naming**: new persistent helper files should be prefixed with `Auto_` unless they are part of the agreed pipeline step contract.
6. **Modifying existing files**: when extending existing R code, wrap substantial new blocks in 20-hash comment delimiters.
7. **No destructive cleanup**: do not delete, move, or overwrite user data outside this working directory.
8. **Ephemeral test scripts**: temporary debug scripts should use disposable names and should not become part of the pipeline contract.
9. **Max concurrent PBS jobs**: throttle per-sample submission to 46 jobs.

### Workspace-Specific Rules

1. **Raw-only regeneration**: stages 1 to 6 must be regenerated from raw counts and the rules in upstream `pipeline.R`.
2. **No precomputed annotations**: do not rely on any existing upstream `celltypist` outputs or prefiltered RDS files.
3. **Source manifest**: sample discovery must come from `names_tmdata_sn.txt`.
4. **Input reading**: raw count matrices must be loaded exactly from the directories listed in that manifest, using the same `Read10X`-based logic as the upstream script.
5. **Current sample scope**: downstream stages must use the `sn_outs/sample_manifest.csv` produced by step 1 as the authoritative sample list for the current run.
6. **Stage status gates**: per-sample downstream submission should use current status files written by the previous stage (for example `expr_filter_status.txt`, `infercna_status.txt`) rather than stale directory contents.
7. **QC thresholds**: mitochondrial, minimum genes, and housekeeping-expression thresholds must stay aligned with the upstream snRNA-seq script unless the user explicitly asks for a change.
8. **Output root**: use `sn_outs/`, not `ref_outs/`.
9. **Step 1 annotation**: keep `celltype_manual` in `QC_Pipeline.R`; do not replace it with `celltype_initial`.
10. **Step 2 annotation**: derive `celltype_initial` separately in `Clustering.R` using the 3CA marker workbook `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Marker_Genes.xlsx`.
11. **No lymph in initial annotation**: `celltype_initial` must not introduce a `lymph` label.
12. **Step 3 marker derivation**: `Annotation.R` must derive `anno_markers.rds` from confident `celltype_initial` labels and use those robust markers for `celltype_update`.
13. **Reference grouping**: inferCNA references must be built per raw source batch / technology-aligned batch (`cynthia_sn`, `multiomes`, `gemx`), not per patient.
14. **No fallback references**: if the required per-batch reference is missing or insufficient, fail that sample explicitly; do not substitute a global or patient-level reference.
15. **Cell-level filtering**: expression filtering and malignancy decisions are cell-level operations; sample-level status files are orchestration markers only.
16. **Expression-filter rescue**: step 4 currently allows a resolved lineage to pass marker support with either the standard `>=2` detected markers at the normal threshold or a single very strong assigned-lineage marker.
17. **Summary integrity**: `analysis/summary/summary.R` must stop if `snSeq_merged.rds` contains legacy cell types or if final cell-type counts do not sum to the merged cell count.
18. **Annotation review PDF**: `analysis/annotation/annotation_featureplots.R` must reflect the current implemented scoring logic from `Clustering.R`, `Annotation.R`, and the current expression-filter marker panel; do not replace it with ad hoc marker averages.
19. **Annotation sample selection**: the annotation plotting wrapper defaults to a random subset of 10 samples; use PBS variables to switch to `all` or a different subset size/seed.
20. **Annotation color stability**: lineage colors in annotation plots must be fixed and reused consistently across `celltype_initial`, `celltype_update`, and lineage score panels.
21. **Downstream malignant merge scope**: `geneNMF.R` must only use samples whose `expr_filter_status`, `infercna_status`, and `malignancy_status` are all `ok`, and it must merge only `malignant_level_1` and `malignant_level_2` epithelial cells.
22. **State-analysis batch field**: downstream state scripts must use `reference_batch` as the study replacement for normalization and plotting order; do not branch on `technology`.
23. **noreg only for snRNA-seq**: snRNA-seq state mapping uses the `noreg` workflow only. Ignore the `reg` branch from the scRef workspace.
24. **Unresolved relabel constraint**: unresolved relabelling must use the fixed scRef-retained 3CA metaprograms and must not re-derive a new retained 3CA set from the snRNA-seq cohort.

### PBS Job Template

```bash
#!/bin/bash
#PBS -l select=1:ncpus=<N>:mem=<M>gb
#PBS -l walltime=<HH:MM:SS>
#PBS -N <jobname>
echo "$(date +%T)"
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline
cd "${WD}"
Rscript <script>.R
echo "$(date +%T)"
```

## Code Style Guidelines

### R Scripts

- Put `library()` calls at the top of each file.
- Use `setwd()` near the top of each script so output paths stay relative to `sn_outs/`.
- Per-sample scripts should receive the sample name through:

```r
args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
```

- Use `snake_case` for variables and functions.
- Keep Seurat objects in `tmdata`, `tmdata_list`, `merged_obj`, `epi`, or similarly descriptive names.
- Prefer list-of-Seurat-object patterns keyed by sample name.
- Use `.rds` for data outputs, `.csv` for summary tables, and `.png` / `.pdf` for plots.

### File Organization Pattern

Per-sample outputs live under:

```text
sn_outs/by_samples/<sample>/
```

Expected step outputs include:

- `<sample>.rds` — post-QC sample object
- `<sample>_anno.rds` — annotated sample object with per-cell filtering metadata
- `anno_markers.rds` — robust OAC-specific annotation markers derived in step 3
- `annotation_score_markers.csv` — weighted refined annotation markers
- `expression_filter_markers.csv` — current expression-filter marker panel with canonical backfill
- `<batch>_reference.rds` — per-batch inferCNA reference built in step 4
- `<sample>_epi.rds` — epithelial subset with inferCNA metrics
- `<sample>_epi_f.rds` — final malignancy-labelled epithelial object
- `<sample>_outs.rds` — inferCNA output
- `<sample>_signatures.rds` — optional per-sample tumour signature genes
- `no_ref`, `no_epi`, `no_cell` — sentinel markers for skipped samples
- `analysis/summary/*.png|*.pdf` — step-4 acceptance plots
- `analysis/summary/pair3_umap_celltype_technology.png` — publication-style summary UMAP pair showing final cell type and technology
- `sn_outs/Auto_QC_snSeq_prefilter.png` — scRef-style cohort QC heatmap after QC but before expression filtering
- `sn_outs/Auto_QC_snSeq_final.png` — scRef-style cohort QC heatmap after final singlet / marker filtering
- `sn_outs/Auto_QC_snSeq_heatmaps.pdf` — two-page combined QC heatmap report
- `sn_outs/Auto_qc_heatmap_stage_summary.csv` — plotted-cell counts per stage and cell type for the QC heatmaps
- `sn_outs/snseq_annotation_featureplots_<selection>.pdf` — annotation review PDF for all or selected samples
- `sn_outs/selected_samples_<selection>.csv` — sample list used for that render
- `snSeq_malignant_epi.rds` — merged malignant epithelial object across completed step-6 samples
- `snSeq_malignant_epi_meta.rds` — metadata-only copy of the merged malignant epithelial object
- `snSeq_malignant_epi_sample_summary.csv` — per-sample malignant epithelial contribution table
- `snSeq_malignant_epi_status_input.csv` — step-6 status table used for the malignant merge
- `snSeq_malignant_epi_overall_summary.csv` — overall malignant merge summary
- `snSeq_malignant_epi_cc_top50.rds` — top malignant-epithelial consensus cell-cycle genes used for downstream CC scoring
- `snSeq_malignant_epi_cc_score.rds` — per-cell malignant-epithelial cell-cycle score used by downstream state scripts
- `Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds` — local copy of the scRef metaprogram definitions used downstream
- `Metaprogrammes_Results/UCell_nMP19_filtered.rds` — UCell scores for retained scRef metaprograms on the merged malignant epithelial object
- `UCell_3CA_MPs.rds` — UCell scores for 3CA metaprograms on the merged malignant epithelial object
- `Auto_final_states.rds` — final malignant epithelial state labels after unresolved relabelling
- `task4_unresolved_states/*` — unresolved relabelling tables and plots
- `task3_sample_abundance/*` — per-sample abundance plots and summaries
- `task6_hybrid_pairwise/*` — pairwise hybrid node-plot outputs
- `Auto_malignancy_integrity_by_sample.csv` — per-sample malignancy integrity audit
- `Auto_malignancy_integrity_summary.csv` — summary malignancy integrity audit

### Shell Scripts

- Use `#!/bin/bash` on the first line.
- Keep `echo "$(date +%T)"` at start and end.
- Load `tools/dev`, initialize conda, and activate `dmtcp`.
- Master scripts should read `sn_outs/sample_manifest.csv` for the current sample list.

### Error Handling

- Use sentinel files to mark skipped samples for downstream stages.
- Prefer status files over stale sentinels when deciding downstream submission eligibility.
- Guard empty subsets with `stop()` after writing the appropriate marker.
- Treat missing reference cells or missing epithelial cells as explicit skip states, not silent failures.

## Key Packages

Core packages in scope for this workspace:

- `Seurat`
- `Matrix`
- `dplyr`
- `tidyr`
- `ggplot2`
- `parallel`
- `infercna`
- `gridExtra`
- `UCell`
- `ComplexHeatmap`
- `circlize`
- `patchwork`

## Naming Conventions for Cell Types

Broad cell-type labels are lowercase with dots where needed, for example:

- `b.cell`
- `dendritic`
- `endothelial`
- `epithelial`
- `erythrocyte`
- `fibroblast`
- `macrophage`
- `mast`
- `neutrophil`
- `t.cell`

Ambiguous cells should be labelled `unresolved`. Non-target residual categories may be labelled `others`.
