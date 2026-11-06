# AGENTS.md — snSeq_Pipeline

Single-nucleus / single-cell RNA-seq QC and analysis pipeline for the ITH all-samples dataset. Runs on Imperial College HPC with PBS Pro. The step layout follows the top-level `scRef_Pipeline` workflow, but this workspace regenerates stages 1 to 6 directly from raw inputs.

## Repository Structure

```text
*.R              — Core pipeline scripts for steps 1 to 6, executed with Rscript
*_master.sh      — PBS orchestrators for per-sample steps
N_<Step>.sh      — PBS job scripts for each pipeline step
Auto_*.sh        — Auxiliary PBS helpers for review / diagnostics
analysis/summary/    — Step-4 summary plotting scripts; outputs go under `sn_outs/summary/`
analysis/qc/         — QC heatmap and cohort-level plotting helpers
analysis/annotation/ — Annotation review plotting scripts and outputs
analysis/infercna/   — InferCNA / malignancy audit and plotting scripts
analysis/cell_states/ — Malignant epithelial state-mapping scripts and outputs
analysis/legacy/     — Legacy comparison scripts not used downstream
analysis/lib/        — Shared analysis configuration, constants, logging, and helpers
analysis/methodology/ — Workspace and per-script methodology write-ups
analysis/ANALYSIS_MAP.md — Analysis dependency, run-order, and location map
sn_outs/         — Canonical live pipeline and analysis outputs
job_scripts/     — Local submission-only analysis launchers (gitignored)
temp/            — Historical scheduler logs; new logs go in the project root
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
| 7 | `8_geneNMF.sh` | `analysis/cell_states/prepare_malignant_epithelial_ucell.R` | All malignant epithelial cells | `gnmf` |
| 8 | `cell_states_master.sh` | `analysis/cell_states/assign_states_approach_b_noreg.R` and downstream state plotting scripts | All malignant epithelial cells | `dmtcp` |

All persistent outputs for this workspace go under `/rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/sn_outs/`. The live project root is the canonical working directory.

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
Rscript analysis/annotation/plot_annotation_featureplots.R sample_mode=subset n_samples=10 sample_seed=1

# Cohort QC heatmaps (post-QC/pre-expression-filter and final filtered)
Rscript analysis/qc/qc_heatmap.R

# Step 5 per-sample inferCNA
qsub 5_master.sh

# Step 6 per-sample malignancy labelling
qsub 6_master.sh

# Step-6 integrity verification
Rscript analysis/infercna/verify_malignancy_integrity.R

# Step 7 malignant epithelial merge + UCell scoring
qsub 8_geneNMF.sh

# Step 8 state mapping / visualisation
qsub cell_states_master.sh

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
3. **Compute nodes only**: never run analytical, memory-heavy, CPU-heavy, or substantial I/O work on a login node. Use PBS. Login-node work is limited to lightweight text/file inspections, syntax checks, and job status.
4. **PBS required**: heavy tasks must have PBS `.sh` wrappers with `#PBS` resource headers.
5. **Live Logging**: Always use live streaming log file mode by adding `#PBS -koed` to the submission script. This ensures standard out and standard error are written to their final destination as the job is running, allowing for real-time monitoring from login nodes.
6. **File naming**: outside `analysis/`, new persistent helper files should be prefixed with `Auto_` unless they are part of the agreed pipeline step contract. Inside `analysis/`, follow the Analysis Folder Rules below.
7. **Modifying existing files**: when extending existing R code, wrap substantial new blocks in 20-hash comment delimiters.
8. **File safety**: never delete files without explicit user permission. Never modify or move existing files outside this working directory. External inputs may be read or copied. Preserve unrelated user edits in dirty worktrees.
9. **Ephemeral test scripts**: temporary debug scripts should use disposable names and should not become part of the pipeline contract.
10. **Max concurrent PBS jobs**: throttle per-sample submission to 46 jobs.

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
17. **Summary integrity**: `analysis/summary/expression_filter_summary.R` must stop if `snSeq_merged.rds` contains legacy cell types or if final cell-type counts do not sum to the merged cell count.
18. **Annotation review PDF**: `analysis/annotation/plot_annotation_featureplots.R` must reflect the current implemented scoring logic from `Clustering.R`, `Annotation.R`, and the current expression-filter marker panel; do not replace it with ad hoc marker averages.
19. **Annotation sample selection**: the annotation plotting wrapper defaults to a random subset of 10 samples; use PBS variables to switch to `all` or a different subset size/seed.
20. **Annotation color stability**: lineage colors in annotation plots must be fixed and reused consistently across `celltype_initial`, `celltype_update`, and lineage score panels.
21. **Downstream malignant merge scope**: `analysis/cell_states/prepare_malignant_epithelial_ucell.R` must only use samples whose `expr_filter_status`, `infercna_status`, and `malignancy_status` are all `ok`, and it must merge only `malignant_level_1` and `malignant_level_2` epithelial cells.
22. **State-analysis batch field**: downstream state scripts must use `reference_batch` as the study replacement for normalization and plotting order; do not branch on `technology`.
23. **noreg only for snRNA-seq**: snRNA-seq state mapping uses the `noreg` workflow only. Ignore the `reg` branch from the scRef workspace.
24. **Current centred state contract**: use only the current fixed scRef centred refined MPs and grouping. Current noreg states are not relabelled by 3CA; retained-3CA relabelling is historical provenance only.

### Analysis Folder Rules

1. **Headers required**: every R analysis script must start with a 20-hash-delimited header that states status, script path, methodology path, analysis-map path, description, inputs, outputs, and usage.
2. **Methodology required**: every analysis script must have a matching methodology file under `analysis/methodology/<same-subfolder>/<script>_methodology.md`.
3. **Map required**: any new, renamed, legacy, or deleted analysis script must update `analysis/ANALYSIS_MAP.md` with dependencies, run order, terminal/legacy status, and downstream consumers.
4. **No ambiguous names**: active analysis scripts should use descriptive names without `Auto_`. Scripts retained only for comparison must use `legacy_prefix_`. Scripts proposed for manual deletion must use `delete_prefix_`.
5. **Current state method**: the active workflow is current centred refined scRef projection, Approach B, `noreg`; shared paths and thresholds live in `analysis/lib/config.R`. Canonical states are `sn_outs/cell_states/centred/intermediate/centred_refined_noreg_states.rds`.
6. **No legacy downstream outputs**: legacy scripts must write only under `sn_outs/legacy/...` and must not produce files consumed by active downstream scripts.
7. **Output tiers**: new outputs should be placed under `intermediate/`, `tables/`, `figures/`, `logs/`, or `reports/`. Root-level `sn_outs/` outputs are allowed only when they are existing canonical downstream compatibility files.
8. **Cacheable plotting**: long-running scripts must save computational intermediates under `intermediate/` and provide a documented option such as `plot_only=true` so visual styling can be regenerated without recomputing results.
9. **Logging**: long-running or reporting scripts must write a lightweight log under `logs/` containing start/end time, inputs, outputs, parameters, cache reuse, and session information where relevant.
10. **Shared helpers**: repeated constants/functions must be added to `analysis/lib/` instead of copied between scripts. Current shared files are `config.R`, `io_helpers.R`, `logging.R`, `state_helpers.R`, and `plot_helpers.R`.
11. **Presentation figures**: plots intended for slides must use readable font sizes, legends, row/column labels, and point sizes relative to figure dimensions. Avoid tiny heatmap row names or legends that cannot be read in PPTX.
12. **External dependencies**: scripts that require external files, reference databases, downloads, or upstream project outputs must list those paths in the script header and methodology.

####################
## Visualization readability for presentations and papers

- Wherever possible and space allows, use generous, clearly readable font sizes for figures intended for PowerPoint presentations and inclusion in papers. Apply this to titles, axis titles, tick labels, legends, annotations, panel labels, and heatmap row/column labels.
- Judge readability at the intended final display or printed size, including after a figure is placed on a slide or reduced to a manuscript column. A large export canvas alone does not guarantee readable text after scaling.
- Use available space to enlarge text and improve spacing without overlapping or clipped labels and legends. Increasing font size must never cause a figure to be split across additional pages. Preserve the intended comparison on one page; for clinical associations, keep all MPs on one page per clinical variable, as in OCCAMS. Adjust layout, label rotation/wrapping, margins, canvas dimensions and font sizes within that page. Only split pages when the user explicitly requests it, not to accommodate larger text. Keep point sizes and line widths clear at the intended final size.
- Choose font sizes appropriate to the output dimensions and intended use; follow journal requirements when supplied. Check the exported figure at its intended presentation or publication size before handoff.
####################

####################
## Stable current-result directories

- Each maintained analysis has one stable, descriptive output directory under the project's canonical output root. Put current files directly in `figures/`, `tables/`, `reports/`, `logs/`, and, when needed, `intermediate/` and `reference/`. Do not create empty tiers without a purpose.
- Do not create timestamped `Auto_run_*` directories, numbered attempts, `v2`/`v3` result directories, or competing final/latest copies for routine reruns or plot repairs. Timestamps and job IDs belong in logs, not in the current result path.
- Authorized reruns update the analysis's documented current outputs. Keep a single selected final figure set and one compact current-results index with working relative links. Do not overwrite unrelated files or bypass an explicit user restriction on modification.
- Current reports describe the implemented method, results, validation and unresolved scientific limitations. Exclude obsolete run inventories, failed-attempt narratives, superseded job IDs and debugging chronology from the current results index and scientific reports.
- When previous attempts must be retained, move them to an ignored project-local session archive outside the active analysis output directory, such as `Auto_session_checks/<analysis>/Auto_archive/`. Do not leave old run directories mixed with current results. Never delete files without explicit permission, and do not move files outside the project without authorization.
- Keep indispensable source data and references in the current live output tiers. Current scripts, reports and downstream analyses must not depend on the session archive or temporary audit scripts.
- When reorganizing outputs, update the producing script, registry, methodology, analysis map and links together. Verify current file paths and figure/table integrity before handoff. This rule governs new work and scoped cleanup; it does not authorize unrelated migrations.
####################

####################
## Commit-ready analysis code and figure organization

- Keep `analysis/` limited to important, maintained scientific workflows and indispensable shared helpers or setup dependencies intended for Git. Consolidate small plotting extensions into their main analysis when practical. Do not keep a separate persistent script for each repair or verification step.
- One-off preflight checks, audits, reference replays, validation probes, PDF render checks, and audit-only package installers are temporary session tools. Run them inline or place them outside `analysis/` in a Git-ignored project-local scratch directory such as `Auto_session_checks/`. They must not become required production dependencies or appear as active pipeline steps in the analysis inventory. Keep useful validation results in the output `reports/` tier.
- This rule takes precedence over older placement or naming guidance for temporary checks. A package installer is an essential tracked setup dependency only when the scientific workflow itself needs it; an audit-only renderer does not qualify. Never delete files without explicit user permission. If deletion is prohibited, retain temporary tools in the ignored session archive outside `analysis/`.
- Organize figures by purpose: a small, explicitly selected `Auto_final/` set; `Auto_supporting/` for secondary evidence; `Auto_diagnostics/` for QC; and, when needed, `Auto_archive/` for retained duplicates or obsolete exports. Use clear descriptive filenames, order final figures numerically, and provide a compact figure index explaining each file's role. Do not present all generated PDFs as necessary final figures.
- Avoid exporting both combined panels and duplicate standalone panels by default. Keep scientific source data and tests reproducible independently of temporary auditing tools. Before handoff, check that production scripts do not source temporary files and that the proposed Git changes contain no session-only tools.
####################

####################
## Analysis job scripts and PBS logs

- Store every new submission-only shell wrapper for an R script under `analysis/` in the project-root `job_scripts/` directory. Create the directory and required subdirectories when needed. Mirror the R script's path relative to `analysis/` and keep its exact basename, replacing `.R` with `.sh`: `analysis/cell_states/example.R` uses `job_scripts/cell_states/example.sh`. Preserve deeper nesting. Do not create submission-only wrappers inside `analysis/`.
- Before running or rerunning an analysis, first check that matching `job_scripts/` path. Reuse the existing launcher after checking its R command, arguments, environment, resources, and log destinations; update it there when needed. If it is absent, create it at that exact path rather than inventing another name. This basename rule takes precedence over generic shell-file naming prefixes.
- Keep these launchers minimal: PBS directives, timestamps, standard module/conda activation, a change to the project root, and the R invocation. Record the resource requirements, environment, arguments, and run command in the R script's registry so the launcher can be recreated. Keep analytical logic in R; it must run directly on an allocated compute node without reading or depending on the launcher.
- Submit from the project root, using the wrapper's `job_scripts/...` path, and change to that same root inside the job. Use `#PBS -koed` for live logging and explicit absolute `#PBS -o <project_root>/` and `#PBS -e <project_root>/` destinations. PBS stdout/stderr belong directly in the project root, alongside `analysis/` and `job_scripts/`; never place them in `analysis/`, `job_scripts/`, or an analysis output/log subdirectory. This rule concerns scheduler stdout/stderr; R-generated scientific logs and reports may still use the documented output tiers.
- Ignore `/job_scripts/` in the project's `.gitignore`. These launchers are local, reconstructable execution aids and should not be staged or committed. Keep project-root PBS stdout/stderr untracked as well.
- Essential shell scripts that implement nontrivial environment/tool installation, required input preparation, Python/helper orchestration, or other indispensable workflow logic must remain in the appropriate `analysis/` subdirectory, tracked by Git and maintained with the analysis. For example, scRef's `analysis/cell_states/Auto_drug_reversal/Auto_setup_drug_reversal_env.sh` is an essential setup script. Standard module/conda activation alone is not an exception. Do not hide indispensable logic or helper scripts solely in ignored `job_scripts/`; a disposable launcher may call a tracked essential script.
- This policy governs future launcher creation and updates. Existing pipeline entry points and historical wrappers are not automatically moved, deleted, or renamed; migrate them only as part of an explicitly scoped change.
####################

### PBS Job Template

```bash
#!/bin/bash
#PBS -l select=1:ncpus=<N>:mem=<M>gb
#PBS -l walltime=<HH:MM:SS>
#PBS -N <jobname>
#PBS -koed
#PBS -o /rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/
#PBS -e /rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/
echo "$(date +%T)"
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline
cd "${WD}"
Rscript "analysis/<subdirectory>/<script>.R"
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
- `sn_outs/summary/expression_filter/figures/*.png` and `sn_outs/summary/expression_filter/reports/*.pdf` — step-4 acceptance plots
- `sn_outs/summary/expression_filter/figures/pair3_umap_celltype_technology.png` — publication-style summary UMAP pair showing final cell type and technology
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
- `cell_states/centred/intermediate/merged_refined_mp_genes.rds` — current centred refined scRef definitions
- `cell_states/centred/intermediate/merged_refined_ucell_scores.rds` — current full-gene-ranked UCell scores
- `cell_states/centred/intermediate/centred_refined_noreg_mp_adj.rds` — current adjusted MP matrix
- `cell_states/centred/intermediate/centred_refined_noreg_group_max.rds` — current state-group scores
- `UCell_3CA_MPs.rds` — UCell scores for 3CA metaprograms on the merged malignant epithelial object
- `cell_states/centred/intermediate/centred_refined_noreg_states.rds` — canonical current centred state labels; historical `Auto_final_states.rds` is not current
- `task4_unresolved_states/*` — historical unresolved relabelling only
- `task3_sample_abundance/*` — historical uncentred per-sample outputs
- `task6_hybrid_pairwise/*` — historical uncentred hybrid outputs
- `cell_states/*/{intermediate,tables,figures,logs,reports}` — tiered cell-state outputs and logs
- `infercna/*/{intermediate,tables,figures,logs,reports}` — tiered InferCNA/malignancy review outputs and logs
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

####################
## Shared policy synchronization (2026-10)

The common governance sections are synchronized across scRef_Pipeline,
PDOs_Pipeline, Parse_Pipeline, and snSeq_Pipeline. The following rules retain
project-specific details and summarize those shared policies. They take precedence over older conflicting guidance above.

### Storage and execution

- Use `/rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline` as the
  working directory. All canonical pipeline and analysis outputs belong in its
  `sn_outs/`, including downstream inputs, per-cell scores, states, gene lists,
  normalized matrices, source data, run summaries, final figures/tables, and
  intermediates needed to replot. Ephemeral storage is allowed only for
  exceptionally large, fully reconstructable same-script caches. If another
  script reads a file, its canonical copy must be live. When uncertain, save both
  and document the canonical copy. Create output directories recursively.
- Submit analytical jobs yourself and monitor until exit status and required
  outputs are verified. A zero exit code alone is insufficient. PBS binaries may
  require `/opt/pbs/bin/qsub` and `/opt/pbs/bin/qstat`. Keep at most 46 concurrent
  jobs. Use `dmtcp` for general Seurat/analysis and `gnmf` for GeneNMF; check each
  existing workflow's documented environment before submission.
- Strictly forbid analytical fallbacks. Missing or failed required inputs must
  stop the analysis; do not silently substitute another reference or method.
- Require CPU, memory, walltime, name directives, and `#PBS -koed`. Submit from
  the live project root and change to that root inside the job. Set explicit
  absolute `#PBS -o` and `#PBS -e` destinations to that project root. Scheduler
  stdout/stderr belong directly in the root; scientific logs retain output tiers.
  `temp/` is historical log storage only.

### Launchers and maintained analysis

- Store new submission-only analysis wrappers under project-root `job_scripts/`,
  mirroring the R script's path relative to `analysis/` and exact basename with
  `.sh`. Check and reuse the matching launcher first. Keep launchers minimal;
  analytical logic belongs in tracked R scripts. Existing pipeline entry points
  are not automatically moved, renamed, or deleted.
- Keep `job_scripts/` and scheduler stdout/stderr gitignored. Indispensable
  setup/orchestration logic stays tracked under `analysis/`.
- One-off audits and validation probes run inline or in an ignored project-local
  `Auto_session_checks/` directory. They must not become production dependencies
  or active analysis-map entries. Retain temporary files when deletion is
  forbidden. Keep useful validation results in the output `reports/` tier.
- Keep `analysis/` focused on maintained workflows and shared helpers. New
  persistent helpers use `Auto_` unless an established contract specifies a
  basename. Superseded scripts use `legacy_`; manual-removal candidates use
  `delete_`. Do not rename historical files without a scoped request.
- New or substantially updated analysis registries state exact live inputs,
  outputs, downstream use, cache/replot behavior, environment, and run command.
  Complex workflows require detailed mirrored methodology. Simple deterministic
  exports/plots may explain why no separate methodology is required. Update
  `analysis/ANALYSIS_MAP.md` when workflow dependencies or outputs change.
- Organize presentation figures into a small selected `Auto_final/` set,
  `Auto_supporting/`, `Auto_diagnostics/`, and optional `Auto_archive/`. Provide a
  compact figure index. Avoid duplicate combined and standalone exports by default.
- Sample identity is `orig.ident`; never infer samples by trimming cell barcodes.

### Scientific contracts and file safety

- The current snSeq state contract projects the fixed current scRef centred
  refined 17-MP panel and grouping onto malignant epithelial cells. Use the full
  assayed-gene expression matrix for UCell ranking; never rank signature-only
  expression. Centre by orig.ident and scale by reference_batch SD, with no
  regression. Call Unresolved below 0.5 and Hybrid for a top-group gap below 0.3.
  Cell-cycle MPs are diagnostic only. No current labels receive 3CA relabelling.
  Current scRef reference files are under live ref_outs/Metaprogrammes_Results/
  centred/mp_refinement/{intermediate,tables}; current snSeq outputs are under
  sn_outs/cell_states/centred/{intermediate,tables,figures,logs,reports}.
- The updated PDO matched-pair comparison retains its historical legacy_ basename
  but uses current centred native states and adjusted scores in each dataset.
  PDO and scRef MP IDs and state spaces are independent; do not equate same-number
  MPs or silently rename native state labels. See both analysis maps and methods.
- Auto_run_standard_pipeline.sh is the tracked PBS controller for stages 2–8
  and the matched comparison after successful QC. It validates stage results,
  throttles total user jobs to 46, and stops on unexpected errors or missing
  references. Its job ledger and current status are in sn_outs/reports and logs.
- For execution requests that prohibit edits, change only explicitly authorized
  documentation. If major code errors or scientific concerns block execution,
  stop and ask before repairing code. Never modify external workspaces.
- Never delete files without explicit user permission. Never modify or move
  existing files outside this workspace. Preserve unrelated user changes.
####################


The user explicitly authorized the current centred migration and the related
PDO comparison script, methodology, map, and local launcher updates in the live
PDOs_Pipeline workspace. This is a scoped exception to the default prohibition
on external modifications; unrelated external files remain untouched.

## Full-cohort UMAP and matched native-CNA visualization

- The stage-4 `snSeq_merged.rds` embedding is an uncorrected PCA/UMAP of every
  retained filtered cell. `analysis/summary/plot_final_uncorrected_umaps.R` reuses
  it and separately builds an uncorrected embedding of every final malignant
  epithelial cell. No batch integration, regression, or cell subsampling is
  applied. Current state labels annotate, rather than define, this geometry.
  Reusable coordinates/reductions and figures are under
  `sn_outs/summary/final_uncorrected_umaps/`.
- `analysis/infercna/compare_matched_pdo_sn_cna.R` reads the explicitly matched
  PDO native InferCNA and production subclone caches from the live PDO workspace,
  validates their source provenance, preserves PDO calls, and applies the exact
  cached/current calling method independently to matched biopsy malignant cells.
  It reads the maintained PDO driver and hg38 gene order without sourcing the
  external workflow or modifying external outputs. New results stay under
  `sn_outs/infercna/matched_pdo_sn/`.
- Matched CNA heatmaps use common genomic bins and the same display colour scale,
  but preserve native reference baselines (PDO Carroll; biopsy source batch).
  Independent clone labels/colours do not imply shared ancestry. Small unassigned
  groups remain explicit in complete-cell abundance denominators. These are
  RNA-inferred relative signals, not absolute DNA copy numbers.

## Matched native-state concordance

- `analysis/cell_states/compare_matched_pdo_sn_state_concordance.R` produces
  SUR680, SUR791 and pooled reports under the stable live
  `sn_outs/cell_states/matched_pdo_sn_concordance/` tiers. The workflow adapts the
  current PDO/scRef reference concordance evidence layers without external edits.
- Each matched patient has one sample per modality. Per-patient expression
  contrasts are descriptive, without replicate-based DGE significance or cells
  treated as independent biological replicates. States with fewer than 20 cells
  are unavailable for pseudobulk effects; missing results remain NA, never zero.
- Pooled effect pairs average the same contributing matched patients equally in
  both modalities. Report one-patient versus two-patient support explicitly and
  use cell-pooled UCell means/projection percentages as in the reference script.
- Current MP signatures are scored against full assayed gene matrices. Native
  PDO and sn/scRef state spaces remain independent. Verify all 17 local sn MP
  gene vectors/grouping against the current scRef sources before scoring.
- The user specified the exact 3,538-gene effect pool from the original PDO/scRef
  concordance cache. Use only genes assayed in all matched samples (currently
  3,268); export all skipped names. Do not reselect informative genes from matched
  data or treat unavailable genes as expression zeros. Effects are descriptive.
- Each of the three reports has 10 pages: effect/marker concordance; reciprocal
  PDO-to-sn and current-scRef-to-PDO signature matrices plus gene overlap;
  match specificity; current PDO MPs projected onto sn; then one page per PDO
  state versus every sn state. No additional pooled weighting-sensitivity page.
- The projection on page 4 uses current individual PDO MPs, centres by sn sample
  and scales by reference_batch SD. Native confident states define receiving
  groups as in the original; projected Hybrid/Unresolved remain in denominators.
  All source cells remain in audit tables. The 20-cell gate applies to expression
  effects only, not UCell receiving-state means. Canonical native labels are unchanged.

## Relative matched-PDO cohort concordance

- `analysis/cell_states/rank_matched_pdo_state_concordance.R` complements the
  matched-pair analysis by ranking each biopsy state against the corresponding
  state in every eligible current PDO sample. Reuse the fixed original reference
  gene pool, restrict only by common assay availability, and keep the 20-cell gate.
- Normalize all candidate PDOs together before within-sample effects. Never mix
  old two-PDO-normalized vectors with newly normalized comparator PDOs. Biopsy
  profiles are normalized separately within their modality. Native labels stay fixed.
- Rank PDO samples, not independent patients; current cultures include repeated
  donors and treatment conditions. Report rank denominators, ties, state cell
  counts and matched-minus-best-other margins without inferential significance.
- The single selected combined PDF belongs under live
  `sn_outs/cell_states/matched_pdo_relative_concordance/figures/Auto_final/`.
  Combine both biopsies in the heatmap and rank overview; examples compare
  Matched PDO, Highest matching (highest other PDO) and Median matching (distinct
  median-rank example). Omit unmapped immune mimicry rows from figures but retain
  their unavailable source-table entries. Keep undersized biopsy states visibly
  unavailable. Heatmap row axis is snSeq state; red boxes (1.4) mark row maxima,
  black boxes (0.7) mark true matches. Omit stars and page-footer paragraphs.
  Preserve the previous analysis.

## Pan-cancer 3CA relative concordance

- `analysis/cell_states/rank_matched_pdo_3ca_concordance.R` adds an independent
  mean-MP endpoint under `sn_outs/cell_states/matched_pdo_3ca_concordance/`.
  Recompute all 67 fixed pan-cancer MPs consistently on current malignant RNA
  counts, using identical common-assayed signature genes and full assayed-gene
  UCell ranking, maxRank 1500. Preserve canonical native labels and PDO exclusions.
- Primary ranks correlate whole-sample mean UCell profiles across all 67 MPs.
  Prespecified sensitivity excludes MP1–MP5 only. Supporting native-state
  comparisons retain the existing four biological correspondences and 20-cell
  gate. Report observed ranks without tuning towards true matched identities.
- Celligner is distinct from this correlation endpoint. Its published pan-cancer
  reference alignment and 70-PC distances cannot be reproduced by relabelling a
  local PCA or MP correlation. Document feasibility and obtain the user's scope
  choice before a scientifically different adaptation or public-reference run.

## Public Celligner reference acquisition

- The user selected obtaining a pan-cancer reference for the Celligner analysis.
  `analysis/cell_states/celligner/prepare_celligner_reference.R` downloads and
  validates the original public Figshare release 11965269 v5 under live
  `sn_outs/cell_states/matched_pdo_celligner_concordance/`. Verify published MD5,
  13,485 samples, 19,188 genes, sample identities and cPC/gene alignment before use.
- The exact 2026 HCMI aligned expression and cPC files are controlled-access at
  GDC; do not mislabel the original public reference as the exact HCMI fit.
  Do not attempt restricted downloads without authorization or bypass controls.
- Public aligned expression alone is not a complete pretrained model. Keep
  reference acquisition distinct from a validated projection. Document any
  UMI/log-CPM input adaptation and never tune parameters using true-match ranks.

- `analysis/cell_states/celligner/rank_matched_pdo_celligner_concordance.R`
  specifies an exploratory projection using that public reference: common
  uniquely mapped genes, UMI log2(CPM+1), published original cPCs 1–4, reference-only
  top-1,000 differential genes and the pinned authors' unmodified MNN function
  (k1=20, k2=50, ndist=3). This reconstructed feature state is explicit and is
  not described as an unavailable complete serialized model or the HCMI fit.
  Freeze the reference-only 70-PC distance basis, export direct MNN support and
  before/after ranks, and use UMAP only as visual context. Matched identities
  enter only after alignment. Do not claim the study's FDR for local UMI pairs.
