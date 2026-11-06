# Exploratory Celligner public-reference projection

Script: `analysis/cell_states/celligner/rank_matched_pdo_celligner_concordance.R`.
Status: active terminal analysis, following verified public-reference preparation.

## Scope and distinction from the HCMI publication

The user requested obtaining a pan-cancer reference to enable Celligner-based
PDO–biopsy comparison. This workflow uses the original authors' public aligned
reference (Figshare 11965269 v5), containing 12,236 tumours and 1,249 cell lines
across 19,188 genes. It is a frozen-reference projection adapted to UMI data.
It does not reproduce the exact newer HCMI reference/integration or its FDR calls.
The HCMI aligned expression and cPC files are controlled-access at GDC; the
historical complete pretrained-model URL also returns HTTP 403.

The available aligned matrix is not a complete serialized Celligner object.
Projection state is reconstructed explicitly: use the published reference
clusters to identify reference-only differential features and the four original
public cPC loading vectors. The latter are frozen and are not recomputed from
the small local cohort. This avoids attempting to train a pan-cancer alignment
from the two matched biopsies and makes the reconstructed feature state explicit.

The exact pinned authors' `run_lm_stats_limma_group` and `modified_mnnCorrect`
functions are used unmodified. Sources:
[official code](https://github.com/broadinstitute/celligner),
[original reference study](https://doi.org/10.1038/s41467-020-20294-x), and
[HCMI study](https://doi.org/10.1038/s41586-026-10806-y).

## Inputs and eligibility

- Live local `sn_outs/snSeq_malignant_epi.rds`, `sample_manifest.csv`, and each
  represented biopsy's current expression-filter, infercna and malignancy status.
  All three statuses must be `ok`; all cells must carry malignancy level 1 or 2.
- Live read-only `PDOs_Pipeline/PDOs_outs/PDOs_merged.rds` and shared PDO config.
  Retain the canonical exclusion of SUR843T3_PDO.
- Canonical `intermediate/public_celligner_reference.rds` and pinned official
  helper source produced by `celligner/prepare_celligner_reference.R` under the
  same stable analysis output root.
- Tracked local configuration, I/O, logging and count-extraction helpers.

Every retained malignant biopsy contributes an unlabelled local profile to
reference projection, alongside every retained PDO. Require at least 20 cells
per sample, as in the current concordance gate; export excluded low-cell samples.
Samples are identified only by `orig.ident`. The true matches are SUR680 /
H_post_T1_biopsy / SUR680T3_PDO and SUR791 / L_post_T1_biopsy / SUR791T3_PDO.
Those pair identities enter scoring and visualization only after alignment.

## UMI pseudobulk and identifier mapping

Sum raw RNA counts across all retained malignant cells of each sample, including
Hybrid/Unresolved native states. Retain full-assay library totals and the raw
source vectors. Normalize to log2(1 + UMI CPM), using the full-assay library
denominator even when genes are later restricted. Do not call this bulk TPM.

Map exact approved HGNC symbols through the original public HGNC table. Use only
unambiguous one-to-one symbol/Ensembl mappings; direct Ensembl IDs can have their
version suffix removed. No guessed aliases, gene-name substitutions or zero
imputation for unassayed genes are used. Stop if multiple assay features map to
one retained reference ID. Intersect genes assayed in every included local source
with the reference, preserving reference order. Require at least 10,000 common
genes before attempting projection, and export all mapping decisions and common
IDs. This threshold is fixed before inspecting any matched-rank result.

Mean-centre each gene within the combined local cohort, and centre the retained
reference columns. This follows Celligner's input-centering logic but cannot
eliminate all UMI-versus-bulk or nucleus-versus-cell differences. Low-expression
pseudocount effects, state composition, sample quality, treatment and culture
remain important limitations.

## Frozen public-reference projection

1. Reconstruct differential features using published reference cluster labels
   and the pinned official limma helper, with trend moderation as in its R
   implementation. Rank by F statistic and retain the top 1,000 reference genes.
   Do not select features from pair identities or test alternative feature pools.
2. Restrict the four original public cPC loading vectors to the shared genes.
   Remove their span from each local profile by no-intercept least-squares
   residualization. Keep reference expression fixed. These are the original four
   public cPCs, not the HCMI paper's three newly estimated components.
3. Apply the authors' unmodified `modified_mnnCorrect` to the public reference
   and local residuals, using k1=20, k2=50 and ndist=3. The reference is large
   enough for 50 neighbours. Export every MNN pair and each sample's direct MNN
   support. Stop on missing/non-finite outputs or no global MNN support.
4. Fit a frozen 70-component, centred, unscaled PCA to the reference only. Project
   corrected local expression onto that basis. Euclidean distance in these 70
   coordinates defines the primary similarity ranking. Reference fitting does
   not use the matched sample labels or add local profiles to its PCA training.
5. Fit a reference-only UMAP of reference PCs (10 neighbours, min_dist=0.5,
   Euclidean metric, seed 1, one SGD thread), then transform local PCs. UMAP is
   visual context and never defines concordance scores.

Official MNN interpolates corrections from neighbouring anchors for target
samples without direct mutual pairs. Such samples are explicitly flagged as
extrapolated; a rank is not evidence of independent biological validation.
Audit the ten nearest reference profiles and direct MNN-anchor lineages for
every local sample. Export these neighbors and per-sample lineage context as QC.
Per-gene centring of a disease-specific UMI cohort can alter reference lineage
context; do not infer an appropriate tumour type from a high within-PDO rank
when nearest-reference context disagrees. This QC never changes scoring or
selects an alternative alignment.

## Relative ranks and diagnostic comparison

For each matched biopsy, compare its projected profile with all included PDOs.
Rank lower Euclidean distance first. Ties within 1e-10 receive the same minimum
rank; show every tied minimum with a red outline (linewidth 1.4) and the true
match with a black outline (0.7). A descriptive percentile gives half credit for
distance ties against other PDOs. It is a within-cohort rank measure, not a
probability, P value, or the paper's cross-tissue FDR.

Also project the centred local profiles before the local cPC/MNN step through
the same frozen reference PCA and show before/after relative ranks. The public
reference itself was already aligned by the original authors, so this baseline
is labelled "Before local alignment", never "raw uncorrected public data".
Do not choose between before/after methods on the basis of true-match rank.

One four-page PDF contains the all-PDO distance heatmap, matched rank contexts,
the pan-cancer reference embedding with a local view, and before/after ranks.
Keep readable labels and the comparison on its intended page; no footer
paragraphs, stars, or duplicate standalone exports. Reference GI highlighting
uses explicit upper-GI lineage annotations (esophagus or gastric/stomach), not
inferred cell names. Colorectal is retained as its own distinct lineage.

## Execution, outputs and caches

Run from the live root:

```bash
qsub job_scripts/cell_states/celligner/rank_matched_pdo_celligner_concordance.sh
qsub -v plot_only=true job_scripts/cell_states/celligner/rank_matched_pdo_celligner_concordance.sh
```

The ignored launcher uses dmtcp, 4 CPUs, 64 GB and four hours, with BLAS/OpenMP
thread limits of four, live `-koed` logging and absolute project-root stdout and
stderr. No installer or temporary audit is a production dependency.

Stable output root: `sn_outs/cell_states/matched_pdo_celligner_concordance/`.

- `intermediate/local_UMI_pseudobulks.rds`: raw/common pseudobulks, full source
  vectors, metadata, library totals, normalized profiles, gene mapping and input
  status/parameter fingerprints.
- `intermediate/celligner_projection_cache.rds`: corrected local expression,
  reference PCA/UMAP models, all local coordinates, public reference metadata,
  exact differential features, cPCs, MNN pairs and support.
- `intermediate/frozen_reference_UMAP_model.tar`: portable fitted UMAP model,
  written with `uwot::save_uwot` rather than serializing an external Annoy pointer.
- `tables/`: complete local gene mapping, shared IDs, alignment features, cells,
  excluded low-cell samples, MNN support/pairs, 70-PC coordinates, candidate
  distances, ranks and matched summaries.
- `figures/Auto_final/01_matched_PDO_Celligner_public_reference_concordance.pdf`.
- `reports/`: linked current results, page/figure indices and validation outputs.
- `logs/public_reference_projection.log`: inputs, outputs, start/end, parameters,
  reuse, session information and failures.

`plot_only=true` validates scientific input and parameter fingerprints and
reuses scores and embeddings. It does not reload RNA or rerun alignment. Public
reference preparation remains a separate indispensable input step. No external
files are modified and no files are deleted.
