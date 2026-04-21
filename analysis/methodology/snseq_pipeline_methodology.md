# snRNA-seq Pipeline Methodology

This document describes how the current `snSeq_Pipeline` code is working for:

1. Initial filtering
2. Initial annotation
3. Final annotation
4. Expression filtering

It is written against the current code in:

- `QC_Pipeline.R`
- `Clustering.R`
- `Annotation.R`
- `Expr_filtering.R`

The description below is intentionally operational. It documents the implemented logic, not a simplified conceptual version.

## 1. Initial Filtering (`QC_Pipeline.R`)

### 1.1 Raw input discovery

The pipeline reads sample paths from:

- `source_root = /rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples`
- `names_tmdata_sn.txt`

Each entry in `names_tmdata_sn.txt` is treated as a relative path under `source_root`. The sample ID is taken as `basename(path)` and the source/study is taken as `dirname(path)`.

### 1.2 Raw data loading

For each sample:

1. `Read10X()` is used to load the raw matrix.
2. If the 10x object is a list, `"Gene Expression"` is used when present, otherwise the first element is taken.
3. A Seurat object is created from the raw counts.

The following metadata are attached immediately:

- `orig.ident = sample_id`
- `input_source = source_batch`
- `reference_batch = source_batch`
- `technology = classify_technology(source_batch)`
- `sample_id = sample_id`

Technology classification is implemented as:

- `gemx -> GEMX`
- `cynthia_sn -> snSeq`
- everything else -> `Multiome`

### 1.3 Mitochondrial filtering

Mitochondrial reads are quantified using:

- `PercentageFeatureSet(pattern = "^MT-")`

Cells are retained only if:

- `percent.mt < 15`

This is the first hard filter.

### 1.4 Library normalization

After mitochondrial filtering, each sample is normalized manually:

1. Library size per cell is computed from raw counts.
2. CPM is computed as `counts / library_size * 1e6`
3. A log-like expression layer is then defined as:

   `log2((CPM / 10) + 1)`

These matrices are stored as:

- `RNA@assays$CPM`
- `RNA@assays$data`

This `data` layer is the expression layer used by later annotation and filtering steps.

### 1.5 Gene-count and housekeeping filtering

Two QC statistics are computed from the `data` layer:

1. `n_genes`:
   the number of genes with `data > 0`
2. `hk_mean`:
   the mean expression of the housekeeping gene set

Housekeeping genes:

- `ACTB`, `GAPDH`, `RPS11`, `RPS13`, `RPS14`, `RPS15`, `RPS16`, `RPS18`, `RPS19`, `RPS20`, `RPL10`, `RPL13`, `RPL15`, `RPL18`

Cells are retained only if:

- `n_genes >= 200`
- `hk_mean >= 0.5`

The housekeeping filter is applied after the gene-count filter.

### 1.6 Legacy `celltype_manual` field

After QC, the code still computes a legacy marker-based label named `celltype_manual`.

This uses a broad canonical/contaminant marker panel containing:

- `b.cell`
- `dendritic`
- `endothelial`
- `epithelial`
- `erythrocyte`
- `fibroblast`
- `keratinocyte`
- `lymph`
- `macrophage`
- `mast`
- `neutrophil`
- `t.cell`

For each cell type, the score is:

- the mean log-expression of the available marker genes in that cell

`celltype_manual` is then assigned as the single highest-scoring lineage only if:

- the maximum score is unique
- the maximum score is at least `1`

otherwise the cell is labelled:

- `unclassified`

Important: this `celltype_manual` field is only a legacy QC-era label. It is not the confident cluster-based annotation used later for `celltype_initial` or `celltype_update`, and it is not the marker panel now used for step-4 expression filtering.

### 1.7 Step-1 outputs

Per sample:

- `sn_outs/by_samples/<sample>/<sample>.rds`

Global outputs:

- `sample_manifest.csv`
- `filtering_summary_sn.csv`
- `Inspections_sn.pdf`
- `cells_filtering_sn.pdf`

## 2. Initial Annotation (`Clustering.R`)

This step creates `celltype_initial`.

### 2.1 Per-sample clustering

Each sample is processed independently.

If `n_cells > 50`:

1. Select `3000` highly variable genes with `FindVariableFeatures()`
2. Scale the data with `ScaleData()`
3. Run PCA
4. Build a neighbour graph
5. Run Leiden clustering with:
   - `resolution = 0.8`
   - `algorithm = 4`
   - `method = "igraph"`
6. Run UMAP

The PCA helper retries with `approx = FALSE` if `RunPCA()` fails because of a null-space issue.

Neighbourhood settings are dataset-size dependent:

- If `n_cells > 20000`:
  - `k.param = 10`
  - PCs `1:min(20, pca_npcs)`
  - UMAP neighbours `15`
- Otherwise:
  - `k.param = min(30, max(10, round(sqrt(ncells) - 2)))`
  - PCs `1:min(30, pca_npcs)`
  - UMAP neighbours `30`

If `n_cells <= 50`, the sample is still clustered, but a smaller-data branch is used:

- `FindVariableFeatures(nfeatures = 3000)`
- `ScaleData(features = VariableFeatures(tmdata))`
- PCA with `npcs = max(2, min(length(HVGs), ncells - 1))`
- `k.param = max(5, min(15, ncells - 1))`
- Leiden clustering and UMAP on those PCs

However, if `n_cells <= 50`, the final label is still set to:

- `celltype_initial = "unknown"`

### 2.2 3CA pan-cancer marker input

The code reads:

- `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Marker_Genes.xlsx`

Marker rows are restricted to:

- `specificity > 0.2`
- `cell_type != "Malignant"`

Cell type labels are remapped to the pipeline naming convention:

- `Fibroblast -> fibroblast`
- `Macrophage -> macrophage`
- `Mast -> mast`
- `B cell -> b.cell`
- `T cell -> t.cell`
- `Dendritic -> dendritic`
- `Endothelial -> endothelial`
- `Epithelial -> epithelial`
- `NK cell -> nk.cell`
- `Plasma -> plasma`

### 2.3 Marker ranking

Within each cell type, 3CA genes are ranked by a weighted combined score:

- specificity rank weight = `0.2`
- sensitivity rank weight = `0.8`

Ranks are converted to percentile-like values with:

- `rank / (number_of_non_missing + 1)`

The combined score is:

`(0.2 * percentile_rank(specificity) + 0.8 * percentile_rank(sensitivity)) / 1.0`

The top `50` genes per cell type are retained for initial annotation.

### 2.4 Initial annotation method 1: cluster DEG enrichment

Cluster DEGs are computed with:

- `FindAllMarkers(only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)`

These DEGs are then further restricted to:

- `pct.1 >= 0.5`

For each cluster and each cell type:

1. Overlap between cluster DEGs and the top-50 3CA list is counted
2. A one-sided Fisher exact test is performed
3. P-values are adjusted with BH correction

Only cell types with:

- `padj <= 0.05`

are considered significant.

Calling rule:

- If only one significant cell type remains, assign it
- If multiple significant cell types remain, the top one is called only if:
  - `top_overlap - second_overlap >= 3`
- Otherwise, all cell types within `3` overlap genes of the top hit and with `overlap > 3` are retained and concatenated with `|`
- If nothing qualifies, the cluster is `unknown`

This result is stored as `step1`.

### 2.5 Initial annotation method 2: 3CA module-score signatures

For the same top-50 marker sets:

1. `AddModuleScore()` is run per cell type
2. Scores are aggregated to cluster medians
3. Scores are Z-normalized across cell types within each cluster

The helper `gap_signature_call()` then assigns the cluster using:

- `z_cut = 1`
- `z_gap = 0.5`

Calling rule:

- the top cell type must have `top_z >= 1`
- the difference between the top and second score must be `>= 0.5`
- otherwise the cluster is `unknown`

This result is stored as `step2`.

### 2.6 Final `celltype_initial` call

`step1` and `step2` are combined cluster-by-cluster:

- if both are single and identical: keep that label
- if one is single and the other is missing/unknown: keep the single label
- otherwise: store the disagreement explicitly as `step1 <> step2`

Cells inherit the cluster-level result as:

- `celltype_initial`

If a cluster was not mapped, it is filled as:

- `unknown`

### 2.7 Step-2 outputs

Per sample:

- updated `by_samples/<sample>/<sample>.rds`
- `by_samples/<sample>/<sample>_cluster_overview.png`

The step-2 overview currently shows:

- UMAP colored by `seurat_clusters`
- UMAP colored by `celltype_initial`

## 3. Final Annotation and OAC Marker Derivation (`Annotation.R`)

This step creates:

- `anno_markers.rds`
- `annotation_score_markers.csv`
- `celltype_update`

### 3.1 Input state

Each QC-passing sample is read from:

- `by_samples/<sample>/<sample>.rds`

If a sample has fewer than `50` cells, it is not used for marker derivation and is initialized as:

- `celltype_update = "unresolved"`

### 3.2 3CA marker preparation

The same 3CA marker table is reloaded.

At this step:

- top `100` genes per cell type are retained after combined ranking

These 100-gene 3CA lists are used only as the candidate universe for deriving the OAC/snRNA marker set.

### 3.3 High-confidence seed-cell construction

The confident seed is not taken from all exact `celltype_initial` cells.

Instead, each sample builds a stricter per-cell seed using the canonical marker panels:

- fibroblast
- macrophage
- mast
- epithelial
- t.cell
- b.cell
- nk.cell
- plasma
- dendritic
- endothelial

For each cell:

1. Canonical marker expression is binarized at:
   - `data > 1`
2. For each lineage, the proportion of canonical markers expressed is computed
3. The top lineage and its gap to the second lineage are recorded

A cell is retained as a confident seed only if:

- `celltype_initial` is a valid lineage
- the canonical top lineage is unique
- the canonical top lineage matches `celltype_initial`
- the top-minus-second gap is at least `0.1`

The resulting per-cell label is stored as:

- `celltype_initial_seed`

Diagnostic outputs:

- `celltype_initial_seed_summary.csv`
- `celltype_initial_seed_retention.csv`

### 3.4 Sample-wise marker exclusivity calculation

For each candidate 3CA marker gene and each seeded lineage:

1. Expression is binarized at:
   - `data > 1`
2. For every sample and lineage, the percentage of cells expressing that gene is computed

This produces a sample-by-lineage gene-expression frequency table.

### 3.5 Allowed lineage exclusions in off-target checks

Off-target expression is not evaluated naively across all other lineages.

The code explicitly allows some related pairs when building exclusivity tables:

- `b.cell` ignores `plasma` and `dendritic`
- `plasma` ignores `b.cell` and `dendritic`
- `t.cell` ignores `nk.cell` and `dendritic`
- `nk.cell` ignores `t.cell` and `dendritic`
- `macrophage` ignores `dendritic`

These exclusions are only for marker derivation, not for the final cell labels themselves.

### 3.6 Exclusive-marker derivation

For each gene and home lineage, two quantities are computed:

1. `home_pct`
   the percentage of samples in which the home lineage passes the home-side criterion
2. `off_pct_max`
   the maximum percentage of samples in which any allowed off-target lineage passes the off-side criterion

Off-target expression is always tested with:

- off-side expression threshold: `> 15% of cells`
- off-target sample cap: `< 30% of samples`

The home-side criterion is relaxed over a fixed grid:

1. `strict`:
   - home expression threshold `> 30% of cells`
   - must pass in `> 15% of samples`
2. `home_expr25`:
   - home expression threshold `> 25%`
   - sample requirement `> 15%`
3. `home_expr20_samp10`:
   - home expression threshold `> 20%`
   - sample requirement `> 10%`
4. `home_expr15_samp5`:
   - home expression threshold `> 15%`
   - sample requirement `> 5%`
5. `home_expr10_samp1`:
   - home expression threshold `> 10%`
   - sample requirement `> 1%`

A marker is considered exclusive within a given tier only if:

- `home_pct > home_min_pct_samps`
- `off_pct_max < 30`

All tier tables are saved as:

- `exclusive_tbls.rds`
- `exclusive_tbls_long.rds`
- `exclusive_tbls.csv`

### 3.7 Building the final OAC/snRNA marker list

The final marker list begins from the `strict` exclusive tier only.

Then, for every required lineage except dendritic:

- if fewer than `4` markers are present, the code adds more exclusive genes from progressively relaxed tiers
- genes are added in order of:
  - higher `home_pct`
  - lower `off_pct_max`
  - higher combined 3CA ranking

Required lineages are:

- all cell types present in the 3CA list except `dendritic`

If any required lineage is still missing after this relaxation, the script stops with an error.

The final marker object is saved as:

- `anno_markers.rds`

### 3.8 Scoring weights for the refined marker list

Each final marker gene is assigned a weight:

- `score_weight = max(home_pct - off_pct_max, 0.1)`

This weighted marker table is saved as:

- `annotation_score_markers.csv`

For diagnostics, the top four markers per lineage by `home_pct`/`off_pct_max` order are also written to:

- `annotation_score_markers_top4.csv`

### 3.9 Refined cluster scoring

For each sample:

1. Only refined marker genes present in the sample are used
2. Expression is taken from the `data` layer
3. Expression is Z-normalized gene-by-gene across cells
4. For each lineage, a weighted linear score is computed per cell from the refined markers

These per-cell weighted scores are then summarized to cluster level in two separate ways:

- median per cluster
- mean per cluster

For each aggregation, cluster scores are Z-normalized across lineages within the cluster.

### 3.10 Refined cluster calling

Cluster calls are made from these Z-scores using:

- `gap_cut = 0.8`

Allowed close pairs are:

- `t.cell | nk.cell`
- `nk.cell | t.cell`
- `b.cell | plasma`
- `plasma | b.cell`

Calling rule:

- if the top lineage is at least `0.8` above all others, call the top lineage
- if exactly one competing lineage lies within `0.8` and it forms an allowed pair with the top lineage, call the pair
- otherwise the cluster is `unresolved`

This is applied separately to:

- cluster medians
- cluster means

### 3.11 Canonical-marker rescue

Clusters that still fail refined scoring can be rescued using canonical markers.

For each cluster and lineage:

1. canonical marker expression is binarized at `data > 1`
2. the per-gene fraction of expressing cells is computed
3. lineage score = mean of those fractions
4. `ngenes_gt10` = number of canonical markers expressed in more than `10%` of cells

Canonical rescue is allowed only for clusters with:

- at least `50` cells

Calling thresholds:

- Immune lineages (`b.cell`, `t.cell`, `nk.cell`, `plasma`)
  - top score `>= 0.20`
  - second score `<= 0.05`
  - `ngenes_gt10 >= 3`
- Other lineages
  - top score `>= 0.12`
  - second score `<= 0.08`
  - `ngenes_gt10 >= 2`

If these conditions are not met, the canonical rescue call is `unresolved`.

### 3.12 Choosing the final annotation

The code keeps several levels of annotation:

- `celltype_refined_median`
- `celltype_refined_mean`
- `celltype_canonical`
- `celltype_refined`
- `celltype_update`
- `celltype_update_detail`

The decision logic is:

1. `choose_refined_call()`
   combines median and mean refined calls
2. `choose_final_call()`
   combines the refined call with canonical rescue and the original `celltype_initial`

Important implemented behavior:

- if median and mean agree on a resolved label, that label is used
- if only median is resolved, median is preferred
- mean-only calls are only trusted more selectively
- canonical rescue is only allowed to replace unresolved refined calls in sufficiently large clusters
- if refined/canonical still fail, a single confident `celltype_initial` label is preserved
- otherwise the final label is `unresolved`

`celltype_update_detail` keeps disagreement information via:

- `label1 <> label2`

when initial and updated evidence are inconsistent.

### 3.13 Step-3 outputs

Per sample:

- `*_anno.rds`
- `*_cluster_scores.csv`
- `*_cluster_scores_median.csv`
- `*_cluster_scores_mean.csv`
- `*_cluster_canonical_scores.csv`
- `*_cluster_annotation_summary.csv`
- `*_cluster_overview.png`

The current cluster overview now shows three UMAPs:

- clusters
- `celltype_initial`
- `celltype_update`

Global outputs:

- `anno_markers.rds`
- `annotation_score_markers.csv`
- `exclusive_tbls.*`
- `pct_mat.csv`
- `count_mat.csv`

## 4. Expression Filtering (`Expr_filtering.R`)

This step now uses the derived OAC/snRNA marker panel rather than the old legacy filtering panel.

### 4.1 Marker panel used for filtering

The filtering panel is built from:

- `annotation_score_markers.csv`

This file contains the final refined marker genes from step 3, together with:

- `home_pct`
- `off_pct_max`
- `score_weight`

The filter panel is built lineage-by-lineage as follows:

1. Take all refined markers with:
   - `score_weight >= 5`
2. If none exist for a lineage, fall back to the top derived markers
3. Add a curated canonical backfill order to improve sensitivity
4. Truncate to a lineage-specific target panel size

Target panel sizes:

- epithelial: `7`
- endothelial: `8`
- fibroblast: `8`
- macrophage: `8`
- b.cell: `6`
- plasma: `6`
- t.cell: `6`
- nk.cell: `6`
- mast: `4`
- dendritic: `3`

This guarantees at least three markers per lineage.

Important change:

The old legacy categories:

- `lymph`
- `erythrocyte`
- `keratinocyte`
- `neutrophil`

are no longer part of the active expression-filter or doublet-filter marker panel.

The active filter lineages are now aligned with the annotated cell types:

- epithelial
- endothelial
- fibroblast
- macrophage
- b.cell
- plasma
- t.cell
- nk.cell
- mast
- dendritic

The final panel is written to:

- `expression_filter_markers.csv`
- `expression_filter_marker_counts.csv`

### 4.2 Per-cell marker scoring

For each sample:

1. Expression is taken from the `data` layer
2. Only genes present in the sample are used
3. For each lineage:
   - `detected_markers = number of panel genes with data > 0.5`
   - lineage score = mean of the top up to `4` expressed markers above `0.5`

This scoring is intentionally more sensitive for sparse snRNA-seq data than requiring every marker in the panel to be expressed.

### 4.3 Doublet screening

Two doublet flags are computed.

#### Strict coexpression flag

A lineage is considered strongly present in a cell if:

- lineage score `> 1.0`
- at least `2` lineage markers are detected

The cell is marked as:

- `coexpression = "doublet"`

if more than one lineage is strongly present, except for allowed non-exclusive lineage groups:

- macrophage / mast / dendritic
- t.cell / nk.cell / dendritic
- b.cell / plasma

Otherwise:

- `coexpression = "singlet"`

#### Loose coexpression flag

The same high-scoring lineages are used, but the code also checks whether the highest lineage is sufficiently separated from the next competing lineage.

The cell is marked as:

- `coexpression_loose = "doublet"`

if:

- more than one lineage is strongly present
- and the highest lineage is less than `1` score unit above the next non-exempt competitor

This is the stricter flag used for the final keep/remove decision.

### 4.4 Expression-quality classification

Each cell is then classified using its current `celltype_update`.

For resolved cells:

- the assigned lineage is accepted as `good` if either:
  - score `>= 1.25` and detected markers `>= 2`
  - or score `>= 5` and detected markers `>= 1`

If so:

- `marker_expression = "good"`

If the assigned lineage fails but another lineage passes the same gate:

- `marker_expression = "good_inconsistent"`

If no lineage passes:

- `marker_expression = "poor"`

For unresolved cells:

- if any lineage passes the same gate:
  - `marker_expression = "good_unresolved"`
- otherwise:
  - `marker_expression = "poor_unresolved"`

The code also stores:

- `marker_match_score`
- `marker_match_detected`
- `max_filter_score`
- `max_filter_celltype`
- `max_filter_detected`

### 4.5 Suspicious-cluster reclustering

The script then revisits suspicious clusters among singlet cells that are not already clearly poor.

A cluster is selected for reclustering if either:

1. at least `10%` of eligible cells are `good_inconsistent` and the cluster contains at least `100` eligible cells
2. the cluster contains more than `50` `good_unresolved` cells

For each such cluster:

1. subset singlet cells from that original cluster
2. re-run:
   - `FindVariableFeatures(nfeatures = 3000)`
   - `ScaleData()`
   - PCA
   - neighbour graph
   - Leiden clustering (`resolution = 0.8`, `algorithm = 4`)
   - UMAP
3. re-score the new subclusters with the refined expression-filter marker panel
4. convert cluster medians to within-cluster Z-scores
5. re-call the subclusters using the same `gap_cut = 0.8` and the same allowed lineage pairs as step 3

Cells in these reclustered suspicious regions can therefore receive an updated `celltype_update` before the final keep/remove decision.

The flag:

- `expr_filter_reclustered`

records whether a cell went through this rescue step.

### 4.6 Final keep/remove rule

The final filter decision is still singlet-only, but it is no longer as harsh on sparse nuclei with one highly convincing marker.

A cell is retained only if:

- `coexpression_loose == "singlet"`
- and `marker_expression == "good"`

Everything else is removed.

The final reason is stored in:

- `filter_reason`

with categories:

- `keep`
- `doublet`
- `marker_inconsistent`
- `marker_positive_unresolved`
- `poor_assigned`
- `poor_unresolved`
- `too_few_cells`
- `other`

The boolean keep flag is:

- `filter_keep`

### 4.7 Filtered outputs

Per sample:

- updated `*_anno.rds` with filtering metadata
- `*_filtered.rds` containing only retained cells
- `expr_filter_status.txt`
- `no_cell` when nothing survives

Global outputs:

- `filtered_sample_summary.csv`
- `expression_filter_status_by_sample.csv`
- `expression_filter_reason_by_sample.csv`
- `expression_filter_reason_overall.csv`
- `filtered_celltype_counts.csv`
- `snSeq_merged.rds`
- `snSeq_filtered.png`

### 4.8 Reference construction after filtering

After filtering, batch-specific reference objects are constructed from the retained cells.

For each `reference_batch`:

1. all filtered sample CPM matrices in that batch are merged
2. metadata are merged
3. only cells satisfying both are kept as reference candidates:
   - `coexpression == "singlet"`
   - `celltype_initial == celltype_update`
4. reference lineages are chosen from this ordered priority:
   - macrophage
   - endothelial
   - fibroblast
   - t.cell
5. the first two available lineages are used as the reference set

If fewer than two of these lineages are available, that batch is skipped for reference generation.

## 5. Summary plotting

Step 4 now also runs:

- `analysis/summary/summary.R`

This produces review plots for step-4 acceptance, including:

- total cells after each filter
- final distribution by technology
- expression-filter outcome categories / why cells were filtered
- final cell-type composition
- side-by-side UMAPs for final cell type and technology
- per-sample retention

These plots are saved under:

- `analysis/summary/`
- including `analysis/summary/pair3_umap_celltype_technology.png`

## 6. Cohort QC Heatmaps

The workspace also provides a scRef-style QC heatmap helper generated from:

- `analysis/plotting/qc_heatmap.R`
- `Auto_qc_heatmap.sh`

This helper ports the original `scRef_Pipeline/analysis/plotting/qc_heatmap.R` layout and legends to the snRNA-seq workspace, while keeping the snRNA-seq QC thresholds:

- mitochondrial DNA `< 15`
- detected genes `> 200`
- housekeeping expression `> 0.5`

Because the snRNA-seq cohort is much larger than the original per-sample scRef example, the script builds lightweight marker-only cohort objects and downsamples up to `20,000` cells per stage for plotting, stratified by `celltype_update`, so the figure remains runnable while preserving the original visual style.

Two stage-level heatmaps are written:

1. `sn_outs/Auto_QC_snSeq_prefilter.png`
   - all post-QC cells before expression filtering
   - includes unresolved / low-quality states present at this stage
2. `sn_outs/Auto_QC_snSeq_final.png`
   - final filtered singlets with retained marker support
3. `sn_outs/Auto_QC_snSeq_heatmaps.pdf`
   - two-page combined report
4. `sn_outs/Auto_qc_heatmap_stage_summary.csv`
   - plotted-cell counts per stage and per `celltype_update`

## 7. Annotation Visualization PDF

The workspace also provides a dedicated annotation review PDF generated from:

- `analysis/annotation/annotation_featureplots.R`
- `Auto_annotation_featureplots.sh`

This plot set is built directly from the current annotated per-sample objects `sn_outs/by_samples/<sample>/<sample>_anno.rds`.

### 7.1 Sample selection

The plotting script can run on:

- all samples
- or a random subset of `n` samples

The PBS wrapper currently defaults to:

- `sample_mode = subset`
- `n_samples = 10`
- `sample_seed = 1`

The sample selection logic is controlled by variables at the top of the R script:

- `sample_mode`: Choices are `"subset"`, `"all"`, or a specific sample ID (e.g., `"patient_A_pre"`).
- `n_samples`: Number of samples to pick if mode is `"subset"`.
- `sample_seed`: Seed for subset reproducibility.

To change which samples are plotted, edit these values in `analysis/annotation/annotation_featureplots.R` and then run:

```bash
qsub Auto_annotation_featureplots.sh
```

For each sample, four pages are produced, followed by a blank separator page before the next sample:

1. **Initial 3CA Top-50 Scores**
   - one UMAP panel coloured by `seurat_clusters`
   - ten feature plots showing per-cell top-50 3CA lineage Z-scores for all ten lineages
   - these scores are built as:
     - gene-wise Z-normalized expression across cells
     - unweighted mean across the top 50 lineage markers
     - per-cell Z-normalization across the ten lineages
   - one UMAP panel coloured by `celltype_initial`
   - diverging blue–white–red colour scale, fixed at `[-2.5, 2.5]`
   - **interpretation note**: the top-50 genes are selected with 0.8 sensitivity / 0.2 specificity weighting, so many genes are broadly expressed. This is intentional — `Clustering.R` uses these for cluster-level methods (DEG enrichment and AddModuleScore with cluster-median aggregation) where broad sensitivity ensures robust overlap. At the cell level, individual cells may appear similar across lineages; this confirms why the pipeline correctly uses cluster-level aggregation rather than per-cell thresholds.

2. **Refined Weighted Scores**
   - one UMAP panel coloured by `seurat_clusters`
   - ten feature plots showing per-cell weighted refined Z-scores using the current expression-filter marker panel
   - this panel uses:
     - derived OAC/snRNA markers from `annotation_score_markers.csv`
     - canonical backfill genes from `expression_filter_markers.csv`
     - per-gene weights from `score_weight`
     - canonical backfill genes without an empirical `score_weight` are included with the floor weight `0.1`
   - these scores are then per-cell Z-normalized across the ten lineages
   - one UMAP panel coloured by `celltype_update`
   - diverging blue–white–red colour scale, fixed at `[-2.5, 2.5]`

3. **Canonical Marker Expression**
   - one UMAP panel coloured by `seurat_clusters`
   - ten feature plots showing per-cell mean log-expression of canonical markers for each lineage
   - one UMAP panel coloured by `celltype_update`
   - sequential grey-to-red colour scale, auto-scaled to the data maximum (shared across all ten lineages)

4. **Expression-Filter Panel Scores**
   - one UMAP panel coloured by `seurat_clusters`
   - ten feature plots showing per-cell expression scores using the step-4 filter panel
   - score is the mean of the top up-to-4 expressed markers above `0.5`, matching the `Expr_filtering.R` scoring logic
   - one UMAP panel coloured by **filter status** (Kept / Low Expr / Others) instead of `celltype_update`
   - sequential grey-to-red colour scale, auto-scaled to the data maximum (shared across all ten lineages)

### 7.2 Visualization conventions

The annotation PDF enforces:

- a fixed lineage order across all pages
- **diverging blue–white–red** colour scale for Z-score pages (pages 1 and 2): negative values appear blue, zero is white, positive values are red — this clearly separates below-average from above-average lineage scores
- **sequential grey-to-red** colour scale for expression score pages (pages 3 and 4): grey is zero, red is the shared maximum
- a colour bar legend on each score panel
- shared scale across all ten lineage panels within each page for direct comparability
- increased adaptive point sizes (~2× larger than original) for better colour visibility
- a visible grey/white baseline so low-score cells remain visible
- a separate legend on the `celltype_initial` / `celltype_update` / filter status panel

The output file is:

- `sn_outs/snseq_annotation_featureplots_<selection>.pdf`
- `sn_outs/annotation_featureplots_status_<selection>.csv`
- `sn_outs/selected_samples_<selection>.csv`

## 8. Practical interpretation

### `celltype_manual`

`celltype_manual` is still generated in step 1, but it is a legacy broad marker-based QC label. It is not the cluster-based confident annotation.

### `celltype_initial`

`celltype_initial` is the first structured, cluster-level annotation produced from two complementary 3CA-based methods in step 2.

### `celltype_update`

`celltype_update` is the final step-3 cluster-level annotation after:

- seed-cell refinement
- OAC/snRNA marker derivation
- weighted refined-marker scoring
- canonical rescue
- preservation of confident initial singlets when appropriate

### Step-4 filtering

Step 4 does not simply check canonical markers. It now uses the OAC/snRNA-derived marker list, with controlled canonical backfill only where needed for sensitivity, and then filters at the cell level for:

- marker support
- inconsistency
- unresolved marker-positive cells
- doublet-like coexpression

## 9. Malignancy Integrity Verification

The workspace now includes:

- `analysis/cell_states/Auto_verify_malignancy_results.R`

This is a post-step-6 audit. It does not modify any malignancy objects. It verifies that the current malignancy outputs are internally consistent before any merged malignant-epithelial downstream analysis is run.

### 9.1 Samples included in the audit

Samples are included only if all of the following are present and marked complete:

- `expr_filter_status.txt = ok`
- `infercna_status.txt = ok`
- `malignancy_status.txt = ok`
- `<sample>_epi_f.rds`
- `<sample>_malignancy_summary.csv`
- `<sample>_malignancy_cluster_summary.csv`

The audit uses `sn_outs/sample_manifest.csv` as the authoritative sample list.

### 9.2 Checks performed

For each completed sample, the verifier loads `<sample>_epi_f.rds` and checks:

- required metadata columns are present:
  - `orig.ident`
  - `reference_batch`
  - `classification`
  - `malignant_clus`
  - `malignancy`
- epithelial cell count matches `<sample>_malignancy_summary.csv`
- total `classification` counts equal the object cell count
- total `classification` counts equal the summary-table classification totals
- total `malignancy` counts equal the object cell count
- total `malignancy` counts equal the summary-table malignancy totals
- total cells across `<sample>_malignancy_cluster_summary.csv` equal the object cell count
- each object contains a single `orig.ident`
- each object contains a single `reference_batch`

The script stops if any completed sample fails one or more checks.

### 9.3 Outputs

The verification step writes:

- `sn_outs/Auto_malignancy_integrity_by_sample.csv`
- `sn_outs/Auto_malignancy_integrity_summary.csv`

## 10. Malignant Epithelial Merge and UCell Scoring (`geneNMF.R`)

The workspace now includes:

- `geneNMF.R`
- `8_geneNMF.sh`

This stage creates the merged malignant epithelial snRNA-seq object and scores both the scRef metaprograms and the external 3CA metaprograms on that merged object.

### 10.1 Input scope

Only samples with:

- `expr_filter_status = ok`
- `infercna_status = ok`
- `malignancy_status = ok`

are eligible.

For each eligible sample, only epithelial cells with:

- `malignancy = malignant_level_1`
- or `malignancy = malignant_level_2`

are retained for the merged malignant-epithelial object.

Cells labelled `unresolved` or non-malignant are not included in the merged object.

### 10.2 scRef metaprogram definitions

The metaprogram definitions are copied from:

- `/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds`

They are then stored locally in:

- `sn_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds`

Only metaprograms with:

- `silhouette >= 0`

are retained for downstream UCell scoring.

### 10.3 Merge construction

The script loads each eligible `<sample>_epi_f.rds`, subsets malignant cells, and merges them into:

- `sn_outs/snSeq_malignant_epi.rds`

If cell names are duplicated across samples, cells are renamed with the sample prefix before merging.

The metadata-only companion file is:

- `sn_outs/snSeq_malignant_epi_meta.rds`

### 10.4 UCell scoring of scRef metaprograms

The merged malignant epithelial object is a Seurat v5 object with split per-sample `data.*` layers after merge. To avoid ambiguous multi-layer `GetAssayData()` access, the implementation does not score directly on the Seurat object.

Instead:

1. genes from the retained scRef metaprograms, the 3CA metaprograms, and the consensus cell-cycle panel are pooled
2. the required `data.*` layers are read from the merged Seurat object
3. a sparse gene-by-cell matrix is rebuilt for only the genes needed downstream
4. `ScoreSignatures_UCell()` is run on that explicit sparse matrix

The resulting per-cell scores are stored in:

- `sn_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds`

The signature sizes used are written to:

- `sn_outs/Metaprogrammes_Results/UCell_nMP19_filtered_summary.csv`

### 10.5 UCell scoring of 3CA metaprograms

The external 3CA metaprograms are read from:

- `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv`

The metaprogram names are converted to the same naming convention used in the scRef downstream scripts:

- `MP... -> 3CA_mp_...`
- then `make.names()` is applied, which yields names such as `X3CA_mp_30.Respiration.1`

After intersecting signature genes with the merged object, `ScoreSignatures_UCell()` is run again on the same explicit sparse expression matrix.

The resulting per-cell scores are stored in:

- `sn_outs/UCell_3CA_MPs.rds`

The signature sizes used are written to:

- `sn_outs/UCell_3CA_summary.csv`

### 10.6 Merge-stage outputs

The stage also writes:

- `sn_outs/snSeq_malignant_epi_sample_summary.csv`
- `sn_outs/snSeq_malignant_epi_status_input.csv`
- `sn_outs/snSeq_malignant_epi_overall_summary.csv`
- `sn_outs/snSeq_malignant_epi_cc_top50.rds`
- `sn_outs/snSeq_malignant_epi_cc_score.rds`

## 11. State Mapping on the snRNA-seq Malignant Epithelial Dataset

The downstream state scripts are:

- `analysis/cell_states/states_topmpB_reg_noreg.R`
- `analysis/cell_states/states_unresolved_relabel.R`
- `analysis/cell_states/overall_state_proportions.R`
- `analysis/cell_states/sample_abundance.R`
- `analysis/cell_states/states_hybrid_pairwise_nodeplot.R`
- `Auto_cell_states.sh`

For the snRNA-seq workspace, only the `noreg` branch is implemented. The `reg` branch from the scRef workspace is intentionally ignored.

### 11.1 Batch variable

In the snRNA-seq pipeline:

- `reference_batch` is used as the study replacement

This field is used for:

- score normalization
- batch-wise proportion plots
- batch ordering in downstream sample-level visualizations

`technology` is not used as the normalization or grouping variable for state mapping.

### 11.2 Primary state assignment (`states_topmpB_reg_noreg.R`)

The script reads:

- `snSeq_malignant_epi.rds`
- `Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds`
- `Metaprogrammes_Results/UCell_nMP19_filtered.rds`

The retained scRef metaprograms are grouped into the same state structure used downstream:

- `Classic Proliferative = MP2`
- `Basal to Intestinal Metaplasia = MP17 + MP14 + MP5 + MP10 + MP8`
- `Stress-adaptive = MP13 + MP12`
- `SMG-like Metaplasia = MP18 + MP16`
- `Immune Infiltrating = MP15`

Cell-cycle metaprograms are treated as:

- `MP1`
- `MP7`
- `MP9`

and are excluded from the primary state competition.

The cell-cycle score used in the boxplots and heatmap annotations is not recomputed from the merged Seurat assay slots inside the state scripts. Instead it is read from:

- `sn_outs/snSeq_malignant_epi_cc_score.rds`

This score was precomputed in `geneNMF.R` from the merged malignant epithelial `data.*` layers using the consensus cell-cycle panel.

Non-cell-cycle UCell scores are normalized by:

1. centering each score within `orig.ident`
2. dividing by the per-`reference_batch` standard deviation

The primary state call is then made as follows:

- take the maximum adjusted score within each state group
- assign the best-scoring state
- set the state to `Unresolved` if the best score is `< 0.5`
- set the state to `Hybrid` if the gap between the top two state-group scores is `< 0.3`

The outputs are:

- `sn_outs/Auto_topmp_v2_noreg_states_B.rds`
- `sn_outs/Auto_topmp_v2_noreg_mp_adj.rds`
- `sn_outs/Auto_topmp_v2_noreg_group_max.rds`
- `sn_outs/Auto_topmp_v2_noreg_proportion_B_withpie.pdf`
- `sn_outs/Auto_topmp_v2_noreg_ccscore_boxplot_B.pdf`
- `sn_outs/Auto_topmp_v2_noreg_heatmap_B_cconly.pdf`
- `sn_outs/Auto_topmp_v2_noreg_summary.csv`

### 11.3 Unresolved relabelling with retained 3CA metaprograms (`states_unresolved_relabel.R`)

Only cells initially labelled:

- `Unresolved`

are reconsidered.

Among 3CA metaprograms, the cell-cycle-like 3CA programs are excluded before choosing the top per-cell 3CA program.

The snRNA-seq pipeline does not derive a new retained 3CA panel. Instead it uses the retained panel already chosen from scRef:

- `X3CA_mp_30.Respiration.1 -> Classic Proliferative`
- `X3CA_mp_12.Protein.maturation -> 3CA_EMT_and_Protein_maturation`
- `X3CA_mp_17.EMT.III -> 3CA_EMT_and_Protein_maturation`

If an unresolved cell’s top 3CA program is one of these retained programs, the cell is relabelled to the mapped state above. Otherwise it remains `Unresolved`.

The final state labels are stored in:

- `sn_outs/Auto_final_states.rds`

Additional outputs are:

- `sn_outs/task4_unresolved_states/Auto_task4_unresolved_relabel_states.rds`
- `sn_outs/task4_unresolved_states/Auto_task4_unresolved_relabel_mp_coverage.csv`
- `sn_outs/task4_unresolved_states/Auto_task4_unresolved_relabel_proportion.pdf`
- `sn_outs/task4_unresolved_states/Auto_task4_unresolved_relabel_cc_boxplot.pdf`
- `sn_outs/task4_unresolved_states/Auto_task4_unresolved_relabel_heatmap.pdf`
- `sn_outs/task4_unresolved_states/Auto_task4_unresolved_relabel_summary.csv`

## 12. Downstream State Visualisation and Sample-Level Analysis

### 12.1 Overall state proportions (`overall_state_proportions.R`)

This script reads `Auto_final_states.rds` when present, otherwise it falls back to `Auto_topmp_v2_noreg_states_B.rds`.

It writes:

- `sn_outs/Auto_overall_state_proportions.pdf`

### 12.2 Sample abundance plots (`sample_abundance.R`)

This script uses:

- the merged malignant epithelial metadata
- final state labels
- adjusted scRef metaprogram scores
- retained scRef metaprogram definitions

It reproduces the sample-level abundance views for:

- top non-cell-cycle metaprogram
- top metaprogram including cell-cycle programs
- final state labels
- the internal metaprogram breakdown of `Basal to Intestinal Metaplasia`

Sample ordering is generated in two ways:

- by state-diversity-style ranking
- by `reference_batch` then sample name

Outputs:

- `sn_outs/task3_sample_abundance/Auto_task3_sample_abundance.pdf`
- `sn_outs/task3_sample_abundance/Auto_sample_abundance_summary.csv`

### 12.3 Hybrid pairwise plots (`states_hybrid_pairwise_nodeplot.R`)

This script uses the `noreg` state labels and adjusted state-group scores to summarize pairwise hybrid relationships among the five resolved scRef-derived states.

Outputs:

- `sn_outs/task6_hybrid_pairwise/Auto_task6_hybrid_pairwise_nodeplot_noreg.pdf`
- `sn_outs/task6_hybrid_pairwise/Auto_task6_hybrid_pairwise_heatmap_noreg.pdf`
- `sn_outs/task6_hybrid_pairwise/Auto_task6_hybrid_pairwise_nodeplot_summary.csv`
