# Matched PDO and snRNA-seq native CNA/subclone visualization

## Scope and inputs

The terminal script analysis/infercna/compare_matched_pdo_sn_cna.R compares
SUR680T3_PDO with H_post_T1_biopsy, and SUR791T3_PDO with L_post_T1_biopsy.
Use current local manifest/statuses, final malignant epithelial metadata and
centred state labels, and the original per-biopsy InferCNA *_outs.rds matrices.
Only malignant_level_1/2 biopsy cells enter the comparison. The existing PDO
Auto_PDO_cnv_subclone_mp/intermediate/Auto_infercna/<sample>/Auto_cna.rds and
Auto_subclone_calls/<sample>/Auto_call.rds supply native PDO CNA matrices and
production subclone labels. Current PDO centred states supply annotation.

External paths are the live PDOs_Pipeline workspace and the hg38 gene-order file
/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt.
The maintained PDO driver is analysis/cnv/Auto_PDO_cnv_subclone_mp_heatmap.R.
No external file is modified. All new outputs belong to live snSeq sn_outs/
infercna/matched_pdo_sn/{intermediate,tables,figures,logs,reports}.

## Method and provenance

Reuse the existing PDO calls unchanged. Validate their CNA/call fingerprints
and the MD5s of their canonical producer inputs. Load only the exact named
calling functions and parameters from the validated cache/maintained driver
into an isolated environment; do not source the driver or run the external PDO
workflow. Current maintained function bodies must match the cached method.
The source centromere definition and hg38 gene order supply genomic arm labels.
The helper analysis/lib/cna_helpers.R performs this validation and restoration.

For each biopsy, apply that same seeded method to its full malignant-cell CNA
matrix, independently of PDO cells. Keep finite autosomal/chrX genes in genomic
order, select genes at or above the sample's lower-tertile mean absolute CNA
signal, run an exact k=15 nearest-neighbour graph and seeded Louvain (seed 42),
then iteratively merge insufficiently distinct arm profiles. Minimum retained
cluster size is max(20, ceiling(0.05*N)). If no initial cluster meets that size,
the unchanged producer algorithm retains its largest candidate; samples with
fewer than 40 retained cells are one group. This degenerate-case behavior is
reported rather than silently replaced with another calling method.
The exact cached implementation remains authoritative for all merge rules:
arm calls +/-0.10, distinct-arm difference 0.08, at least two moderate arms,
merge maximum difference 0.06, same-pattern correlation >=0.80 and mean absolute
difference <=0.015. No retained cluster is identified through gene-expression
clusters or malignant_clus labels. Small excluded clusters remain explicitly
Not assigned in full-cell denominators; they are not silently discarded.

Calls use the native matrices and gene sets. Display bins use all common finite
genes in the paired matrices and hg38 coordinates, with identical chromosome-
wise groups of 100 genes for both modalities. No per-modality z scaling or
recentring changes the displayed values. Clone consensus values average every
cell in each group; abundance uses every input cell. Cell heatmaps use a fixed
seed, at most 1,200 proportionally stratified cells, ordered within each clone
by Ward clustering on displayed bins. Plot subsampling does not affect calls,
clone means, state assignment or proportions. Full binned matrices, native
calls/methods, exact selected cells and metadata are retained in live caches.

## Figure and interpretation

One PDF page per patient places PDO and biopsy cell heatmaps side by side with
common genomic rows and a shared diverging CNA scale (-0.15, 0, +0.15; display
clipping only). A narrow adjacent panel shows clone-mean profiles. Below are
full-cell clone counts/proportions and a native-state annotation key. PDO and
biopsy clone IDs/colours are independent. Genes/bins, all cell assignments,
consensus profiles and abundance source tables are exported separately.
Within each labelled modality panel, short slice labels A/B/C identify its own
clones and U denotes Not assigned. Full modality-prefixed clone names remain in
the consensus panel, abundance labels, and exported tables. This avoids label
overlap for rare displayed groups without merging or hiding them.

These are RNA-inferred relative CNA signals, not absolute DNA copy numbers or
lineage relationships. PDO uses its existing Carroll reference; snRNA uses its
required source-batch cynthia_sn reference. Do not replace either native reference
or interpret shared colour intensity as calibrated absolute copy number. Similar
profiles suggest descriptive agreement only; clone A in one dataset is not the
same clone as A in the other. The large biopsy and small PDO cell numbers and
the 5% calling threshold limit rare-clone comparisons. Sampling and organoid
culture can alter representation; this figure does not establish selection.

## Execution and caching

Submit job_scripts/infercna/compare_matched_pdo_sn_cna.sh: dmtcp, 4 CPUs, 128 GB,
12 hours. The larger biopsy graph is computed on all malignant cells and can be
slow. Missing/stale inputs stop without fallbacks. plot_only=true validates
input size/mtime/compute signatures and reuses cached computations. The driver,
calling helper, reference files, status files, states and input matrices enter
that signature. Complete run/session logs and a compact figure index are saved.
