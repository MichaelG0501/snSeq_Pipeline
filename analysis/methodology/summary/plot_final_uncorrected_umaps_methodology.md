# Full-cell uncorrected UMAP visualization

## Existing cohort map

Step 4 already computes a full filtered-cohort UMAP in Expr_filtering.R from
variable genes, scaled expression and PCA dimensions 1:20. It does not run
Harmony, integration, batch covariate regression, or batch correction. Its saved
sn_outs/snSeq_merged.rds contains the complete embedding; snSeq_filtered.png and
summary/expression_filter/figures/pair3_umap_celltype_technology.png already show
it. The new script reuses that exact embedding without recomputing or subsampling.

## Final malignant map

The final malignant-only merge has no cohort-wide embedding. Build one for every
malignant_level_1/2 cell in snSeq_malignant_epi.rds. Assemble all normalized RNA
and count layers with exact barcode matching. Use existing log-normalized data;
select 3,000 variable genes from counts, scale those genes with no covariate
regression, run 30 PCA dimensions and UMAP on dimensions 1:20, seed 1. No batch
correction, integration, neighbour-cell subsampling or state-score normalization
enters this embedding. State labels come from the current centred state vector;
those labels annotate the map but do not define its geometry.

## Outputs and execution

Run job_scripts/summary/plot_final_uncorrected_umaps.sh through PBS: dmtcp,
4 CPUs, 128 GB, 8 hours. All results live under
sn_outs/summary/final_uncorrected_umaps/ output tiers. The persistent cache stores
full-cohort coordinates, malignant coordinates/metadata, malignant PCA and UMAP
reductions, variable genes, input file size/mtime fingerprints and compute
parameters. Coordinate CSVs preserve every barcode and sample/batch annotation.
The final two-page PDF shows full-cohort cell type/batch and final malignant
state/batch. Every cell is plotted. Plot order is shuffled deterministically to
reduce layering bias, without changing counts. Batch separation is retained.
Point size 0.40 and alpha 0.90 preserve visibility when the full vector PDF is
rendered for slides; every plotted point still represents an observed cell.

Use qsub -v plot_only=true job_scripts/summary/plot_final_uncorrected_umaps.sh to
reuse unchanged caches for styling. Neither original Seurat input is overwritten.
Count/coordinate integrity, parameter signatures and complete state coverage are
checked before export; run logs include session information. Source file
fingerprints cannot detect content edits that preserve both size and mtime.
