# Merged CNV Profile Methodology

## Goal

`analysis/infercna/plot_merged_cnv_profile.R` renders a merged InferCNA heatmap profile from final epithelial malignancy outputs.

## Inputs

- Final epithelial objects: `sn_outs/by_samples/<sample>/<sample>_epi_f.rds`
- InferCNA output matrices: `sn_outs/by_samples/<sample>/<sample>_outs.rds`
- Gene order file:
  `/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt`

## Current Safety Rule

The script requires final `_epi_f.rds` objects with a `malignancy` column. It no longer applies fallback malignancy logic from step 5.

## Method

1. Iterate over sample directories.
2. Load final epithelial metadata and InferCNA output matrices.
3. Keep cells present in both the final metadata and InferCNA output matrix.
4. Map final malignancy labels to plotting status groups.
5. Add reference-batch metadata as the study field.
6. Merge epithelial and sampled reference cells.
7. Order genes by chromosome coordinates.
8. Render a ComplexHeatmap CNV profile with chromosome boundaries.

## Outputs

- `sn_outs/infercna/cnv_profile/figures/snseq_merged_cnv_profile.pdf`
- `sn_outs/infercna/cnv_profile/logs/plot_merged_cnv_profile.log`

## Downstream Dependencies

None. This is a terminal figure-generation script.
