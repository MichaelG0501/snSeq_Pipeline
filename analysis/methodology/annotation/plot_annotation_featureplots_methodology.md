# Annotation Featureplot Methodology

## Goal

`analysis/annotation/plot_annotation_featureplots.R` renders annotation review PDFs that match the current implemented logic in `Clustering.R`, `Annotation.R`, and `Expr_filtering.R`.

## Inputs

- `sn_outs/sample_manifest.csv`
- `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Marker_Genes.xlsx`
- `sn_outs/annotation_score_markers.csv`
- `sn_outs/expression_filter_markers.csv`
- `sn_outs/by_samples/<sample>/<sample>_anno.rds`

## Sample Selection

Default mode is:

- `sample_mode=subset`
- `n_samples=10`
- `sample_seed=1`

Parameters can be passed as command-line key/value pairs or environment variables:

- `sample_mode=all`
- `sample_mode=<sample_id>`
- `n_samples=<N>`
- `sample_seed=<seed>`

## Method

1. Load the 3CA marker workbook and recreate the top-50 initial marker scores.
2. Load refined annotation markers and expression-filter markers from current pipeline outputs.
3. For each selected annotated sample object:
   - verify UMAP availability
   - calculate initial 3CA score matrices
   - calculate refined weighted score matrices
   - calculate canonical marker expression summaries
   - calculate expression-filter panel scores
   - derive filter-status annotations
4. Render four pages per sample:
   - initial 3CA scores
   - refined weighted scores
   - canonical marker expression
   - expression-filter panel scores

## Outputs

- `sn_outs/snseq_annotation_featureplots_<selection>.pdf`
- `sn_outs/annotation_featureplots_status_<selection>.csv`
- `sn_outs/selected_samples_<selection>.csv`
- Tiered copies under `sn_outs/annotation/featureplots/`

## Downstream Dependencies

None. This is a terminal review script.
