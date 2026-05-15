# Malignancy Signature Validation Methodology

## Goal

`analysis/infercna/validate_malignancy_signatures.R` validates cancer-signature and cell-cycle scores against InferCNA classifications before final malignancy labelling.

## Inputs

- `sn_outs/cancer_signatures.txt`
- `sn_outs/cancer_signatures_summary.csv`
- Per-sample step-5 epithelial objects: `<sample>_epi.rds`
- Cell-cycle reference: `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv`

## Method

1. Load the global cancer signature genes from compiled InferCNA summaries.
2. Select samples with signature support.
3. Load each step-5 epithelial object.
4. Compute cancer signature score from RNA data.
5. Compute cell-cycle score using consensus cell-cycle genes.
6. Compare scores against cell-level InferCNA classifications and cluster-level consensus labels.
7. Render violin/boxplot panels into a validation PDF.

## Outputs

- `sn_outs/infercna/signature_validation/reports/validation_cancer_signatures_vs_cna.pdf`
- `sn_outs/infercna/signature_validation/logs/validate_malignancy_signatures.log`

## Downstream Dependencies

None. This is a terminal audit/report.
