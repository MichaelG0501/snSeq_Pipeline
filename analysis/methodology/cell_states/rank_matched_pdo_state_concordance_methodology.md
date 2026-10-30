# Relative concordance of matched biopsies against the PDO cohort

## Question and prespecified state correspondence

This terminal analysis asks whether H_post_T1_biopsy/SUR680 and
L_post_T1_biopsy/SUR791 match their own PDO sample more closely than other PDO
samples when comparing the corresponding transcriptional state. The true pairs
are SUR680T3_PDO and SUR791T3_PDO. This extends, without replacing, the previous
matched-pair concordance reports. Both biopsies share one selected combined PDF.

Prespecified native snRNA to native PDO correspondence is Classic proliferation
to Classic proliferation, Squamous-to-intestinal to Columnar-to-intestinal,
Glandular-to-intestinal to Glandular differentiation, and Stress-adaptive to
Stress-adaptive. Immune mimicry has no prespecified PDO counterpart and remains
explicitly unavailable rather than being forced to another state. Native
Hybrid/Unresolved cells remain in source counts but do not define state effects.

## Exact live inputs

Local inputs under snSeq_Pipeline/sn_outs are snSeq_malignant_epi.rds,
cell_states/centred/intermediate/centred_refined_noreg_states.rds,
cell_states/centred/tables/centred_refined_mp_state_grouping.csv,
sample_manifest.csv, and each matched biopsy's expression-filter, infercna and
malignancy status files. All three status gates must be ok. Only final malignant
level-1/2 epithelial biopsy cells enter the analysis.

External inputs under live PDOs_Pipeline/PDOs_outs are PDOs_merged.rds,
centred_mp_refinement/centred_refined_noreg_states.rds, its current MP/state
grouping table, and Auto_PDO_scRef_centred_state_concordance/intermediate/
Auto_state_concordance_results.rds. The shared live PDO analysis config supplies
current native state names and the established SUR843T3_PDO exclusion. The
current merged object has 24 included PDO samples and 72,801 cells; the excluded
sample is not present in this object. Scope is discovered from orig.ident,
never by trimming barcodes. Missing required canonical labels stop execution.

All external inputs are read only. Nothing from the forbidden upstream derived
sn annotation workspace is used. No network access, package installation, UCell
rescoring, new state calling or canonical object modification is required.

## Fixed comparison genes and effect estimation

Read the original PDO/scRef derived$informative_genes vector (3,538 genes) from
the documented reference cache. Preserve its order and restrict only to genes
assayed in every included PDO sample and both biopsies. Export full per-sample
availability and the actual comparison list, expected to contain the same 3,268
available genes as the previous matched analysis. Do not reselect genes from
rank results, filter by detected expression, map aliases or pad unavailable genes
with expression zeros. Assayed zero counts remain genuine measured zeros.

Read current raw RNA counts and fixed native labels. Sum counts within each
sample/state using all its cells; retain effect pseudobulks with at least 20 cells.
At least two eligible native states are needed for a state-vs-rest contrast.
The shared concordance_effects helper performs TMM normalization and log2CPM
with prior.count=2. Subtract the equal mean of the other eligible native states
within the same sample. No biological-replicate DGE significance is inferred.

All cohort PDO pseudobulks are normalized together, including the two true
matches; both biopsy samples are normalized together separately. This applies
identical processing to every PDO candidate. Expanding the PDO normalization
cohort can slightly change coefficients from the earlier two-PDO-only run;
do not mix cached old matched-PDO vectors with newly normalized alternatives.
The complete count matrices, normalization factors and all-gene effects remain
in the computation cache. The shared helper's historical patient-name key is
adapted internally to actual sample IDs so repeated cultures remain distinct.

## Relative ranking and eligibility

For each of ten biopsy/state queries, compute Spearman rho between its effect
vector and each PDO's corresponding-state effect vector on the identical ordered
comparison genes. All samples remain visible as columns. Missing counterparts,
undersized states, unavailable contrasts and constant vectors yield explicit
NA/reason fields, never zero correlations.

Rank only eligible PDO samples within each query. Highest rho is rank 1;
scores within 1e-12 receive the same minimum rank and all best ties are marked.
The denominator counts eligible PDO samples including the true match. Matched
percentile measures how often the matched rho exceeds other eligible PDO rhos,
giving half credit to ties. The matched-minus-best-other margin is positive
when the true matched PDO exceeds every alternative, zero for a top tie, and
negative when another PDO scores higher. Missing scores do not enter denominators.

The fixed 20-cell gate makes SUR791 biopsy proliferation unavailable (15 cells).
Both immune-mimicry queries lack a corresponding PDO state. Thus seven shared
state queries are expected to support relative ranking and scatter examples.
Source tables retain all ten queries. Figures omit the two unmapped immune
mimicry rows, retaining eight rows and the unsupported SUR791 proliferation row.

## Selected visualization

One final PDF combines patient-specific comparisons while preserving identity:

1. A heatmap labels its rows "snSeq state", with all PDO samples as columns and
   shared raw-rho colours. Each row uses its corresponding PDO state. Black boxes
   (linewidth 0.7) mark true matched PDOs; red boxes (linewidth 1.4) mark the
   highest rho, including ties. No asterisks are shown. Coincident boxes retain
   visible red edges around the thinner black line. Unavailable cells are grey
   and row labels state the reason. Column order stays fixed.
2. A rank-context plot shows every eligible PDO as a grey point, the true match
   as an orange diamond, and "Highest matching" as a blue point. Row notes
   show matched rank/eligible denominator and the matched-minus-best-other gap.
3. Subsequent example pages, one per supported biopsy/state query, compare its
   "Matched PDO", "Highest matching" and "Median matching". Highest matching
   means the highest-scoring other PDO; Median matching means the distinct
   remaining sample closest to the cohort median rank.
   Tied example selections are alphabetically deterministic. Every score remains
   available in tables; examples are not a significance-selected subset.

Scatter pages use identical genes and symmetric axes within each page, with
snRNA effect on x and PDO effect on y, a descriptive regression line and up to
six strong same-direction gene labels. Figures keep all comparisons on their
intended page with readable labels. Page-footer paragraphs are omitted; method
details remain in this methodology and the current-results report. No duplicate
standalone PDFs are exported.
The 18-by-10-inch pages use approximately 13.6-point cell/rank/gene annotations
and 13–17-point axis/sample text, providing roughly 10-point annotation text when
placed at full width on a standard widescreen slide. Font increases do not add
pages or split the comparisons. Final source tables retain full precision.

## Interpretation limits

Ranks describe this finite PDO-sample cohort, not a random-patient validation
population. Some donors contribute treated and untreated cultures, so candidates
are not independent biological patients. A high relative rank need not imply
high absolute concordance. No population significance, selection or common
clonal ancestry is established. Native state-vs-rest backgrounds differ when
states fail the cell-count gate; all cell counts and numbers of eligible native
states are exported to make that limitation visible. Platform, culture,
treatment and sampling differences can influence both correlations and ranks.

## Execution, stable outputs and reproducibility

Submit qsub job_scripts/cell_states/rank_matched_pdo_state_concordance.sh:
dmtcp, 2 CPUs, 32 GB, 1 hour, explicit live-root -koed logs. The maintained R
script contains all analytical logic and depends only on tracked shared helpers,
not on scratch audits or its launcher. All heavy work runs on PBS.

Stable root: sn_outs/cell_states/matched_pdo_relative_concordance/.
The single selected PDF is figures/Auto_final/01_matched_PDO_relative_state_concordance.pdf.
Tables contain all candidate scores/reasons/ranks, matched-rank summaries,
selected examples, native state counts, original/available gene lists and wide
PDO/sn effect vectors with sample|||state column names. Intermediate caches
retain pseudobulks, normalization, full effects, ranks and input fingerprints.
Reference/ retains the original gene vector. Reports/ has a linked current
result index, page index, independent checks and preview renders; logs/ records
timestamps, exact inputs/outputs/parameters, cache use and session information.

plot_only=true validates unchanged input paths/sizes/mtimes, shared helper files,
PDO config and parameters before regenerating plots without RNA reload. These
fingerprints cannot detect edits preserving both file size and mtime. Finish
source changes before submitting; never edit code while a running R job reads it.
