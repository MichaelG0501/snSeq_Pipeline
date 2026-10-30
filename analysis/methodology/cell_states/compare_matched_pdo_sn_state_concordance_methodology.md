# Matched PDO / snRNA centred-state concordance

## Scope, reference and approved interpretation

The terminal analysis compares SUR680T3_PDO with H_post_T1_biopsy and
SUR791T3_PDO with L_post_T1_biopsy, plus a combined report. Its page order and
concordance evidence follow the maintained read-only reference:
/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/cell_states/Auto_PDO_scRef_centred_state_concordance.R.

The user explicitly requested relative effects without biological-replicate
significance, the original PDO/scRef pooled gene list, and skipping genes absent
from snRNA annotation. The reference supplies gene selection, while the matched
samples supply their own descriptive state-vs-rest expression effects. No
matched-data reselection of informative genes is performed. Neither cells nor
state pseudobulks are treated as independent patient replicates.

## Exact live inputs and current MP verification

Local inputs under snSeq_Pipeline/sn_outs are snSeq_malignant_epi.rds,
sample_manifest.csv, matched expression-filter/infercna/malignancy statuses,
cell_states/centred/intermediate/centred_refined_noreg_states.rds,
cell_states/centred/intermediate/merged_refined_mp_genes.rds and
cell_states/centred/tables/centred_refined_mp_state_grouping.csv.

Current scRef genes and grouping are read directly from the live scRef_Pipeline/
ref_outs/Metaprogrammes_Results/centred/mp_refinement/{intermediate,tables}/
merged_refined_mp_genes.rds and centred_refined_mp_state_grouping.csv respectively.
Before scoring, all 17 named gene vectors and the mp/state grouping must exactly
match the local sn copies after aligning MP order. Cell-cycle MPs remain
excluded from state calling. This verifies that sn states use the current scRef
signatures, rather than a separately derived sn state panel.

PDO inputs under live PDOs_Pipeline/PDOs_outs are PDOs_merged.rds and
centred_mp_refinement/{centred_refined_noreg_states.rds,
merged_refined_mp_genes.rds,tables/centred_refined_mp_state_grouping.csv}.
The PDO shared configuration supplies native colours/canonical paths. All
external files are read only and no network download is required.

## Frozen reference effect-gene pool

Read derived$informative_genes from live PDOs_Pipeline/PDOs_outs/
Auto_PDO_scRef_centred_state_concordance/intermediate/Auto_state_concordance_results.rds.
Validate its exact ordered reconstruction from the original PDO/scRef DGE tables:
union each state's top 500 absolute moderated-t genes, restricted to the original
common tested genes. Record the source cache MD5 and freeze the full original
vector under the current output reference/ tier.

The original vector contains 3,538 unique genes. Of these, 270 are absent from
both biopsy count-layer annotations, leaving 3,268 assayed genes in all four
matched datasets. Keep this exact ordered assayed intersection for every
supported correlation and scatter in all three reports. Do not drop assayed
zero-count genes, refilter by detection, map aliases or substitute zeros for
unassayed genes. Export all original/retained/skipped genes with per-sample
availability. The 3,268 value is a verified current result, not a hard-coded
selection quota; changed inputs invalidate the cache.

## Effects, markers and pooled summaries

Sum raw RNA counts within sample and native state. Groups require at least 20
cells; at least two eligible native states are necessary for a within-sample
contrast. TMM-normalize eligible pseudobulks separately within each modality,
then compute log2CPM with prior.count=2. For each native state, subtract the equal
mean of the other eligible native-state logCPM vectors from its own vector.
These are descriptive relative expression effects. No moderated DGE t statistic
or biological-replicate significance is estimated for these matched samples.

Correlate effect vectors by Spearman rho on the fixed available reference pool.
For marker overlap, independently rank positive effects over common assayed
genes and take at most 200 per state; report Jaccard and contextual one-sided
Fisher/BH statistics. Those marker statistics do not define or alter the effect
correlation pool. Missing or undersized state effects stay NA, never zero.

The combined effect comparison averages within-sample effect vectors equally
over patients with both compared native states supported. Identical contributing
patients enter both modalities for each state pair; matrices expose n=1 versus
n=2. The pooled UCell means and projected proportions pool native-state cells,
matching the original reference. There is no extra weighting-sensitivity page.

## Reciprocal signatures and PDO-to-sn diagnostic labels

Group the current MPs into native state-signature gene unions independently
for PDO and scRef. Score both sets and current individual MPs with UCell 2.2.0,
maxRank=1500, on each full assayed RNA count matrix. Never rank signature-only
expression. Coverage and exact missing signature genes are exported. UCell's
original missing-signature behavior is retained; the assayed-only effect-gene
restriction is a separate analysis rule.

Signature means include every observed confident native state, including the
seven-cell PDO motile-cilia group and fifteen-cell sn proliferation group.
The 20-cell gate applies only to expression-effect pseudobulks. Standardize
within each source signature across receiving native-state means. Page 2 shows
PDO state signatures by sn receiving state, current scRef state signatures by
PDO receiving state (with its original source/target orientation), and direct
MP-state Jaccard overlap. That fixed overlap matrix is shared reference context,
not independently replicated patient evidence.

For page 4, take current individual PDO MP scores on confident-native sn cells,
centre by orig.ident, divide by reference_batch SD across the matched confident
sn cohort, maximize within each PDO state and apply the fixed rules: best<0.5
is Unresolved; otherwise top-two gap<0.3 is Hybrid. This uses shared checked
state helpers and does not modify canonical native states. Native Hybrid and
Unresolved cells are excluded from receiving groups, exactly as in the original
reference; projected Hybrid/Unresolved cells remain in every displayed denominator.
All original cells remain in the audit tables. The explicit four-state
correspondence supplies match labels; unmatched native immune mimicry is not
forced to a PDO state.

Prespecified PDO-to-sn correspondences are Classic proliferation to Classic
proliferation, Columnar-to-intestinal to Squamous-to-intestinal, Glandular
differentiation to Glandular-to-intestinal, and Stress-adaptive to Stress-adaptive.
PDO ECM-remodelling and Motile-cilia differentiation remain alternative states.
Expected-match margins and exact label permutations are conditional descriptive
checks, without claims of replicated biological significance.

## Ten-page presentation contract

Each of SUR680, SUR791 and pooled has exactly ten PDF pages:

1. Transcriptional-effect Spearman correlation and positive-effect marker overlap.
2. PDO-to-sn and current-scRef-to-PDO signature projections, plus MP gene overlap.
3. Prespecified-match margin heatmap and expected-versus-alternative specificity.
4. Current PDO MPs projected to native sn states, as percentage stacked bars.
5. Classic proliferation PDO state versus all five sn states.
6. Columnar-to-intestinal PDO state versus all five sn states.
7. Glandular differentiation PDO state versus all five sn states.
8. Stress-adaptive PDO state versus all five sn states.
9. ECM-remodelling PDO state versus all five sn states.
10. Motile-cilia differentiation PDO state versus all five sn states.

Scatter orientation matches the reference: sn effect on x, PDO effect on y.
All 30 supported state pairs in a report use the identical available reference
genes and common symmetric axes. Draw a descriptive regression line and label
up to eight strongest same-direction genes, ranked by absolute effect sum;
no moderated-t ranking is invented. Unsupported comparisons are explicitly
marked unavailable. Heatmaps show actual values and complete expected-match
outlines; page assembly uses grid/gridExtra to avoid local patchwork incompatibility.

## Execution, outputs and validation

Submit qsub job_scripts/cell_states/compare_matched_pdo_sn_state_concordance.sh:
dmtcp, 2 CPUs, 32 GB, 1 hour. All analytical work runs on PBS with live -koed logs.
Stable output root: sn_outs/cell_states/matched_pdo_sn_concordance/.
The figures tier contains exactly the three selected named PDFs. Tables contain
complete cell audits, signature coverage, frozen/available genes, effect/marker
matrices and vectors, directional signature matrices, matched-patient support,
contextual specificity and *_pdo_signatures_on_sn_concordance.csv. Intermediate
caches preserve all scores, pseudobulks, effects, projection labels and input
fingerprints. Reports contain relative-linked current results, provenance,
MP identity, integrity checks and independent PDF/source-data validation.

Set plot_only=true through the launcher to reuse a signature-matched cache;
it performs no RNA reload or UCell scoring. Signatures include current raw/state
inputs, official scRef source definitions, local copies, reference gene-pool
cache, helpers and parameters. Styling edits must finish before job submission.
Source-size/mtime fingerprints cannot detect edits preserving both values.
No files are deleted. Obsolete tables from the previous projection direction
are retained for safety and are not listed as current result inputs.
