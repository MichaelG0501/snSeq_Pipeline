# Matched concordance shared helpers

`analysis/lib/concordance_helpers.R` supports the terminal matched-state workflow.
It performs no standalone file writes and takes objects/paths from its caller.

The helper validates current grouped MP genes, extracts exactly requested cells
from their native RNA count layers, and intersects assayed features across layers
within each sample. Cells must occur exactly once, counts must be finite and
nonnegative, and every cell must have a positive library size. Unmeasured features
are not padded as if observed zeros. Sample identity comes from orig.ident.
The downstream UCell algorithm separately applies its documented/reference-method
zero-expression treatment to missing signature genes; missing-gene coverage is
exported explicitly rather than changing the observed count matrix or gene names.

TMM normalization uses eligible sample/state raw-count pseudobulks separately
within each modality. Log2 CPM uses prior.count=2. A state effect is its logCPM
minus the equal mean of other eligible states in that same sample. At least two
eligible native states are necessary. No model residual variance, moderated t
statistic or replicate-based DGE significance is estimated from single samples.

Pairwise pooled effects average matched patients equally and require exactly the
same contributing patients in both modalities for each compared state pair.
Missing effects stay NA. Signature z scores require at least two finite receiving
states; constant finite score rows become zero, while missing states stay NA.

Overlap tests use one-sided Fisher exact tests on the explicitly supplied gene
universe. Empty sets yield NA Jaccard and p-value. Exact specificity tests enumerate
all permutations of the available complete shared-state submatrix, using matched
mean minus off-diagonal mean. When availability prevents a complete submatrix,
states with the most missing entries are removed deterministically, without
looking at score magnitudes. Fewer than two remaining states yields no test.
These are state-label/gene-set checks, not independent biological replication.

The current caller obtains its fixed transcriptional-effect gene selection from
the original PDO/scRef result cache, then restricts that exact ordered vector to
genes assayed in all matched samples. It does not use the helper's expression
screen to select or prune this effect pool. The original 3,538-gene list and
available 3,268-gene subset are exported independently of marker ranking.
The 20-cell threshold applies only to effect pseudobulks; current UCell receiving
state means include all observed confident native groups. Page 4 now projects
current individual PDO MPs onto snRNA, using reference_batch normalization.

See `analysis/methodology/cell_states/compare_matched_pdo_sn_state_concordance_methodology.md`
for full inputs, interpretation, execution and output contracts.
