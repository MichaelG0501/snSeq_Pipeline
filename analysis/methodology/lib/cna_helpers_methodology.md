# Shared native CNA comparison helpers

The helper is consumed by analysis/infercna/compare_matched_pdo_sn_cna.R; the
complete scientific thresholds, inputs and limitations are documented in
analysis/methodology/infercna/compare_matched_pdo_sn_cna_methodology.md.

It restores the exact named calling functions recorded in validated production
PDO Auto_call.rds files into an isolated environment. Every restored function
must also match the maintained PDO driver's parsed function definition. Only
specified top-level assignments are evaluated: calling functions and the
centromere vector. The external driver is never sourced or executed, so loading
its method cannot run PDO analyses or write external outputs. Parameters come
from the cached method; gene coordinates come from the explicit hg38 input.
The calling script separately validates original producer input MD5s.

Shared display bins group the same ordered genes into chromosome-wise sets of
100 genes, retaining shorter terminal bins. Both matrices use identical row
indices; arithmetic cell-wise means preserve native signal units. Plot-cell
selection uses fixed seed 42 and proportional allocation with at least two cells
per nonempty group where the total permits. Sampling changes display only; all
cells remain available for clone consensus, assignments and abundance. The
helper has no standalone outputs and must not create directories or run jobs.
