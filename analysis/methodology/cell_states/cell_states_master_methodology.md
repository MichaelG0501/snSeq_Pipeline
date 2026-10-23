# Current cell-state PBS entry point

`cell_states_master.sh` runs only
`analysis/cell_states/assign_states_approach_b_noreg.R` after step 7 has prepared
current centred scRef MP scores, grouping, and malignant-cell metadata. See its
mirrored methodology for normalization, thresholds, storage, and validation.

The wrapper uses dmtcp, 4 CPUs, 64 GB, 8 hours, live PBS logging, and explicit
project-root scheduler destinations. R failure stops the wrapper immediately.
Current outputs are under live `sn_outs/cell_states/centred/`. Historical 3CA
unresolved relabelling and historical state plotting scripts are not invoked.
