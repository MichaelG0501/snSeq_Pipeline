# Shared Analysis Library Methodology

## Goal

The `analysis/lib/` files centralize constants and reusable helper functions so workflow scripts do not copy/paste state definitions, thresholds, status readers, logging helpers, or plotting defaults.

## Files

- `analysis/lib/config.R`
- `analysis/lib/io_helpers.R`
- `analysis/lib/logging.R`
- `analysis/lib/state_helpers.R`
- `analysis/lib/plot_helpers.R`

## Configuration

`config.R` defines:

- project paths
- output-tier helpers
- preferred state definition: Approach B, `noreg`
- metadata column names
- QC and state thresholds
- plot defaults for PPTX-readable figures
- scRef MP descriptions
- state groups and colors
- fixed retained 3CA state map

## I/O Helpers

`io_helpers.R` defines:

- status-file reading
- optional CSV loading
- required-file checks
- checked CSV and RDS writers
- Seurat v5 layer matrix extraction for targeted gene sets

## Logging

`logging.R` defines lightweight run-summary functions. Calling scripts should write logs under their `logs/` tier and include:

- script path
- start/end time
- duration
- inputs
- outputs
- parameters
- cache reuse status
- session information

## State Helpers

`state_helpers.R` defines:

- sample-centered / `reference_batch`-scaled score normalization
- 3CA label cleanup
- MP labelling
- final state ordering and color resolution
- state-group maximum score calculation
- preferred final state loading

## Plot Helpers

`plot_helpers.R` defines presentation-oriented ggplot/heatmap text defaults. New figure scripts should prefer these defaults or explain deviations in their methodology file.

## Usage

Workflow scripts should source only the helper files they actually need. Shared helpers are not standalone analysis workflows and should not write outputs by themselves.
