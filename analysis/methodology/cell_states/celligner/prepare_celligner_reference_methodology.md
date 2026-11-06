# Public Celligner reference preparation

Script: `analysis/cell_states/celligner/prepare_celligner_reference.R`.
Status: active input preparation for reference-based concordance analysis.

## Scope and reference provenance

The user authorized obtaining the public pan-cancer reference needed to assess
PDO–biopsy concordance with Celligner. This workflow downloads the original
authors' Figshare release, article 11965269 version 5. It includes aligned
expression for 13,485 samples and 19,188 Ensembl genes, sample metadata, five cPC
loading vectors, original HGNC mappings and the release README. Published counts
are 12,236 tumours and 1,249 cell lines. Source:
[original dataset](https://figshare.com/articles/dataset/Celligner_data/11965269)
and [original Celligner study](https://doi.org/10.1038/s41467-020-20294-x).

The 2026 HCMI paper uses a newly fitted reference and additional HCMI integration.
Its deposited aligned expression, cPC loadings and full 70-PC distances are
controlled-access according to the [GDC publication listing](https://gdc.cancer.gov/about-data/publications/HCMI-CMDC-2026).
The original public reference acquired here is a distinct, fully documented
source, not the exact HCMI result. The author's historical pretrained model URL
returns HTTP 403; no attempt is made to bypass access controls. The
[HCMI pipeline](https://github.com/NCICCGPO/HCMI-Flagship-Manuscript/blob/main/pipelines/pipelines_launchers/run_pipeline_celligner.py)
and official Celligner source are public.

## Download and validation

Query the immutable version-specific Figshare API manifest and select only the
five named essential files. Download directly into the stable live reference
directory. Existing downloads are reused only when the published size and MD5
agree; partial downloads can resume through curl. Oversized or checksum-invalid
files stop the workflow without deletion. Do not silently switch data versions.

Read the expression matrix on PBS, check its exact dimensions, unique sample and
gene IDs, numeric finite values, agreement with all metadata sample IDs, and
published tumour/cell-line counts. Match metadata by explicit sample IDs.
Validate that cPC gene IDs exactly match expression columns and reorder cPC rows
explicitly to that order. Verify HGNC symbol and Ensembl mapping columns.
Preserve the exact original source files and the original manifest.

Prepare an uncompressed canonical RDS containing aligned expression, metadata,
cPCs, HGNC mapping and file/parameter fingerprints. It is an input cache, not an
ephemeral or archive file. Reuse it only with an identical verified fingerprint.
Preserve the official R methods/helpers and README from commit
`d9c9246f8a1b6885d07f2f28bbdca24253e57cf1` with a source manifest. No package
installation or external-workspace modification occurs in this preparation.

## Outputs, execution and downstream use

All outputs are under `sn_outs/cell_states/matched_pdo_celligner_concordance/`:

- `reference/original_celligner/`: original source data and pinned source code.
- `intermediate/public_celligner_reference.rds`: canonical prepared input.
- `tables/reference_file_manifest.csv` and `reference_software_manifest.csv`.
- `reports/reference_validation.json` and `reference_acquisition.md`.
- `logs/reference_preparation.log`: start/end, parameters, inputs, outputs,
  cache reuse, session information and failures when applicable.

Run `qsub job_scripts/cell_states/celligner/prepare_celligner_reference.sh` from
the live root. The ignored minimal launcher uses dmtcp, 2 CPUs, 32 GB, two hours,
`-koed` and explicit project-root PBS stdout/stderr paths. Internet and curl are
required. All substantial download, checksum and matrix parsing occurs on PBS.
No plotting option is necessary for deterministic reference preparation.

The reference is intended for a subsequent, explicitly specified Celligner
projection. A published aligned matrix is not a full serialized fitted model;
the projection must initialize clustering, differential-feature and cPCA/MNN
state correctly. Local UMI pseudobulk normalization also needs explicit
documentation because UMI counts are not bulk TPM. This workflow does not infer
concordance, calibrate significance, or tune parameters to true matched pairs.
