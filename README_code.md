# Code guide

The main working code is in [scripts/](scripts/README.md): cleaned and annotated Bash/R entry points with input checks, explicit parameters, output tables and vector PDFs. [CLEANUP_NOTES.md](docs/CLEANUP_NOTES.md) records repairs, adaptations, tests and limits. The [paper analysis index](README.md#analysis-index) links the annotated workflow descriptions.

Original scripts are preserved under `archive/original/`; received text attachments are under `archive/received/`. Their hashes and transcription/extraction steps are recorded in [received_source_manifest.tsv](metadata/received_source_manifest.tsv). These originals are evidence for the cleaned versions and are not the recommended entry points.

The cleaned scripts separate coherent steps from interactive session histories. Numerical choices, assay choices, sample metadata and contrast membership are explicit. Historical reconstruction examples are labelled; an incomplete exploratory branch is not presented as final science. Once final inputs and parameter choices arrive, the remaining work is to run and compare the scripts with published outputs.

RDS is the default object format; H5Seurat inputs and optional integration export are supported. Scripts emit settings and sessionInfo; use [envs/README.md](envs/README.md) to supply the actual environment. Working-directory paths are documented in the run guide. The example TOML is a planning aid, not a hidden dependency.

Use [FILE_REQUESTS.md](docs/FILE_REQUESTS.md) for missing inputs and [INTAKE_REVIEW.md](docs/INTAKE_REVIEW.md) for final scientific choices still to resolve. The earlier `analysis/` planning folders redirect to `analyses/` to preserve existing links.
