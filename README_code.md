# Code guide

The reader-facing workflows are the numbered Markdown documents in [analyses/](analyses/), linked from the [analysis index](README.md#analysis-index). Their Inputs, Outputs, Software, Execution and interpretation, and Code sections follow a consistent format.

## Source treatment

Original scientific scripts are kept under `archive/original/`; received text attachments are under `archive/received/`. Their provenance and transformations are recorded in [received_source_manifest.tsv](metadata/received_source_manifest.tsv). Workflow code fences reproduce complete supplied files or exact line excerpts. No missing parentheses, thresholds, object names or session setup have been guessed or silently repaired.

The long R paste is a session history with exploratory branches and prose, not a runnable final script. Its excerpts are labelled accordingly. Once final branches and inputs are confirmed, portable scripts can be added alongside these annotated documents, with every adaptation recorded.

## Working directories and dependencies

Historical scripts read relative paths and assume objects or variables from earlier sessions. Read each workflow before copying a block into R or Bash. Record the intended working directory, input checksums, sample ordering, active assay, package versions and random seeds when recovering a run. The example path configuration is a planning aid; original scripts do not consume it.

Use [envs/README.md](envs/README.md) to supply actual environment exports. Guessed lockfiles should not replace executed environment records.

## Collection and review

See [FILE_REQUESTS.md](docs/FILE_REQUESTS.md) for the full intake checklist and the collection helper usage. [INTAKE_REVIEW.md](docs/INTAKE_REVIEW.md) documents unresolved thresholds, metadata, statistics, annotation and enrichment choices. [FIGURE_MAP.md](docs/FIGURE_MAP.md) connects received source to paper panels.

The earlier `analysis/` planning folders now point to the numbered workflows. They are retained to preserve existing links; `analyses/` is the main entry point.
