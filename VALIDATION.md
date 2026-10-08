# Validation and provenance

## Checks performed

- The supplied source snapshots are tracked by SHA-256 in [SOURCE_SHA256SUMS.txt](SOURCE_SHA256SUMS.txt) and the [received source manifest](metadata/received_source_manifest.tsv).
- The two supplied shell scripts pass Bash syntax checking. This does not validate their inputs, scheduler settings, output files or scientific results.
- Numbered workflow code excerpts are copied directly from their linked source ranges. Source-needed workflows contain no invented analysis code.
- Relative links in the new reader-facing Markdown documentation were checked locally.

## Limits

No scientific analysis has been rerun. The full dataset, final object builder, environment records and final plot/source tables are incomplete. R was unavailable in the intake workspace. The long R session paste includes syntax errors and exploratory branches, so it is not claimed to pass parsing or execution. The cleaned R scripts now separate those steps but still need R/package execution checks in the actual environment. Shell orchestration checks with mock tools passed; see [cleanup notes](docs/CLEANUP_NOTES.md).

## Required release validation

Recover final inputs and environments; reconcile sample stages, integration parameters, high/low thresholds, DE methods and GO correction settings; run each selected workflow; compare numerical outputs and vector plots with final figure source data. Record exact commands, versions, seeds, exclusions and checksums. Then create a versioned archive and add its DOI to the citation metadata.

See [INTAKE_REVIEW.md](docs/INTAKE_REVIEW.md), [REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) and [PROVENANCE.md](docs/PROVENANCE.md) for the detailed evidence and remaining work.

## miR-193 predictions and comparative annotations

The executed parser reproduces 2,079 predictions across 1,265 genes and the supplied gene summaries exactly. Comparative annotations reproduce all 733 S8 hit counts/Boolean fields and 15 conserved candidates. One missing S8 energy is reported without modifying the original. Twelve text checksums and mock launcher checks pass. RNAhybrid and R scientific workflows were not rerun. See [target validation](docs/TARGET_VALIDATION.md) for scope and binary-input limitations.
