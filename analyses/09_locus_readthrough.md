# Strand-specific ivory locus coverage

**Paper:** Figure 6A.

**Status:** cleaned scripts supplied; final inputs/settings and result matching pending.

Recover per-barcode exonic/intronic coverage and the locus ribbon plot.

## Reusable single-cell plotting tutorial

For a guided workflow that can be applied to other loci, start with **[the scsplice walkthrough](../tools/scsplice/WALKTHROUGH.md)**. The toolkit accepts an indexed BAM and existing sample/barcode/cell-type metadata, extracts the selected region, and produces labelled read-start ribbons and per-cell aligned-block/splice-gap panels. See the [command reference](../tools/scsplice/README.md) for the preserved historical binning and explicit display options.

This reusable tool has synthetic test coverage. Matching the final Figure 6A inputs and biological results remains pending, as described below.

## Inputs

- Aligned BAMs with retained nucleus barcode tags and indices.
- Final plotted barcode selections, strand convention, annotation and plotting coordinates.
- Existing strand-separated CortexLocus BAMs for the historical script below.

## Outputs

- Historical strand-specific spliced BAMs; final per-cell coverage tables and Figure 6A vector plot still needed.

## Software

samtools, awk and original coverage/plotting software; actual versions pending.

## Execution and interpretation

The historical list assignment is commented, so its loop does nothing unless list is already set in the shell. Earlier extraction, filtering and strand splitting are commented as well. The active commands retain reads whose CIGAR contains N from pre-existing strand-specific BAMs. This is only part of the coverage workflow, and it does not distinguish all exonic/intronic coverage or make the final ribbon plot.

See [intake review](../docs/INTAKE_REVIEW.md) for unresolved decisions and [file requests](../docs/FILE_REQUESTS.md) for the complete acquisition checklist.

## Cleaned code and execution

- [extract_spliced_reads.sh](../scripts/locus/extract_spliced_reads.sh)

This runnable script extracts spliced alignments from supplied strand-separated BAMs. It does not fill in the missing barcode-selection or complete coverage/ribbon pipeline.

See [run commands and input formats](../scripts/README.md), [cleanup notes](../docs/CLEANUP_NOTES.md) and the [original source archive](../archive/README.md).

[Return to analysis index](../README.md#analysis-workflows).
