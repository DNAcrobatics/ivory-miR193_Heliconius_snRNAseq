# Reference genome and annotation

**Paper:** Alignment and all gene-based analyses.

**Status:** Source needed.

Prepare the exact Hmel reference and edited ivory annotation used for counting.

## Inputs

- Hmel genome FASTA and mitochondrial sequence, with accession/version and checksums.
- `Hmel.v3.1.annotation.CAT.MOD.plusMito.reducedIvory.gtf`; annotation edit scripts and original GTF.
- `Hmel_ribosomal_genes.txt` and `Hmel_mito_genenames.txt`.

## Outputs

- Reconstructable `Hmel2.5_PLUS_MITO_GENOME` reference and STAR index.

## Software

STAR and the original annotation editing tools; versions and index-generation arguments still needed.

## Execution and interpretation

Recover FASTA/GTF provenance, ivory feature edits, and index-building commands before generating a reference. The supplied alignment script points to an existing index; it does not build one.

See [intake review](../docs/INTAKE_REVIEW.md) for unresolved decisions and [file requests](../docs/FILE_REQUESTS.md) for the complete acquisition checklist.

## Code

No reference-building code has been supplied. See the exact reference filenames in [the alignment workflow](02_starsolo_counts.md).

[Return to analysis index](../README.md#analysis-index).
