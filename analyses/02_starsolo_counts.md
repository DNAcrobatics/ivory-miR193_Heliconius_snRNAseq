# Single-nucleus alignment and counting

**Paper:** Upstream of Figures 5–7 and S7–S8.

**Status:** cleaned scripts supplied; final inputs/settings and result matching pending.

Align snRNA-seq reads and generate STARsolo gene-full count matrices.

## Inputs

- Paired gzipped FASTQs; historical script passes R2 first and R1 second.
- The ordered `jobs` sample list; it has not been supplied.
- The reference from workflow 01.

## Outputs

- STARsolo GeneFull matrices and alignment logs under each `redo.<sample>` prefix.

## Software

STAR, samtools, micromamba and Slurm. Reported versions need confirmation against run logs; see [environments](../envs/README.md).

## Execution and interpretation

The source is a five-task Slurm array. Confirm which run supplied the final filtered matrices and how `redo.*` outputs became `iv_*.filtered`/`filtered.Hm*`. Preserve CB length 16, UMI length 12, no whitelist, EM multimappers, and EmptyDrops_CR unless the final run records say otherwise. BAM and barcode attribute options are commented out, so downstream per-barcode coverage cannot yet be reproduced from this script alone.

See [intake review](../docs/INTAKE_REVIEW.md) for unresolved decisions and [file requests](../docs/FILE_REQUESTS.md) for the complete acquisition checklist.

## Cleaned code and execution

- [starsolo.sh](../scripts/alignment/starsolo.sh)

Counting parameters and R2/R1 order are preserved; scheduler settings, paths, threads and sample selection are explicit.

See [run commands and input formats](../scripts/README.md), [cleanup notes](../docs/CLEANUP_NOTES.md) and the [original source archive](../archive/README.md).

[Return to analysis index](../README.md#analysis-index).
