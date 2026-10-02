# Define ivory-high and ivory-low scale nuclei

**Paper:** Figures 5E/F and 6–7.

**Status:** cleaned scripts supplied; final inputs/settings and result matching pending.

Make the classification used by the early and late expression comparisons explicit.

## Inputs

- Confirmed early/late scale-nucleus subsets, race metadata and RNA assay.
- Final per-cell classification table or final threshold code.

## Outputs

- Per-cell `Iv` high/low labels and assay/threshold provenance.

## Software

R and Seurat; package version and active assay must be recorded.

## Execution and interpretation

An earlier branch uses RNA Ivory >1 for race iv and >4 for other races. A later branch creates RNA_by_Ivory by dividing linearized normalized RNA expression by sample mean Ivory, followed by an incomplete >0.5 classification. These alternatives are preserved; neither is silently declared final. The normalization also needs a check for zero sample means before it can become production code.

See [intake review](../docs/INTAKE_REVIEW.md) for unresolved decisions and [file requests](../docs/FILE_REQUESTS.md) for the complete acquisition checklist.

## Cleaned code and execution

- [04_classify_ivory.R](../scripts/snrna/04_classify_ivory.R)

The complete RNA branch is implemented with explicitly supplied thresholds. The incomplete RNA_by_Ivory experiment is retained in the archive and is not silently repaired into a scientific choice.

See [run commands and input formats](../scripts/README.md), [cleanup notes](../docs/CLEANUP_NOTES.md) and the [original source archive](../archive/README.md).

[Return to analysis index](../README.md#analysis-index).
