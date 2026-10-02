# Cell-cycle scoring and phase plots

**Paper:** Figure 5E/F.

**Status:** cleaned scripts supplied; final inputs/settings and result matching pending.

Score S and G2/M phases with the supplied cross-species gene mapping.

## Inputs

- Annotated `Hm.sops.sct` and early/late objects `test` and `test1`.
- [CellCycleGenes.txt](../archive/original/ivory_proj_transfer_aug2025/CellCycleGenes.txt).

## Outputs

- S.Score, G2M.Score and Phase metadata; early/late phase proportion plots.

## Software

R, Seurat and ggplot2, with the scoring assay and package version still to confirm.

## Execution and interpretation

The mapping contains 66 rows: 29 S and 37 G2/M. Retain its original NO_HIT entries and duplicate mappings in the archive. Confirm how the final scoring run handled absent and repeated features; the snippet passes the HmID column directly. The phase plot fragments occur before scoring in the pasted session, so their dependencies must be restored before making a runnable workflow.

See [intake review](../docs/INTAKE_REVIEW.md) for unresolved decisions and [file requests](../docs/FILE_REQUESTS.md) for the complete acquisition checklist.

## Cleaned code and execution

- [03_cell_cycle.R](../scripts/snrna/03_cell_cycle.R)
- [lineage_panels.R](../scripts/plots/lineage_panels.R)

The cleaned scorer writes an effective-marker audit and explicitly filters absent features and repeated markers within phase. Verify that adaptation against historical scoring. The panel script can plot Phase proportions on a frozen scale-nucleus subset.

See [run commands and input formats](../scripts/README.md), [cleanup notes](../docs/CLEANUP_NOTES.md) and the [original source archive](../archive/README.md).

[Return to analysis index](../README.md#analysis-index).
