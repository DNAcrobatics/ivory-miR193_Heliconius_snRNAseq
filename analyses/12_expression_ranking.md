# Aggregate expression and ivory ranking

**Paper:** Supplementary Figure S8.

**Status:** cleaned scripts supplied; final inputs/settings and result matching pending.

Rank aggregate expression across lineage cell types and recover the ivory histogram.

## Inputs

- Final `Hm.sops.sct`, cell-type labels and WT/mutant metadata.
- Final Figure S8 PDF and underlying values.

## Outputs

- Cell-type aggregate CPM rankings and `Fig_Sx_gene_counts_in_scale_cells_histogram.pdf` in the recovered branch.

## Software

R, Seurat, edgeR and ggplot2; exact package environment pending.

## Execution and interpretation

The paste includes both all-nucleus and mutant-excluded aggregation branches. Confirm the subset used for Figure S8 and the final ivory/heph annotations. AggregateExpression plus CPM here supports expression ranking; it is not evidence of a replicate-level DE model.

See [intake review](../docs/INTAKE_REVIEW.md) for unresolved decisions and [file requests](../docs/FILE_REQUESTS.md) for the complete acquisition checklist.

## Cleaned code and execution

- [06_expression_ranking.R](../scripts/snrna/06_expression_ranking.R)

Supply the intended cohort explicitly. The cleaned aggregate CPM ranking is independent of DE testing and writes both ranked values and a vector histogram.

See [run commands and input formats](../scripts/README.md), [cleanup notes](../docs/CLEANUP_NOTES.md) and the [original source archive](../archive/README.md).

[Return to analysis index](../README.md#analysis-index).
