# Early and late differential expression

**Paper:** Figure 6D–G and Figure 7A/B/D/E.

**Status:** cleaned scripts supplied; final inputs/settings and result matching pending.

Recover four main contrasts and their volcano/dot-plot source tables.

## Inputs

- Final early/late subsets with `Iv` and `mut` labels and confirmed assay.
- Gene labels: `Hi_wt_v_mut_for_dotplot.txt`, `WtVMutInHiLate.Markers_genenames.txt`, `HiVsLoLate.markers_naems.txt`.
- Saved final DE tables and figure PDFs.

## Outputs

- Early/late ivory-high versus low and WT versus mutant in ivory-high nuclei; volcano and dot plots.

## Software

R, Seurat, ggplot2 and plotting packages; exact environment pending.

## Execution and interpretation

The recovered main calls use Seurat FindMarkers on nucleus subsets and omit test.use. The current figure legends specify Wilcoxon/Bonferroni, while other manuscript text describes pseudobulk; confirm the final analysis before preparing release code. ROC exploratory calls elsewhere are separate. In the early WT/mutant branch the computed object is WtVMutInHi.Markers2, but plotting reads WtVMutInHi.Markers. Do not infer that these are interchangeable. Confirm both identities, assay, covariates, filtering, and biological replicate handling for each contrast.

See [intake review](../docs/INTAKE_REVIEW.md) for unresolved decisions and [file requests](../docs/FILE_REQUESTS.md) for the complete acquisition checklist.

## Cleaned code and execution

- [05_find_markers.R](../scripts/snrna/05_find_markers.R)
- [volcano.R](../scripts/plots/volcano.R)
- [dotplot.R](../scripts/plots/dotplot.R)

Use frozen membership tables for all four comparisons, with explicit group direction and assay. Wilcoxon is explicit in this cleaned nucleus-level branch; the pseudobulk wording discrepancy still needs resolution. These plots are clean equivalents awaiting final labels/palette/output matching.

See [run commands and input formats](../scripts/README.md), [cleanup notes](../docs/CLEANUP_NOTES.md) and the [original source archive](../archive/README.md).

[Return to analysis index](../README.md#analysis-index).
