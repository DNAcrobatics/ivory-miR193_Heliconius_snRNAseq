# Gene mapping, GO enrichment and plotting

**Paper:** Figure 7C/F.

**Status:** cleaned scripts supplied; final inputs/settings and result matching pending.

Recover gene-to-FlyBase mapping, submitted gene lists and enrichment settings.

## Inputs

- Final DE tables, `HomologGeneNames.txt`, `fbgn_fbtr_fbpp_expanded_fb_2024_01.tsv`.
- UpGenes/DownGenes and late equivalents, including FBpp submission lists.
- Four DirectGoEnrichment_PANGEA CSV exports named in the plotting code.

## Outputs

- GO enrichment result tables and early/late PDF plots.

## Software

R and ggplot2 for plotting; enrichment service/database version and submission settings need recovery.

## Execution and interpretation

These blocks plot previously downloaded results; they do not execute the enrichment test. Recover the universe, submitted identifiers, term library, correction and export settings. Early filtering mixes Benjamini–Yekutieli and Benjamini–Hochberg, while later filters use BY. Reconcile this with the paper’s permutation/FDR description before choosing final code.

See [intake review](../docs/INTAKE_REVIEW.md) for unresolved decisions and [file requests](../docs/FILE_REQUESTS.md) for the complete acquisition checklist.

## Cleaned code and execution

- [go_enrichment.R](../scripts/plots/go_enrichment.R)

The script plots the exported enrichment results. Correction filter/color columns are explicit, so differing historical BY/BH choices are recorded rather than guessed. The underlying enrichment run is still needed.

See [run commands and input formats](../scripts/README.md), [cleanup notes](../docs/CLEANUP_NOTES.md) and the [original source archive](../archive/README.md).

[Return to analysis index](../README.md#analysis-index).
