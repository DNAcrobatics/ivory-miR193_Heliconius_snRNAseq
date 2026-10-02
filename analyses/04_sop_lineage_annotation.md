# SOP lineage integration and annotation

**Paper:** Figure 5A–D; Supplementary Figure S7.

**Status:** cleaned scripts supplied; final inputs/settings and result matching pending.

Recover the SOP/scale/socket lineage objects, annotations and UMAP panels.

## Inputs

- Final full integrated Seurat object.
- SOP selection, cell-type labels, race and developmental metadata; palettes and helper functions.

## Outputs

- `Hm.sops.sct`, lineage UMAP/feature panels and race/cell-type proportions.

## Software

R, Seurat, SeuratDisk, ggplot2 and plotting helpers from the original session.

## Execution and interpretation

The recovered code contains several alternative UMAP/integration trials. The 7,000-feature, three-PC branch below is historical, contains an incomplete FindNeighbors call, and has not been selected as the final paper workflow. Recover the exact cluster-to-cell-type labels and final selection. Figure 5 proportions must use confirmed denominators and exclusions.

See [intake review](../docs/INTAKE_REVIEW.md) for unresolved decisions and [file requests](../docs/FILE_REQUESTS.md) for the complete acquisition checklist.

## Cleaned code and execution

- [02_integrate.R](../scripts/snrna/02_integrate.R)
- [lineage_panels.R](../scripts/plots/lineage_panels.R)

Supply the final SOP barcode selection via --cells and the confirmed integration settings. A complete integration step repairs the missing FindNeighbors parenthesis. Final cell-type labels must already be in the input metadata; they are not inferred from exploratory cluster trials.

See [run commands and input formats](../scripts/README.md), [cleanup notes](../docs/CLEANUP_NOTES.md) and the [original source archive](../archive/README.md).

[Return to analysis index](../README.md#analysis-index).
