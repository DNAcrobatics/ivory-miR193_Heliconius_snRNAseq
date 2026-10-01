# Aggregate expression and ivory ranking

**Paper:** Supplementary Figure S8.

**Status:** Session fragments.

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

## Code

### All-nucleus aggregate rankings

**Source:** [ivory_reanalysis.expanded.pasted.r, lines 2453–2470](../archive/original/ivory_proj_transfer_aug2025/ivory_reanalysis.expanded.pasted.r#L2453-L2470).

**Archive treatment:** exact excerpt; no syntax repairs, parameter changes, or inferred session setup. This is a recovered session fragment, not a standalone script.

```r
pb <- AggregateExpression(Hm.sops.sct,
                          group.by = "celltypes",
                          slot = "counts")
pb <-  edgeR::cpm(pb$RNA)
rankedSc <- sort(pb[, "scale"], decreasing=TRUE)
rankedSo <- sort(pb[,"socket"], decreasing=TRUE)
rankedSp <- sort(pb[,"sop"], decreasing=TRUE)

rankedSc <- rankedSc[rankedSc > 1]
rankedSo <- rankedSo[rankedSo > 1]
rankedSp <- rankedSp[rankedSp > 1]



which(names(rankedSc)=='Ivory')
which(names(rankedSo)=='Ivory')
which(names(rankedSp)=='Ivory')
```

### Mutant-excluded histogram branch

**Source:** [ivory_reanalysis.expanded.pasted.r, lines 2473–2493](../archive/original/ivory_proj_transfer_aug2025/ivory_reanalysis.expanded.pasted.r#L2473-L2493).

**Archive treatment:** exact excerpt; no syntax repairs, parameter changes, or inferred session setup. This is a recovered session fragment, not a standalone script.

```r
pb <- AggregateExpression(subset(Hm.sops.sct, mut!='mut'),
                          group.by = "celltypes",
                          slot = "counts")
pb <-  edgeR::cpm(pb$RNA)
rankedSc <- sort(pb[, "scale"], decreasing=TRUE)
rankedSo <- sort(pb[,"socket"], decreasing=TRUE)
rankedSp <- sort(pb[,"sop"], decreasing=TRUE)

rankedSc <- rankedSc[rankedSc > 1]
rankedSo <- rankedSo[rankedSo > 1]
rankedSp <- rankedSp[rankedSp > 1]


ggplot(as.data.frame(rankedSc),aes((rankedSc)))+
    geom_histogram()+
    scale_x_log10()+
    geom_vline(xintercept=rankedSc[names(rankedSc)=='Ivory'])+
    labs(x='gene counts in scale-building cells')+
    theme_classic()
dev.copy(pdf, 'Fig_Sx_gene_counts_in_scale_cells_histogram.pdf', width=5,height=2)
dev.off()
```


[Return to analysis index](../README.md#analysis-index).
