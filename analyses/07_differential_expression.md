# Early and late differential expression

**Paper:** Figure 6D–G and Figure 7A/B/D/E.

**Status:** Session fragments; final tests unresolved.

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

## Code

### Early WT versus mutant call

**Source:** [ivory_reanalysis.expanded.pasted.r, lines 1044–1045](../archive/original/ivory_proj_transfer_aug2025/ivory_reanalysis.expanded.pasted.r#L1044-L1045).

**Archive treatment:** exact excerpt; no syntax repairs, parameter changes, or inferred session setup. This is a recovered session fragment, not a standalone script.

```r
Idents(test) <- test$mut
WtVMutInHi.Markers2 <- FindMarkers(subset(test,Iv =='hi'), ident.1='wt',recorrect_umi=FALSE)
```

### Early ivory-high versus low

**Source:** [ivory_reanalysis.expanded.pasted.r, lines 1124–1126](../archive/original/ivory_proj_transfer_aug2025/ivory_reanalysis.expanded.pasted.r#L1124-L1126).

**Archive treatment:** exact excerpt; no syntax repairs, parameter changes, or inferred session setup. This is a recovered session fragment, not a standalone script.

```r
Idents(test) <- test$Iv
DefaultAssay(test) <- "RNA"
HiVsLo.markers <- FindMarkers(test, ident.1="hi", ident.2="lo",recorrect_umi=FALSE)
```

### Late WT versus mutant call

**Source:** [ivory_reanalysis.expanded.pasted.r, lines 1415–1416](../archive/original/ivory_proj_transfer_aug2025/ivory_reanalysis.expanded.pasted.r#L1415-L1416).

**Archive treatment:** exact excerpt; no syntax repairs, parameter changes, or inferred session setup. This is a recovered session fragment, not a standalone script.

```r
Idents(test1) <- test1$mut
WtVMutInHiLate.Markers <- FindMarkers(subset(test1,Iv=='hi'), ident.1='wt',recorrect_umi=FALSE)
```

### Late ivory-high versus low

**Source:** [ivory_reanalysis.expanded.pasted.r, lines 1510–1512](../archive/original/ivory_proj_transfer_aug2025/ivory_reanalysis.expanded.pasted.r#L1510-L1512).

**Archive treatment:** exact excerpt; no syntax repairs, parameter changes, or inferred session setup. This is a recovered session fragment, not a standalone script.

```r
Idents(test1) <- test1$Iv
DefaultAssay(test1) <- "RNA"
HiVsLoLate.markers <- FindMarkers(test1, ident.1="hi", ident.2="lo",recorrect_umi=FALSE)
```


[Return to analysis index](../README.md#analysis-index).
