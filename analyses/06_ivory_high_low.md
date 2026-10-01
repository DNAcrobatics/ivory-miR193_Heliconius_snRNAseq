# Define ivory-high and ivory-low scale nuclei

**Paper:** Figures 5E/F and 6–7.

**Status:** Alternative branches; final choice needed.

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

## Code

### Earlier late-sample classification

**Source:** [ivory_reanalysis.expanded.pasted.r, lines 1400–1412](../archive/original/ivory_proj_transfer_aug2025/ivory_reanalysis.expanded.pasted.r#L1400-L1412).

**Archive treatment:** exact excerpt; no syntax repairs, parameter changes, or inferred session setup. This is a recovered session fragment, not a standalone script.

```r
test1 <- subset(Hm.sops.sct, orig.ident %in% c('iv_55','Hmm_40','Hmr_40') & celltypes=='scale')
DefaultAssay(test1) <- 'RNA'
barcodes <- subset(x = test1, Ivory > 1 & race=='iv')
barcodes2<- subset(x = test1, Ivory > 4 & race!='iv')
test1$Iv <- ""
test1$Iv[WhichCells(test1) %in%  WhichCells(barcodes)] <- "hi"
test1$Iv[WhichCells(test1) %in%  WhichCells(barcodes2)] <- 'hi'
test1$Iv[test1$Iv == ''] <- 'lo'
#test$Iv[WhichCells(test) %ni%  WhichCells(barcodes)] <- "lo"
#test$Iv[WhichCells(test) %ni%  WhichCells(barcodes2)] <- 'lo'

Idents(test1) <- test1$Iv
VlnPlot(test1,'Ivory',split.by='race')
```

### Later sample normalization branch

**Source:** [ivory_reanalysis.expanded.pasted.r, lines 2352–2385](../archive/original/ivory_proj_transfer_aug2025/ivory_reanalysis.expanded.pasted.r#L2352-L2385).

**Archive treatment:** exact excerpt; no syntax repairs, parameter changes, or inferred session setup. This is a recovered session fragment, not a standalone script.

```r

#### creating a new layer where i normalise to ivory in each sample. 

# PARAMETERS
ref_gene <- "Ivory"
sample_column <- "orig.ident"
orig_assay <- "RNA"
new_assay <- "RNA_by_Ivory"

#test <- Hm.sops.sct
# 1. linearize log1p data (if needed)
logdata <- GetAssayData(Hm.sops.sct, assay = orig_assay, slot = "data")
lin <- expm1(logdata)

# 2. compute per-sample mean of ref_gene and attach
cell2sample <- Hm.sops.sct@meta.data[[sample_column]]
per_sample_mean <- tapply(as.numeric(lin[ref_gene, , drop = TRUE]), cell2sample, mean, na.rm = TRUE)
Hm.sops.sct@meta.data$norm_factor <- per_sample_mean[cell2sample]

# 3. divide and re-log, store in new assay
lin_div <- sweep(lin, 2, Hm.sops.sct@meta.data$norm_factor, FUN = "/")
new_logdata <- log1p(lin_div)
Hm.sops.sct[[new_assay]] <- CreateAssayObject(data = new_logdata)

DefaultAssay(Hm.sops.sct) <- "RNA_by_Ivory"

test <- subset(test, orig.ident %in% c('iv_20','iv_61','Hmm_25','Hmr_25','Hmr_30') & celltypes=='scale')
barcodes <- subset(x = test, Ivory > 0.5 &)
test$Iv <- ""
test$Iv[WhichCells(test) %in%  WhichCells(barcodes)] <- "hi"
test$Iv[test$Iv == ''] <- 'lo'

Idents(test) <- test$Iv
HiVsLo.markers <- FindMarkers(test, ident.1="hi", ident.2="lo",recorrect_umi=FALSE)
```


[Return to analysis index](../README.md#analysis-index).
