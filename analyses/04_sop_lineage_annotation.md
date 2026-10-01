# SOP lineage integration and annotation

**Paper:** Figure 5A–D; Supplementary Figure S7.

**Status:** Session fragments.

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

## Code

### Recovered lineage integration branch

**Source:** [ivory_reanalysis.expanded.pasted.r, lines 145–164](../archive/original/ivory_proj_transfer_aug2025/ivory_reanalysis.expanded.pasted.r#L145-L164).

**Archive treatment:** exact excerpt; no syntax repairs, parameter changes, or inferred session setup. This is a recovered session fragment, not a standalone script.

```r
Hm_list <- SplitObject(Hm_sops, split.by = "race")
#######Hm_list <- lapply(X = Hm_list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = Hm_list, nfeatures = 7000)
Hm_list <- PrepSCTIntegration(object.list = Hm_list, anchor.features = features)

Hm.anchors <- FindIntegrationAnchors(object.list = Hm_list, normalization.method = "SCT",
    anchor.features = features)
Hm.sops.sct <- IntegrateData(anchorset = Hm.anchors, normalization.method = "SCT")

##SaveH5Seurat(Hm.sops.sct, 'Hm.sops.sct', overwrite = TRUE, verbose = TRUE)
##Hm.sops.sct <- LoadH5Seurat('Hm.sops.sct.h5seurat')

### so far, 1:3 works best
#Hm.sops.sct <- subset(Hm.sops.sct, seurat_clusters %ni% c(16,19))
DefaultAssay(Hm.sops.sct) <- 'integrated'
Hm.sops.sct <- RunPCA(Hm.sops.sct,  npcs = 50, verbose = F)
Hm.sops.sct <- RunUMAP(Hm.sops.sct,dims = c(1:3))
Hm.sops.sct <- FindNeighbors(Hm.sops.sct, dims =  c(1:3)
Hm.sops.sct <- FindClusters(Hm.sops.sct, resolution =1)
#Hm.sops.sct <- NormalizeData(Hm.sops.sct, normalization.method = "LogNormalize", assay = 'RNA')
```


[Return to analysis index](../README.md#analysis-index).
