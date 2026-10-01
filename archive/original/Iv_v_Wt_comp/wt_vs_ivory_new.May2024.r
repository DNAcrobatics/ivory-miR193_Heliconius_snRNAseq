library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratDisk)
library(sctransform)

## read data
iv_20.data <- Read10X(data.dir="iv_20.filtered")
iv_33.data <- Read10X(data.dir="iv_33.filtered")
iv_38.data <- Read10X(data.dir="iv_38.filtered")
iv_55.data <- Read10X(data.dir="iv_55.filtered")
iv_61.data <- Read10X(data.dir="iv_61.filtered")
Hmm_40.data <- Read10X(data.dir='filtered.Hmm_40')
Hmr_35.data <- Read10X(data.dir='filtered.Hmr_35')
Hmr_40.data <- Read10X(data.dir='filtered.Hmr_40')
Hmr_45.data <- Read10X(data.dir='filtered.Hmr_45')
Hmr_30.data <- Read10X(data.dir='filtered.Hmr_30')
##create objects
iv_20 <- CreateSeuratObject(counts = iv_20.data, project = "iv_20", min.cells = 3, min.features = 2000)
iv_33 <- CreateSeuratObject(counts = iv_33.data, project = "iv_33", min.cells = 3, min.features = 2000)
iv_38 <- CreateSeuratObject(counts = iv_38.data, project = "iv_38", min.cells = 3, min.features = 2000)
iv_55 <- CreateSeuratObject(counts = iv_55.data, project = "iv_55", min.cells = 3, min.features = 2000)
iv_61 <- CreateSeuratObject(counts = iv_61.data, project = "iv_61", min.cells = 3, min.features = 2000)
Hmm_40 <- CreateSeuratObject(counts= Hmm_40.data, project= 'Hmm_40',min.cells= 3, min.features= 2000)
Hmr_35 <- CreateSeuratObject(counts= Hmr_35.data, project= 'Hmr_35',min.cells= 3, min.features= 2000)
Hmr_40 <- CreateSeuratObject(counts= Hmr_40.data, project= 'Hmr_40',min.cells= 3, min.features= 2000)
Hmr_45 <- CreateSeuratObject(counts= Hmr_45.data, project= 'Hmr_45',min.cells= 3, min.features= 2000)
Hmr_30 <- CreateSeuratObject(counts= Hmr_30.data, project= 'Hmr_30',min.cells= 3, min.features= 2000)

iv_20$time <- 32
iv_33$time <- 35
iv_38$time <- 58
iv_55$time <- 46
iv_61$time <- 36
Hmm_40$time <- 40
Hmr_35$time <- 35
Hmr_40$time <- 40
Hmr_45$time <- 45
Hmr_30$time <- 30

iv_20$zygosity <- 'wt'
iv_33$zygosity <- 'het'
iv_38$zygosity <- 'mut'
iv_55$zygosity <- 'mut'
iv_61$zygosity <- 'mut'
Hmm_40$zygosity <- 'wt'
Hmr_35$zygosity <- 'wt'
Hmr_40$zygosity <- 'wt'
Hmr_45$zygosity <- 'wt'
Hmr_30$zygosity <- 'wt'


print(2)

#reduced time series
Iv.merged <- merge(Hmr_30, y = c(iv_20,iv_33,iv_38,iv_55,iv_61,Hmm_40,Hmr_35,Hmr_40,Hmr_45),
                  add.cell.ids = c('Hmr_30','iv_20','iv_33','iv_38','iv_55','iv_61','Hmm_40','Hmr_35','Hmr_40','Hmr_45'), project = "Iv_merged_series")

unique(sapply(X = strsplit(colnames(Iv.merged), split = "_"), FUN = "[", 1))
table(Iv.merged$orig.ident)
Iv.list <- SplitObject(Iv.merged, split.by = "orig.ident")


#print(3)

## make list of ribosomal genes
GeneNames <- Reduce(intersect, lapply(Iv.list, rownames))
ribo_genes <- read.table('../Iv_v_Wt_comp/Hmel_ribosomal_genes.txt',h=F)
ribo_genes$V2 <- NULL
ribo_genes <- GeneNames[GeneNames %in% as.vector(ribo_genes$V1)]

## make list of mitochondrial genes
mito_genes <- read.table('../Hmel2.5_PLUS_MITO_GENOME/Hmel_mito_genenames.txt', h=F)
mito_genes <- GeneNames[GeneNames %in% as.vector(mito_genes$V1)]

### Flagging mitochondrial and ribosomal genes
Iv.merged  <- PercentageFeatureSet(Iv.merged, features=as.vector(ribo_genes), col.name = "percent.Rb")

Iv.merged  <- PercentageFeatureSet(Iv.merged, features=as.vector(mito_genes), col.name = "percent.MT")

##Vlnplots
#png(file="hm_ts.violin_plot_nFeat.nCount.pcMTs.png")
#VlnPlot(Iv.merged, features = c("nFeature_RNA", "nCount_RNA","percent.Rb","percent.MT"), group.by = 'time',  ncol = 4, pt.size = 0)
#dev.off()

Iv.merged <- subset(Iv.merged, subset = nFeature_RNA > 200 & percent.MT < 20 )
Iv.merged <- SCTransform(Iv.merged, verbose = T, variable.features.n = 5000,vars.to.regress = 'time')
Iv.merged <- RunPCA(Iv.merged,  npcs = 50, verbose = T)
Iv.merged <- RunUMAP(Iv.merged,  dims = 1:8)
Iv.merged <- FindNeighbors(Iv.merged, dims = 1:8)
Iv.merged <- FindClusters(Iv.merged, resolution = 1)

#png(file="DimPlot_5000_8")
# print(DimPlot(Iv.merged, group.by='time'))
# dev.off()
Iv.list <- SplitObject(Iv.merged, split.by = "zygosity")
#######Hm_list <- lapply(X = Hm_list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = Iv.list, nfeatures = 5000)
Iv.list <- PrepSCTIntegration(object.list = Iv.list, anchor.features = features)

Iv.anchors <- FindIntegrationAnchors(object.list = Iv.list, normalization.method = "SCT",
    anchor.features = features)
Iv.combined.sct <- IntegrateData(anchorset = Iv.anchors, normalization.method = "SCT")

Iv.combined.sct <- RunPCA(Iv.combined.sct,  npcs = 50, verbose = T)
Iv.combined.sct <- RunUMAP(Iv.combined.sct,  dims = 1:8)
Iv.combined.sct <- FindNeighbors(Iv.combined.sct, dims = 1:8)
Iv.combined.sct <- FindClusters(Iv.combined.sct, resolution = 1)

SaveH5Seurat(Iv.combined.sct, 'Iv_combined.full.sct', overwrite = TRUE, verbose = TRUE)
