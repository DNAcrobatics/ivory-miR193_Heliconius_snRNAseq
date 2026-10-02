library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratDisk)
library(sctransform)

#
#Hmr_5i.data <- Read10X(data.dir='filtered.Hmr_5i')
#Hmm_10.data <- Read10X(data.dir="filtered.Hmm_10")
#Hmm_25.data <- Read10X(data.dir="filtered.Hmm_25")
#Hmm_40.data <- Read10X(data.dir="filtered.Hmm_40")
#Hmr_10.data <- Read10X(data.dir="filtered.Hmr_10")
#Hmr_15.data <- Read10X(data.dir="filtered.Hmr_15")
#Hmr_20.data <- Read10X(data.dir="filtered.Hmr_20")
#Hmr_25.data <- Read10X(data.dir="filtered.Hmr_25")
#Hmr_30.data <- Read10X(data.dir="filtered.Hmr_30")
#Hmr_35.data <- Read10X(data.dir="filtered.Hmr_35")
#Hmr_40.data <- Read10X(data.dir="filtered.Hmr_40")
#Hmr_45.data <- Read10X(data.dir="filtered.Hmr_45")
#Hmr_50.data <- Read10X(data.dir="filtered.Hmr_50")
#iv_20.data <- Read10X(data.dir="redo.iv_20")
#iv_33.data <- Read10X(data.dir="redo.iv_33")
#iv_38.data <- Read10X(data.dir="redo.iv_38")
#iv_55.data <- Read10X(data.dir="redo.iv_55")
#iv_61.data <- Read10X(data.dir="redo.iv_61")


print(1)

#Hmr_5i <- CreateSeuratObject(counts = Hmr_5i.data, project = "Hmr_5i", min.cells = 3, min.features = 2000) 
#Hmm_10 <- CreateSeuratObject(counts = Hmm_10.data, project = "Hmm_10", min.cells = 3, min.features = 2000)
#Hmm_25 <- CreateSeuratObject(counts = Hmm_25.data, project = "Hmm_25", min.cells = 3, min.features = 2000)
#Hmm_40 <- CreateSeuratObject(counts = Hmm_40.data, project = "Hmm_40", min.cells = 3, min.features = 2000)
#Hmr_10 <- CreateSeuratObject(counts = Hmr_10.data, project = "Hmr_10", min.cells = 3, min.features = 2000)
#Hmr_15 <- CreateSeuratObject(counts = Hmr_15.data, project = "Hmr_15", min.cells = 3, min.features = 2000)
#Hmr_20 <- CreateSeuratObject(counts = Hmr_20.data, project = "Hmr_20", min.cells = 3, min.features = 2000)
#Hmr_25 <- CreateSeuratObject(counts = Hmr_25.data, project = "Hmr_25", min.cells = 3, min.features = 2000)
#Hmr_30 <- CreateSeuratObject(counts = Hmr_30.data, project = "Hmr_30", min.cells = 3, min.features = 2000)
#iv_20 <- CreateSeuratObject(counts = iv_20.data, project = "iv_20", min.cells = 3, min.features = 2000)
#iv_33 <- CreateSeuratObject(counts = iv_33.data, project = "iv_33", min.cells = 3, min.features = 2000)
#iv_38 <- CreateSeuratObject(counts = iv_38.data, project = "iv_38", min.cells = 3, min.features = 2000)
#iv_55 <- CreateSeuratObject(counts = iv_55.data, project = "iv_55", min.cells = 3, min.features = 2000)
#iv_61 <- CreateSeuratObject(counts = iv_61.data, project = "iv_61", min.cells = 3, min.features = 2000)
#Hmr_35 <- CreateSeuratObject(counts = Hmr_35.data, project = "Hmr_35", min.cells = 3, min.features = 2000)
#Hmr_40 <- CreateSeuratObject(counts = Hmr_40.data, project = "Hmr_40", min.cells = 3, min.features = 2000)
#Hmr_45 <- CreateSeuratObject(counts = Hmr_45.data, project = "Hmr_45", min.cells = 3, min.features = 2000)
#Hmr_50 <- CreateSeuratObject(counts = Hmr_50.data, project = "Hmr_50", min.cells = 3, min.features = 2000)

#Hmr_5i$time <- 0
#Hmm_10$time <- 10
#Hmm_25$time <- 25
#Hmm_40$time <- 40
#Hmr_10$time <- 10
#Hmr_15$time <- 15
#Hmr_20$time <- 20
#Hmr_25$time <- 25
#Hmr_30$time <- 30
#Hmr_35$time <- 35
#Hmr_40$time <- 40
#Hmr_45$time <- 45
#Hmr_50$time <- 50 
#iv_20$time <- 32
#iv_33$time <- 35
#iv_38$time <- 56
#iv_55$time <- 46
#iv_61$time <- 36

print(2)

#reduced time series
#Hm_merged <- merge(Hmr_10, y = c(Hmm_10,Hmm_25,Hmr_15,Hmr_20,Hmr_25,Hmr_30),
#                  add.cell.ids = c('Hmr_10','Hmm_10','Hmm_25','Hmr_15','Hmr_20','Hmr_25','Hmr_30'), project = "Hm_merged_series")

#unique(sapply(X = strsplit(colnames(Hm_merged), split = "_"), FUN = "[", 1))
#table(Hm_merged$orig.ident)
#Hm.list <- SplitObject(Hm_merged, split.by = "orig.ident")

#full time series , removing hmr_5i
#Hm_merged <- merge(Hmm_10, y = c(Hmm_25,Hmm_40,Hmr_10,Hmr_15,Hmr_20,Hmr_25,Hmr_30,Hmr_35,Hmr_40,Hmr_45,Hmr_50,iv_20,iv_33,iv_38,iv_55,iv_61),
#                  add.cell.ids = c('Hmm_10','Hmm_25','Hmm_40','Hmr_10','Hmr_15','Hmr_20','Hmr_25','Hmr_30','Hmr_35','Hmr_40','Hmr_45','Hmr_50','iv_20','iv_33','iv_38','iv_55','iv_61'), project = "Hm_merged_series")

#unique(sapply(X = strsplit(colnames(Hm_merged), split = "_"), FUN = "[", 1))
#table(Hm_merged$orig.ident)
#Hm.list <- SplitObject(Hm_merged, split.by = "orig.ident")


print(3)

## make list of ribosomal genes
#GeneNames <- Reduce(intersect, lapply(Hm.list, rownames))
#ribo_genes <- read.table('../Iv_v_Wt_comp/Hmel_ribosomal_genes.txt',h=F)
#ribo_genes$V2 <- NULL
#ribo_genes <- GeneNames[GeneNames %in% as.vector(ribo_genes$V1)]

## make list of mitochondrial genes
#mito_genes <- read.table('../Hmel2.5_PLUS_MITO_GENOME/Hmel_mito_genenames.txt', h=F)
#mito_genes <- GeneNames[GeneNames %in% as.vector(mito_genes$V1)]

### Flagging mitochondrial and ribosomal genes
#Hm_merged  <- PercentageFeatureSet(Hm_merged, features=as.vector(ribo_genes), col.name = "percent.Rb")

#Hm_merged  <- PercentageFeatureSet(Hm_merged, features=as.vector(mito_genes), col.name = "percent.MT")

##Vlnplots
#png(file="hm_ts.violin_plot_nFeat.nCount.pcMTs.png")
#VlnPlot(Hm_merged, features = c("nFeature_RNA", "nCount_RNA","percent.Rb","percent.MT"), group.by = 'time',  ncol = 4, pt.size = 0)
#dev.off() 

print(3)
########
# i will run with some filtering on mitochondria, that should be enough to clean things right up
##Subsetting the data for all samples
#Hm_merged <- subset(Hm_merged, subset = nFeature_RNA > 200 & percent.MT < 20 )
#dim(Hm_merged)
#dimensions <- c(3,4,5,6,7,8)
#vfeatures <- c(3000,4000,5000,6000)
#for (x in vfeatures) {
#   Hm_merged <- SCTransform(Hm_merged, verbose = T, variable.features.n = x, vars.to.regress = 'time')
#   Hm_merged <- RunPCA(Hm_merged,  npcs = 50, verbose = T)
#   for (y in dimensions) {
#     print(x)
#     print(y)
#     Hm_merged <- RunUMAP(Hm_merged,  dims = 1:y)
#     Hm_merged <- FindNeighbors(Hm_merged, dims = 1:y)
#     Hm_merged <- FindClusters(Hm_merged, resolution = 1)
#     png(file=paste("dimplot_", x, "_vf_", y , "_noTimeReg_noLarva.png", sep = ""))
#      print(DimPlot(Hm_merged, group.by='time', label=T))
#      dev.off()
#     png(file=paste("dimplot_", x, "vf", y, "ivory.png", sep = ""))
#      print(FeaturePlot(Hm_merged, "Ivory"))
#      dev.off()
#   }
#} 

#Hm_merged <- subset(Hm_merged, subset = nFeature_RNA > 200 & percent.MT < 20 )
#Hm_merged <- SCTransform(Hm_merged, verbose = T, variable.features.n = 6000,vars.to.regress = 'time')
#Hm_merged <- RunPCA(Hm_merged,  npcs = 50, verbose = T)
#Hm_merged <- RunUMAP(Hm_merged,  dims = 1:4)
#Hm_merged <- FindNeighbors(Hm_merged, dims = 1:4)
#Hm_merged <- FindClusters(Hm_merged, resolution = 1)

Hm_merged$race <- ifelse(Hm_merged$orig.ident == "Hmm_10" | Hm_merged$orig.ident == "Hmm_25" | Hm_merged$orig.ident == 'Hmm_40', "mal", ifelse(Hm_merged$orig.ident == 'iv_20' | Hm_merged$orig.ident == 'iv_33' | Hm_merged$orig.ident == 'iv_38' | Hm_merged$orig.ident == 'iv_55' | Hm_merged$orig.ident == 'iv_61','iv','ros'))


Hm_list <- SplitObject(Hm_merged, split.by = "race")
#######Hm_list <- lapply(X = Hm_list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = Hm_list, nfeatures = 6000)
Hm_list <- PrepSCTIntegration(object.list = Hm_list, anchor.features = features)

Hm.anchors <- FindIntegrationAnchors(object.list = Hm_list, normalization.method = "SCT",   anchor.features = features)
Hm.combined.sct <- IntegrateData(anchorset = Hm.anchors, normalization.method = "SCT")

Hm.combined.sct <- RunPCA(Hm.combined.sct,  npcs = 50, verbose = T)
Hm.combined.sct <- RunUMAP(Hm.combined.sct,  dims = 1:4)
Hm.combined.sct <- FindNeighbors(Hm.combined.sct, dims = 1:4)
Hm.combined.sct <- FindClusters(Hm.combined.sct, resolution = 1)




## plot pc1 vs time
#plt_dat <- FetchData(Hm.combined.sct, c("PC_1","PC_2", "time","orig.ident","seurat_clusters"))
#png(file='pc1_vs_time.png')
#ggplot(plt_dat, aes(time, PC_1)) +
#  geom_jitter(aes(color = orig.ident)) +
#  labs(y = "PC1")
#dev.off()

Hm.combined.sct <- PrepSCTFindMarkers(Hm.combined.sct)

SaveH5Seurat(Hm.combined.sct, 'Hm_combined.sct', overwrite = TRUE, verbose = TRUE)
