## load libraries
library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(sctransform)
library(SeuratObject)

#import ivory + melp ts
Hm_merged <- LoadH5Seurat('Hm_combined.sct.h5seurat')

DefaultAssay(Hm_merged)  <- "RNA"
Hm_merged <- NormalizeData(Hm_merged, normalization.method = "LogNormalize", assay = 'RNA')

## look at some plots
 a <- DimPlot(Hm_merged, label=T,repel=T, shuffle=T)
 b <- DimPlot(Hm_merged, label=T,repel=T, group.by='orig.ident', shuffle=T)
 c <- DimPlot(Hm_merged, label=T,repel=T, group.by='time', shuffle=T)
 d <- FeaturePlot(Hm_merged, 'Ivory',order=T)
 e <- FeaturePlot(Hm_merged, 'Hmel.Hmel218003oG779',order=T) #sv
 f <- FeaturePlot(Hm_merged, 'Hmel.Hmel212001oG165',order=T) #sox15
 a+b+c+d+e+f

#making dotplot to select sop lineage
cd_genes <- c("Ivory", 
              "Hmel.Hmel211001oG677", # notch
			  "Hmel.Hmel209001oG161", # delta
              "Hmel.Hmel201001oG607", #Hr38
              "Hmel.Hmel211001oG691", #psGEF
              "Hmel.Hmel201001oG635", #ncad
              "Hmel.Hmel212001oG609", # pdm3
              "Hmel.Hmel207001oG974", #fz2
              "Hmel.Hmel201001oG32",# sns
              "Hmel.Hmel201001oG577",#sens
              "Hmel.Hmel218003oG779", # sv
              "Hmel.Hmel209001oG81", #ss
              "Hmel.Hmel212001oG165", #Sox15
              "Hmel.Hmel201001oG21", #nrm
              "Hmel.Hmel201001oG412") #Su(H)
               

labels <- c('ivory','notch','delta','hr38','psGEF', 'ncad', 'pdm3', 'fz2', 'sns','sens', 'sv', 'ss', 'sox15','nrm','suh')
DotPlot(object = subset(Hm_merged, time==10), features = cd_genes,dot.min=0.1,col.min=0,group.by='sop')+
  theme(axis.text.x = element_text(angle = 30))+
  scale_color_batlow(reverse=T)+
  scale_x_discrete(labels=labels)

sop - 13, some of 12,3,0,17
scale - 1,2,4,6,7,11,19,22,27 (mabye 24, 25)
socket - 5,20,23

## based on gene expression, splitting like this
#scales = c(19,11,1,7,6,22,2,27,4,25,24)
#sockets = c(5,10,18,20,23,14)
##early = c(0,14,8,10,4,15,31,24)
##epi = c(1,9,22,13,17,28,23,21,29)
#
#Hmel.Hmel201001oG635 = ncad
#Hmel.Hmel213001oG861 = sox102f
#Hmel.Hmel201001oG577 = sens
## label the early set into sop and non-sop.

### labelling code
#barcodes <- subset(x = Hm_merged, Hmel.Hmel201001oG635 > 0.5 & Hmel.Hmel213001oG861 < 2 &  time < 25)
#Hm_merged$sop <- ""
#Hm_merged$sop[WhichCells(Hm_merged) %in%  WhichCells(barcodes)] <- "sop"
#Hm_merged$sop[!(WhichCells(Hm_merged) %in%  WhichCells(barcodes)) & Hm_merged$time < 21] <- "epi"
#barcodes <- subset(Hm_merged, seurat_clusters %in% sockets)
#Hm_merged$sop[WhichCells(Hm_merged) %in%  WhichCells(barcodes)] <- "sockets"
#barcodes <- subset(Hm_merged, seurat_clusters %in% scales)
#Hm_merged$sop[WhichCells(Hm_merged) %in%  WhichCells(barcodes)] <- "scales"
#barcodes <- subset(Hm_merged, seurat_clusters %in% c(1,9,22,13,17,28,23,21,29))
#Hm_merged$sop[WhichCells(Hm_merged) %in%  WhichCells(barcodes)] <- "epi"
#Hm_merged$sop[Hm_merged$sop == ""] <- 'unknown'


## trying a new approach to subsetting, becuase i think this previous iteration was bad. 
DefaultAssay(Hm_merged) <- 'RNA'
barcodes <- subset(Hm_merged, seurat_clusters %in% c(13,12,3,0,17,14) & time < 25 & Hmel.Hmel209001oG161 > 0.5)
Hm_merged$sop <- ""
Hm_merged$sop[WhichCells(Hm_merged) %in%  WhichCells(barcodes)] <- "sop"
socket <- c(5,20,23)
scale <- c(1,2,4,6,7,11,19,22,27)
Hm_merged$sop <- ifelse(Hm_merged$seurat_clusters %in% socket, "socket", ifelse(
                         Hm_merged$seurat_clusters %in% scale, "scale", Hm_merged$sop))
Hm_merged$sop[Hm_merged$sop == ""] <- 'unknown'

Idents(Hm_merged) <- 'orig.ident'
iv61_vs_Hm30 <- FindMarkers(Hm_merged, ident.1="iv_61", ident.2=c('Hmr_30'),recorrect_umi=FALSE)



#############################
##### realnalysing sops #####
#############################
Hm_merged$celltypes <- Hm_merged$sop
Idents(Hm_merged) <- Hm_merged$celltypes

Hm_sops <- subset(Hm_merged, celltypes %in% c('sop','socket','scale'))
DimPlot(Hm_sops)

#clusters <- c(0,7,11,18,27,29,2,6,13,24,14,23,25,20,8,28,7,21)
#a <- subset(Hm_merged, idents = c(0,7,11,18,27,29,2,6,13,24,14,23,25,20,8,28,7,21) )
#b <- subset(Hm_merged, sop=='sop' & time < 25)
#Hm_sops <- merge(a,b,merge.data=T, merge.dr = T)
#rm(a,b)
#DimPlot(Hm_sops)

DefaultAssay(Hm_sops) <- 'integrated'

Hm_sops <- FindVariableFeatures(Hm_sops)
#features <- FindVariableFeatures(Hm_sops)
Hm_sops <- ScaleData(Hm_sops, verbose = FALSE)

Hm_sops <- RunPCA(Hm_sops, npcs=50, verbose = FALSE)

##based on elbow plot i'm goign to go for 20 dims
DefaultAssay(Hm_sops) <- 'integrated'
Hm_sops <- RunUMAP(Hm_sops, reduction = "pca", dims = 1:8)
Hm_sops <- FindNeighbors(Hm_sops, reduction = "pca", dims = 1:8)
Hm_sops <- FindClusters(Hm_sops, resolution = 2)

DefaultAssay(Hm_sops) <- 'SCT'
a <- DimPlot(Hm_sops,group.by='time')
b <- FeaturePlot(Hm_sops,'Ivory')
c <- FeaturePlot(Hm_sops,'Hmel.Hmel218003oG779')
d <- FeaturePlot(Hm_sops,'Hmel.Hmel207001oG974')
a+b+c+d

### need to reupload to the cluster to run this bit
SeuratDisk::SaveH5Seurat(Hm_sops, 'Hm_sops', overwrite = TRUE, verbose = TRUE)
Hm_sops <- LoadH5Seurat('Hm_sops*')


# it worked really well with 5 dims and then running 1:3 locally, but trying again with 1:4 here
DefaultAssay(Hm_sops) <- 'RNA'
Hm_sops <- SCTransform(Hm_sops, verbose = T, variable.features.n = 5000)
Hm_sops <- RunPCA(Hm_sops,  npcs = 50, verbose = T)
Hm_sops <- RunUMAP(Hm_sops,  dims = 1:5)
Hm_sops <- FindNeighbors(Hm_sops, dims = 1:5)
Hm_sops <- FindClusters(Hm_sops, resolution = 2)

#Hm_sops$race <- ifelse(Hm_sops$orig.ident == "Hmm_10" | Hm_merged$orig.ident == "Hmm_25" | Hm_merged$orig.ident == 'Hmm_40', "mal", "ros")


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

DefaultAssay(Hm.sops.sct) <- 'RNA'
a <- DimPlot(subset(Hm.sops.sct, time != 56),group.by='celltypes',label=T,repel=T,shuffle=T)+theme_void()+theme_void()+scale_colour_manual(values=okabe5[1:3])
b <- DimPlot(subset(Hm.sops.sct, time != 56),group.by='race',shuffle=T)+theme_void()+scale_colour_manual(values=okabe5[1:3])
c <- FeaturePlot(subset(Hm.sops.sct, time != 56),'Hmel.Hmel218003oG779',min.cutoff=0,order=T) # sv
d <- DimPlot(subset(Hm.sops.sct, time != 56), group.by='Phase',shuffle=T) 
e <- DimPlot(subset(Hm.sops.sct, time != 56),pt.size=1.25,alpha=0.8,group.by='time',shuffle=T,label=T,repel=T)+theme_void()+scale_colour_discreterainbow(reverse=T)
f <- DimPlot(subset(Hm.sops.sct, time != 56), group.by='seurat_clusters',shuffle=T,label=T,repel=T)
a+b+c+d+e+f

# plot per race 
RacePlotCluster <- table(Hm.sops.sct@meta.data$celltypes, Hm.sops.sct@meta.data$race)
hsize <- 1
okabe5 <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
RaceClusters <- ggplot(as.data.frame(RacePlotCluster), aes(fill=Var2, x=hsize, y=Freq))+
           geom_bar(position='fill',stat='identity',)+
           coord_polar("y")+
           xlim(c(0.2,hsize+.5))+
           scale_fill_manual(values=okabe5[1:3])+
           theme_void()+
           theme(axis.text=element_text(colour='black'))+
           facet_grid(~Var1)
## saving plots 
a
dev.copy(pdf, 'Hm.sops.sct.celltypes.pdf',width=5,height=5)
dev.off()
b
dev.copy(pdf, 'Hm.sops.sct.race.pdf',width=5,height=5)
dev.off()
e
dev.copy(pdf, 'Hm.sops.sct.time.pdf',width=5,height=5)
dev.off()

RaceClusters
dev.copy(pdf, 'RaceClusters.pdf',width=7,height=3)
dev.off()

### cell cycle scoring
### import list of canonical cell cycle markers
ccgenes <- read.table('CellCycleGenes.txt',h=T)
### make objects for each class
s.genes <- ccgenes$HmID[ccgenes$stage == 'S']
g2m.genes <- ccgenes$HmID[ccgenes$stage == 'G2/M']
### assign each cell a cell cycle score
Hm.sops.sct <- CellCycleScoring(Hm.sops.sct, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

## find markers
Idents(Hm.sops.sct) <- 'seurat_clusters'
markers <- FindAllMarkers(Hm.sops.sct,recorrect_umi=FALSE, only.pos=TRUE)

#labelling wt vs mut
barcodes <- subset(x = Hm.sops.sct, orig.ident %in% c('iv_61', 'iv_38','iv_55'))
Hm.sops.sct$mut <- ""
Hm.sops.sct$mut[WhichCells(Hm.sops.sct) %in%  WhichCells(barcodes)] <- "mut"
Hm.sops.sct$mut[!(WhichCells(Hm.sops.sct) %in%  WhichCells(barcodes))] <- "wt"


Idents(Hm.sops.sct) <- 'mut'
markersMut <- FindMarkers(Hm.sops.sct, ident.1='wt', ident.2='mut',recorrect_umi=FALSE, pseudocount.use=0,test.use='roc')

markersMutTimedAUC <- FindMarkers(subset(Hm.sops.sct, orig.ident %in% c('Hmr_25','Hmm_25', 'Hmr_30', 'Hmr_35','Hmr_40', 'iv_20', 'iv_61','iv_55') & celltypes == 'scale'), ident.1='wt', ident.2='mut',recorrect_umi=FALSE, pseudocount.use=0,test.use='roc')
markersMutTimed <- FindMarkers(subset(Hm.sops.sct, orig.ident %in% c('Hmr_25','Hmm_25', 'Hmr_30', 'Hmr_35','Hmr_40', 'iv_20', 'iv_61','iv_55') & celltypes == 'scale'), ident.1='wt', ident.2='mut',recorrect_umi=FALSE, pseudocount.use=0)

markersMutTimedAUC$gene <- rownames(markersMutTimedAUC)
markersMutTimed$gene    <- rownames(markersMutTimed)

test <- merge(markersMutTimedAUC, markersMutTimed, by='gene')
test$powerSig <- ifelse(test$power > 0.5, 'yes','no')

> markersMutTimedAUC[markersMutTimedAUC$power > 0.6,] 
                     myAUC   avg_diff power avg_log2FC pct.1 pct.2	#
Hmel.Hmel215003oG93  0.022 -2.4031791 0.956 -3.9918185 0.861 0.997	#
Hmel.Hmel221001oG263 0.069 -1.4656976 0.862 -3.1963669 0.406 0.951	# yin
Hmel.Hmel213001oG244 0.104 -1.0866813 0.792 -1.7145142 0.717 1.000	# yellow-c
Hmel.Hmel221001oG563 0.888  1.8267287 0.776  7.7374264 0.783 0.024	#
Hmel.Hmel200219oG4   0.888  1.0714304 0.776  2.9474047 0.893 0.290	# ??
Hmel.Hmel210001oG623 0.875  2.5700351 0.750  4.3658291 0.874 0.603	# CG1273
Hmel.Hmel215003oG103 0.128 -1.2198908 0.744 -3.2127091 0.357 0.860	# RNaseZ
Hmel.Hmel215003oG101 0.130 -1.1555107 0.740 -4.1045382 0.178 0.800	# DNAAF
Hmel.Hmel215003oG92  0.142 -1.3586467 0.716 -4.0933904 0.233 0.794	# Rpn10
Hmel.Hmel219001oG17  0.149 -0.9066942 0.702 -2.7754634 0.260 0.842	# CG11498
Hmel.Hmel210001oG572 0.850  2.0595412 0.700  3.7784237 0.850 0.608	# nw
Hmel.Hmel215003oG104 0.150 -1.1490588 0.700 -4.1496463 0.170 0.762	#
Hmel.Hmel215003oG90  0.151 -1.0902270 0.698 -2.1159694 0.701 0.930	# 
Hmel.Hmel212001oG397 0.842  0.7073431 0.684  1.1585495 0.998 0.939	# Rhp
Hmel.Hmel200218oG6   0.158 -0.8843517 0.684 -1.4243181 0.960 0.995	#
Hmel.Hmel217001oG182 0.159 -0.8593322 0.682 -1.2980660 0.689 0.999	# CG16786
Hmel.Hmel202001oG522 0.163 -0.8590800 0.674 -1.4572563 0.688 0.997	# CG2663
Hmel.Hmel210001oG474 0.165 -1.4551334 0.670 -2.2302235 0.510 0.971	# exp
Hmel.Hmel219001oG138 0.169 -1.1129557 0.662 -2.4480836 0.403 0.877	# - this is an osiris gene need to figure out which 
Hmel.Hmel220003oG342 0.828  1.9927707 0.656  4.5260188 0.778 0.287	# tnc
Hmel.Hmel217001oG882 0.827  1.2861275 0.654  3.1437461 0.825 0.410	# Dhx15
HMEL000025           0.824  1.5721282 0.648  2.4327653 0.955 0.948	# cortex
Hmel.Hmel219001oG15  0.178 -0.8362310 0.644 -2.4230884 0.371 0.830	# Rsod
Hmel.Hmel200185oG1   0.182 -0.5021995 0.636 -0.7304415 1.000 1.000	#
Hmel.Hmel217001oG932 0.183 -0.7406052 0.634 -1.8536301 0.523 0.863	# CG5890
Hmel.Hmel221001oG385 0.185 -0.6817726 0.630 -1.0405486 0.953 0.999	# trio
Hmel.Hmel215003oG438 0.815  1.3722677 0.630  5.8000747 0.656 0.072	# Ca-alpha1T
Hmel.Hmel208001oG189 0.811  0.8201020 0.622  1.3380722 0.992 0.956	# MP1 - melanization proteinase 
Hmel.Hmel201001oG993 0.190 -0.7993129 0.620 -1.4349056 0.520 0.973	# mfas
Hmel.Hmel213001oG820 0.808  0.7703909 0.616  1.9096350 0.881 0.510	#
Hmel.Hmel216002oG22  0.193 -0.9145537 0.614 -1.8165615 0.492 0.910	# CG14636
Hmel.Hmel207001oG518 0.194 -0.8094457 0.612 -1.2067550 0.983 1.000	# Svil
Hmel.Hmel204001oG224 0.194 -0.6087824 0.612 -0.9872698 0.949 0.998	# Fur1
Hmel.Hmel200135oG1   0.194 -0.5002566 0.612 -0.7296494 1.000 1.000	#
Ivory                0.804  2.6795476 0.608  4.1323815 0.940 0.903	# ivory 
Hmel.Hmel210001oG711 0.804  1.2123372 0.608  2.4322759 0.823 0.609	#
Hmel.Hmel212001oG819 0.197 -0.7095305 0.606 -1.7471856 0.458 0.865	#
Hmel.Hmel216002oG153 0.802  1.1671961 0.604  3.0151220 0.743 0.236	#
> markersMutTimedAUC[markersMutTimedAUC$power > 0.5 & markersMutTimedAUC$power < 0.6000001,] 
                      myAUC   avg_diff power avg_log2FC pct.1 pct.2
Hmel.Hmel210001oG964  0.798  3.3784367 0.596  7.3593457 0.687 0.186	# NLaz
Hmel.Hmel205001oG709  0.795  0.7909304 0.590  1.5308420 0.935 0.794	# Alp2
Hmel.Hmel219001oG134  0.206 -1.0335202 0.588 -1.9560295 0.620 0.919	# Osi10a
Hmel.Hmel201001oG933  0.206 -0.8267175 0.588 -1.4764526 0.842 0.960	# trbl
Hmel.Hmel217001oG710  0.793  0.9015591 0.586  1.8972598 0.917 0.600	# hts
Hmel.Hmel219001oG230  0.212 -0.8669181 0.576 -2.2411097 0.349 0.799	# St1
Hmel.Hmel202001oG423  0.212 -1.1906107 0.576 -6.0353861 0.050 0.601	# CG18586
Hmel.Hmel211001oG115  0.787  0.8712442 0.574  2.2570944 0.787 0.364	# ?? (next to hb)
Hmel.Hmel220003oG695  0.214 -0.7722487 0.572 -1.1761053 0.946 1.000	# sprt
Hmel.Hmel218003oG779  0.782  0.6227650 0.564  0.9273779 1.000 0.999	# sv
Hmel.Hmel209001oG223  0.782  1.1672981 0.564  2.5602107 0.814 0.436	# needs to be collapsed into jupiter. 
Hmel.Hmel210001oG965  0.219 -0.5127390 0.562 -0.8033618 0.971 0.999	# bun
Hmel.Hmel211001oG624  0.221 -0.5282520 0.558 -0.8299366 0.637 0.983	# Gasp
Hmel.Hmel203003oG391  0.221 -0.7764261 0.558 -2.1615725 0.398 0.785	# Xrp1 (poor hit)
Hmel.Hmel219001oG127  0.222 -0.6775205 0.556 -2.5395519 0.225 0.710	# Osi7
Hmel.Hmel200217oG5    0.222 -0.6031747 0.556 -1.3998431 0.661 0.873	# 
Hmel.Hmel207001oG199  0.777  1.4676318 0.554  3.0355238 0.811 0.465	# Np
Hmel.Hmel207001oG195  0.777  2.2476750 0.554  3.8112444 0.831 0.587	# Hr3
Hmel.Hmel201001oG1041 0.776  1.1345761 0.552  4.6716156 0.611 0.091	# loaf
Hmel.Hmel213001oG112  0.225 -1.1966952 0.550 -2.6469345 0.555 0.817	# Ctns
Hmel.Hmel218003oG36   0.774  1.0671401 0.548  2.0108193 0.938 0.707	# CG42788
Hmel.Hmel213001oG830  0.773  0.4915099 0.546  0.8275358 0.997 0.979	# CG7974
Hmel.Hmel220003oG414  0.228 -0.4894017 0.544 -0.8189051 0.933 0.996	# Rbfox1
Hmel.Hmel219001oG547  0.231 -0.6002848 0.538 -1.3850249 0.480 0.877	# CG17646
Hmel.Hmel214004oG555  0.231 -0.7474543 0.538 -1.6461153 0.592 0.875	# dpr21
Hmel.Hmel206001oG762  0.232 -0.5561910 0.536 -0.9332123 0.955 0.995	# CycG
Hmel.Hmel201001oG627  0.767  0.8714801 0.534  3.3461130 0.630 0.135	# Phlpp
Hmel.Hmel212001oG886  0.765  1.1702700 0.530  2.1526325 0.897 0.736	# mth
Hmel.Hmel211001oG698  0.235 -0.5715971 0.530 -1.0782770 0.808 0.963	# CG34049_CG42764
Hmel.Hmel214004oG525  0.763  0.9990187 0.526  2.3175891 0.764 0.389	# wun
Hmel.Hmel213001oG298  0.763  0.7544347 0.526  1.5695131 0.881 0.568	# Timp
Hmel.Hmel210001oG713  0.237 -0.5583926 0.526 -1.2393120 0.678 0.909	# CG34347
Hmel.Hmel216002oG565  0.240 -0.5923899 0.520 -1.0027230 0.825 0.987	# C901
Hmel.Hmel221001oG662  0.758  0.5263269 0.516  1.2392825 0.887 0.569	# CG1764_CG6770
Hmel.Hmel221001oG222  0.243 -0.7129067 0.514 -1.4976123 0.662 0.873	# path
Hmel.Hmel210001oG148  0.755  2.0065863 0.510  5.5054412 0.593 0.126	# Mco1
Hmel.Hmel215003oG50   0.245 -0.8337361 0.510 -2.0115011 0.355 0.762	# CG13806
Hmel.Hmel212001oG1007 0.245 -0.5489115 0.510 -0.9263537 0.815 0.996	# Cht6
Hmel.Hmel209001oG229  0.754  1.2435019 0.508  1.8554842 0.992 0.989	# Jupiter
Hmel.Hmel212001oG682  0.752  0.7477344 0.504  7.3274067 0.509 0.006	# CG6204
Hmel.Hmel212001oG74   0.752  0.9971556 0.504  1.9052470 0.853 0.534	# CG32354
Hmel.Hmel213001oG461  0.248 -0.6303377 0.504 -1.2570742 0.674 0.911	# CG13562_CG6847
Hmel.Hmel200218oG19   0.248 -0.6249415 0.504 -1.2601195 0.770 0.890	# 
Hmel.Hmel212001oG1047 0.751  0.5234972 0.502  0.9132516 0.993 0.952	# DIP1_blanks

                     	#
Hmel.Hmel215003oG93  	# LMTK1
Hmel.Hmel221001oG263 	# CG2930_yin
Hmel.Hmel213001oG244 	# yellow-c
Hmel.Hmel221001oG563 	# Glug_CG31100
Hmel.Hmel200219oG4   	# ?? -Tnp_P_element super family
Hmel.Hmel210001oG623 	# CG1273
Hmel.Hmel215003oG103 	# RNaseZ
Hmel.Hmel215003oG101 	# Dnaaf5
Hmel.Hmel215003oG92  	# Rpn10
Hmel.Hmel219001oG17  	# CG11498
Hmel.Hmel210001oG572 	# nw
Hmel.Hmel215003oG104 	# (nested in LMTK, antiparallel. )
Hmel.Hmel215003oG90  	# parn
Hmel.Hmel212001oG397 	# Rhp
Hmel.Hmel200218oG6   	# iris? (poor hit)
Hmel.Hmel217001oG182 	# CG16786
Hmel.Hmel202001oG522 	# CG2663
Hmel.Hmel210001oG474 	# exp
Hmel.Hmel219001oG138 	# Osi12
Hmel.Hmel220003oG342 	# tnc
Hmel.Hmel217001oG882 	# Dhx15
HMEL000025           	# cortex
Hmel.Hmel219001oG15  	# Rsod
Hmel.Hmel200185oG1   	# ??
Hmel.Hmel217001oG932 	# CG5890
Hmel.Hmel221001oG385 	# trio
Hmel.Hmel215003oG438 	# Ca-alpha1T
Hmel.Hmel208001oG189 	# MP1 - melanisation protein 1
Hmel.Hmel201001oG993 	# mfas
Hmel.Hmel213001oG820 	# Pif1
Hmel.Hmel216002oG22  	# CG14636
Hmel.Hmel207001oG518 	# Svil
Hmel.Hmel204001oG224 	# Fur1
Hmel.Hmel200135oG1   	# (poor hit to polybromo)
Ivory                	# Ivory 
Hmel.Hmel210001oG711 	# hui_mtm (poor hit to both)
Hmel.Hmel212001oG819 	# Duox (aka curly)
Hmel.Hmel216002oG153 	# cyr
Hmel.Hmel210001oG964  	# NLaz
Hmel.Hmel205001oG709  	# Alp2
Hmel.Hmel219001oG134  	# Osi10a
Hmel.Hmel201001oG933  	# trbl
Hmel.Hmel217001oG710  	# hts
Hmel.Hmel219001oG230  	# St1
Hmel.Hmel202001oG423  	# CG18586
Hmel.Hmel211001oG115  	# ?? (next to hb)
Hmel.Hmel220003oG695  	# sprt
Hmel.Hmel218003oG779  	# sv
Hmel.Hmel209001oG223  	# needs to be collapsed into jupiter. 
Hmel.Hmel210001oG965  	# bun
Hmel.Hmel211001oG624  	# Gasp
Hmel.Hmel203003oG391  	# Xrp1 (poor hit)
Hmel.Hmel219001oG127  	# Osi7
Hmel.Hmel200217oG5    	# 
Hmel.Hmel207001oG199  	# Np
Hmel.Hmel207001oG195  	# Hr3
Hmel.Hmel201001oG1041 	# loaf
Hmel.Hmel213001oG112  	# Ctns
Hmel.Hmel218003oG36   	# CG42788
Hmel.Hmel213001oG830  	# CG7974
Hmel.Hmel220003oG414  	# Rbfox1
Hmel.Hmel219001oG547  	# CG17646
Hmel.Hmel214004oG555  	# dpr21
Hmel.Hmel206001oG762  	# CycG
Hmel.Hmel201001oG627  	# Phlpp
Hmel.Hmel212001oG886  	# mth
Hmel.Hmel211001oG698  	# CG34049_CG42764
Hmel.Hmel214004oG525  	# wun
Hmel.Hmel213001oG298  	# Timp
Hmel.Hmel210001oG713  	# CG34347
Hmel.Hmel216002oG565  	# C901
Hmel.Hmel221001oG662  	# CG1764_CG6770
Hmel.Hmel221001oG222  	# path
Hmel.Hmel210001oG148  	# Mco1
Hmel.Hmel215003oG50   	# CG13806
Hmel.Hmel212001oG1007 	# Cht6
Hmel.Hmel209001oG229  	# Jupiter
Hmel.Hmel212001oG682  	# CG6204
Hmel.Hmel212001oG74   	# CG32354
Hmel.Hmel213001oG461  	# CG13562_CG6847
Hmel.Hmel200218oG19   	# ??
Hmel.Hmel212001oG1047 	# DIP1_blanks

markersMutTimedEarlyAUC <- FindMarkers(subset(Hm.sops.sct, orig.ident %in% c('Hmr_25','Hmm_25', 'Hmr_30', 'iv_20', 'iv_61')), ident.1='wt', ident.2='mut',recorrect_umi=FALSE, pseudocount.use=0,test.use='roc')
markersMutTimedEarly <- FindMarkers(subset(Hm.sops.sct, orig.ident %in% c('Hmr_25','Hmm_25', 'Hmr_30', 'iv_20', 'iv_61')), ident.1='wt', ident.2='mut',recorrect_umi=FALSE, pseudocount.use=0)

markersMutTimedEarlyAUC$gene <- rownames(markersMutTimedEarlyAUC)
markersMutTimedEarly$gene    <- rownames(markersMutTimedEarly)

markersEarlyMerged <- merge(markersMutTimedEarlyAUC, markersMutTimedEarly, by='gene')
markersEarlyMerged$powerSig <- ifelse(markersEarlyMerged$power > 0.6, 'yes','no')
                gene 		myAUC   avg_diff power avg_log2FC.x pct.1.x pct.2.x          p_val avg_log2FC.y pct.1.y pct.2.y     p_val_adj powerSig
  Hmel.Hmel200185oG1 		0.054 -0.7354552 0.892   -1.0687092   1.000   1.000   0.000000e+00   -1.0687092   1.000   1.000  0.000000e+00      yes
  Hmel.Hmel200219oG4 		0.813  0.8158551 0.626    2.1840945   0.826   0.369  3.430578e-202    2.1840945   0.826   0.369 5.330775e-198      yes
Hmel.Hmel201001oG657 	Hmel.Hmel201001oG657	Osi21	0.127 -1.0883290 0.746   -1.6823174   0.567   1.000  2.015260e-279   -1.6823174   0.567   1.000 3.131512e-275      yes
Hmel.Hmel201001oG790 	Hmel.Hmel201001oG790	CG6055	0.102 -1.1661184 0.796   -1.8291310   0.513   1.000   0.000000e+00   -1.8291310   0.513   1.000  0.000000e+00      yes
Hmel.Hmel201001oG802 	Hmel.Hmel201001oG802	Pka-C3	0.145 -1.0333374 0.710   -2.2933182   0.360   0.889  2.795065e-271   -2.2933182   0.360   0.889 4.343251e-267      yes
Hmel.Hmel201001oG993 	Hmel.Hmel201001oG993	mfas	0.140 -1.0829469 0.720   -2.0537005   0.356   0.940  9.006501e-277   -2.0537005   0.356   0.940 1.399520e-272      yes
Hmel.Hmel202001oG423 	Hmel.Hmel202001oG423	CG18586	0.159 -1.2817286 0.682   -4.3202691   0.094   0.731   0.000000e+00   -4.3202691   0.094   0.731  0.000000e+00      yes
Hmel.Hmel202001oG522 	Hmel.Hmel202001oG522	CG2663	0.177 -0.8906459 0.646   -1.5995358   0.518   0.955  3.336490e-213   -1.5995358   0.518   0.955 5.184571e-209      yes
Hmel.Hmel203003oG148 	Hmel.Hmel203003oG148	dy	0.149 -1.0631101 0.702   -1.7704998   0.438   0.986  5.430269e-255   -1.7704998   0.438   0.986 8.438095e-251      yes
Hmel.Hmel203003oG153 	Hmel.Hmel203003oG153	m	0.157 -0.7870951 0.686   -1.1758853   0.930   1.000  8.837886e-233   -1.1758853   0.930   1.000 1.373319e-228      yes
Hmel.Hmel203003oG386 	Hmel.Hmel203003oG386	CG7800_CG7896	0.145 -0.8854873 0.710   -1.3409073   0.991   0.992  2.968305e-248   -1.3409073   0.991   0.992 4.612449e-244      yes
Hmel.Hmel204001oG224 	Hmel.Hmel204001oG224	Fur1	0.166 -0.6813807 0.668   -1.1000958   0.938   0.993  7.798115e-220   -1.1000958   0.938   0.993 1.211749e-215      yes
Hmel.Hmel204001oG515 	Hmel.Hmel204001oG515	futsch	0.199 -0.8887387 0.602   -1.5200584   0.653   0.941  2.552109e-181   -1.5200584   0.653   0.941 3.965722e-177      yes
Hmel.Hmel205001oG273 	Hmel.Hmel205001oG273	CG3961	0.125 -0.9965298 0.750   -2.3723500   0.353   0.912  1.562460e-301   -2.3723500   0.353   0.912 2.427906e-297      yes
Hmel.Hmel206001oG248 	Hmel.Hmel206001oG248	gukh	0.170 -0.8421710 0.660   -1.9862560   0.394   0.896  1.366557e-231   -1.9862560   0.394   0.896 2.123494e-227      yes
 Hmel.Hmel206001oG48 	Hmel.Hmel206001oG48     Amph	0.181 -0.7663976 0.638   -1.2849203   0.760   0.993  6.186114e-202   -1.2849203   0.760   0.993 9.612602e-198      yes
Hmel.Hmel206001oG635 	Hmel.Hmel206001oG635	LManI	0.115 -1.1145714 0.770   -2.5305330   0.309   0.932   0.000000e+00   -2.5305330   0.309   0.932  0.000000e+00      yes
Hmel.Hmel206001oG777 	Hmel.Hmel206001oG777	Invadolysin	0.153 -1.1720152 0.694   -1.8947930   0.509   0.918  2.072123e-247   -1.8947930   0.509   0.918 3.219872e-243      yes
Hmel.Hmel206001oG841 	Hmel.Hmel206001oG841 	CG5003	0.103 -1.3575504 0.794   -2.0614360   0.544   0.985   0.000000e+00   -2.0614360   0.544   0.985  0.000000e+00      yes
Hmel.Hmel207001oG195 	Hmel.Hmel207001oG195	Hr3	0.806  1.6110267 0.612    2.4888863   0.964   0.756  2.214630e-185    2.4888863   0.964   0.756 3.441313e-181      yes
Hmel.Hmel207001oG199 	Hmel.Hmel207001oG199	Np	0.825  1.3223684 0.650    2.4789326   0.908   0.625  3.744604e-210    2.4789326   0.908   0.625 5.818740e-206      yes
Hmel.Hmel207001oG518 	Hmel.Hmel207001oG518	Svil	0.074 -1.2057516 0.852   -1.7977570   0.963   1.000   0.000000e+00   -1.7977570   0.963   1.000  0.000000e+00      yes
Hmel.Hmel207001oG519 	Hmel.Hmel207001oG519	Svil1	0.179 -0.7894799 0.642   -1.7716205   0.492   0.899  9.712977e-213   -1.7716205   0.492   0.899 1.509300e-208      yes
Hmel.Hmel208001oG189 	Hmel.Hmel208001oG189	MP1 - melanisation protein 1	0.881  1.1455341 0.762    1.9276160   0.991   0.933  1.865145e-285    1.9276160   0.991   0.933 2.898249e-281      yes
Hmel.Hmel209001oG433 							???	0.166 -0.8272209 0.668   -1.3325916   0.726   0.993  1.540525e-221   -1.3325916   0.726   0.993 2.393822e-217      yes
 Hmel.Hmel210001oG19 	Hmel.Hmel210001oG19	CG15080	0.183 -0.8770124 0.634   -1.3926884   0.616   0.983  3.731009e-202   -1.3926884   0.616   0.983 5.797614e-198      yes
Hmel.Hmel210001oG474 	Hmel.Hmel210001oG474	exp	0.090 -1.6564869 0.820   -3.0377067   0.382   0.952   0.000000e+00   -3.0377067   0.382   0.952  0.000000e+00      yes
Hmel.Hmel210001oG572 	Hmel.Hmel210001oG572	nw	0.979  2.5606111 0.958    4.4429830   1.000   0.627   0.000000e+00    4.4429830   1.000   0.627  0.000000e+00      yes
Hmel.Hmel210001oG623 	Hmel.Hmel210001oG623	CG1273	0.960  2.3552981 0.920    3.8885316   0.990   0.749   0.000000e+00    3.8885316   0.990   0.749  0.000000e+00      yes
Hmel.Hmel210001oG642 	Hmel.Hmel210001oG642	dyl	0.166 -0.9235096 0.668   -1.4255329   0.531   1.000  1.638321e-225   -1.4255329   0.531   1.000 2.545786e-221      yes
Hmel.Hmel210001oG711 	Hmel.Hmel210001oG878	CG15239	0.915  1.2595376 0.830    2.3495429   0.989   0.712   0.000000e+00    2.3495429   0.989   0.712  0.000000e+00      yes
Hmel.Hmel210001oG878 	Hmel.Hmel210001oG964	CG31659	0.165 -0.7432737 0.670   -1.4160684   0.679   0.981  4.412138e-223   -1.4160684   0.679   0.981 6.856021e-219      yes
Hmel.Hmel210001oG964 	Hmel.Hmel210001oG964	NLaz	0.889  3.7419827 0.778    7.6039429   0.847   0.196   0.000000e+00    7.6039429   0.847   0.196  0.000000e+00      yes
Hmel.Hmel210001oG965 	Hmel.Hmel210001oG965	bun	0.161 -0.6356279 0.678   -1.0085305   0.956   0.999  2.118523e-226   -1.0085305   0.956   0.999 3.291973e-222      yes
Hmel.Hmel211001oG115 	Hmel.Hmel211001oG115  	# ?? (next to hb)	0.853  1.1148207 0.706    2.3859058   0.933   0.522  7.292576e-249    2.3859058   0.933   0.522 1.133193e-244      yes
Hmel.Hmel211001oG624 	Hmel.Hmel211001oG624	Gasp    0.150 -1.1113910 0.700   -2.0820539   0.535   0.931  3.233709e-249   -2.0820539   0.535   0.931 5.024861e-245      yes
Hmel.Hmel211001oG705 	Hmel.Hmel211001oG705	pot		0.188 -0.7598877 0.624   -1.2157565   0.841   0.990  2.224977e-192   -1.2157565   0.841   0.990 3.457392e-188      yes
Hmel.Hmel211001oG762 	Hmel.Hmel211001oG762	Eip74EF	0.825  2.4016995 0.650    4.9372830   0.803   0.325  5.349257e-221    4.9372830   0.803   0.325 8.312210e-217      yes
Hmel.Hmel212001oG886 	Hmel.Hmel212001oG886 	mth	0.846  1.2965827 0.692    2.2305783   0.983   0.854  3.107434e-236    2.2305783   0.983   0.854 4.828641e-232      yes
Hmel.Hmel213001oG244 	Hmel.Hmel213001oG244	yellow-c	0.121 -1.1265214 0.758   -1.7905746   0.624   0.992  7.087003e-287   -1.7905746   0.624   0.992 1.101249e-282      yes
Hmel.Hmel213001oG325 	Hmel.Hmel213001oG325	Blimp-1	0.157 -0.7498146 0.686   -1.1881506   0.924   0.999  1.077825e-232   -1.1881506   0.924   0.999 1.674832e-228      yes
Hmel.Hmel213001oG502 	Hmel.Hmel213001oG502	MFS9	0.177 -1.0945338 0.646   -1.9151309   0.713   0.927  5.241097e-208   -1.9151309   0.713   0.927 8.144141e-204      yes
Hmel.Hmel213001oG590 		0.806  1.0200931 0.612    2.9699800   0.757   0.244  1.760828e-201    2.9699800   0.757   0.244 2.736151e-197      yes
Hmel.Hmel213001oG624 	Hmel.Hmel213001oG624	NT1	0.197 -0.7643622 0.606   -1.3160990   0.692   0.948  6.934322e-184   -1.3160990   0.692   0.948 1.077524e-179      yes
 Hmel.Hmel213001oG68 	Hmel.Hmel213001oG68	Lmpt	0.131 -1.1683245 0.738   -2.0692033   0.572   0.945  1.413522e-275   -2.0692033   0.572   0.945 2.196472e-271      yes
Hmel.Hmel214004oG555 	Hmel.Hmel214004oG555	dpr21	0.175 -0.9820963 0.650   -2.2906581   0.474   0.847  8.371121e-221   -2.2906581   0.474   0.847 1.300788e-216      yes
Hmel.Hmel215003oG101 	Hmel.Hmel215003oG101	Dnaaf5	0.183 -1.0115997 0.634   -3.6851646   0.158   0.707  3.996392e-266   -3.6851646   0.158   0.707 6.209994e-262      yes
Hmel.Hmel215003oG103 	Hmel.Hmel215003oG103	RNaseZ	0.196 -0.9914441 0.608   -2.7387921   0.335   0.758  1.842334e-209   -2.7387921   0.335   0.758 2.862802e-205      yes
Hmel.Hmel215003oG362 	Hmel.Hmel215003oG362	ect	0.134 -0.8396799 0.732   -1.2559634   0.828   1.000  1.739878e-264   -1.2559634   0.828   1.000 2.703597e-260      yes
Hmel.Hmel215003oG471 	Hmel.Hmel215003oG471	wake	0.834  0.9024580 0.668    2.0942912   0.903   0.454  3.218579e-224    2.0942912   0.903   0.454 5.001350e-220      yes
 Hmel.Hmel215003oG50    Hmel.Hmel215003oG50 	CG13806	0.134 -1.1198747 0.732   -2.5071978   0.315   0.894  1.365502e-293   -2.5071978   0.315   0.894 2.121854e-289      yes
 Hmel.Hmel215003oG55 	Hmel.Hmel215003oG55	Cht7	0.137 -1.1643289 0.726   -1.8164811   0.469   0.983  2.761616e-270   -1.8164811   0.469   0.983 4.291276e-266      yes
 Hmel.Hmel215003oG93 		0.077 -1.6911822 0.846   -2.9153184   0.832   0.977   0.000000e+00   -2.9153184   0.832   0.977  0.000000e+00      yes
Hmel.Hmel216002oG269 	Hmel.Hmel216002oG269	CG13377	0.144 -1.1760905 0.712   -1.8402023   0.582   0.952  4.262217e-255   -1.8402023   0.582   0.952 6.623060e-251      yes
Hmel.Hmel216002oG536 	?????                         0.145 -0.9269306 0.710   -2.0809439   0.355   0.928  1.615534e-269   -2.0809439   0.355   0.928 2.510379e-265      yes
Hmel.Hmel216002oG565 	Hmel.Hmel216002oG565	C901		0.178 -0.8625102 0.644   -1.4518151   0.751   0.983  1.047471e-205   -1.4518151   0.751   0.983 1.627665e-201      yes
Hmel.Hmel217001oG182 	Hmel.Hmel217001oG182	CG16786	0.101 -1.3775188 0.798   -2.1259610   0.604   0.995   0.000000e+00   -2.1259610   0.604   0.995  0.000000e+00      yes
Hmel.Hmel217001oG262 	Hmel.Hmel217001oG262	ovo	0.840  1.6767723 0.680    3.6185707   0.828   0.367  5.612592e-238    3.6185707   0.828   0.367 8.721407e-234      yes
 Hmel.Hmel218003oG36 	Hmel.Hmel218003oG36	CG42788	0.831  1.1913191 0.662    2.2696796   0.934   0.628  9.659259e-218    2.2696796   0.934   0.628 1.500952e-213      yes
Hmel.Hmel218003oG722 	Hmel.Hmel218003oG722	CG33978	0.056 -1.6124094 0.888   -2.7106333   0.518   0.996   0.000000e+00   -2.7106333   0.518   0.996  0.000000e+00      yes
Hmel.Hmel218003oG735 	Hmel.Hmel218003oG735	vri	0.190 -0.7430572 0.620   -1.3919040   0.525   0.976  8.689808e-196   -1.3919040   0.525   0.976 1.350309e-191      yes
Hmel.Hmel218003oG749 	Hmel.Hmel218003oG749	rn	0.803  0.6550289 0.606    0.9851219   0.999   1.000  6.121003e-182    0.9851219   0.999   1.000 9.511427e-178      yes
Hmel.Hmel219001oG125 	Hmel.Hmel219001oG125	Osi4	0.136 -1.0715077 0.728   -1.5947298   0.780   0.998  2.585905e-262   -1.5947298   0.780   0.998 4.018238e-258      yes
Hmel.Hmel219001oG127 	Hmel.Hmel219001oG127	Osi7	0.117 -1.0689073 0.766   -4.0503322   0.113   0.821   0.000000e+00   -4.0503322   0.113   0.821  0.000000e+00      yes
Hmel.Hmel219001oG144 	Hmel.Hmel219001oG144	Osi17	0.166 -0.8561321 0.668   -1.4879408   0.500   0.993  3.830161e-227   -1.4879408   0.500   0.993 5.951688e-223      yes
 Hmel.Hmel219001oG15 	Hmel.Hmel219001oG15	Rsod	0.139 -0.9980891 0.722   -2.8679615   0.288   0.838  1.888549e-295   -2.8679615   0.288   0.838 2.934616e-291      yes
 Hmel.Hmel219001oG17 	Hmel.Hmel219001oG17	CG11498	0.178 -0.7715433 0.644   -1.8666169   0.266   0.836  1.291123e-237   -1.8666169   0.266   0.836 2.006276e-233      yes
Hmel.Hmel219001oG294 	Hmel.Hmel219001oG294	Eip93F	0.874  0.7970271 0.748    1.1772623   1.000   1.000  2.574798e-276    1.1772623   1.000   1.000 4.000979e-272      yes
Hmel.Hmel219001oG657 	Hmel.Hmel219001oG657    CG6142	0.097 -1.1859191 0.806   -2.3580803   0.374   0.982   0.000000e+00   -2.3580803   0.374   0.982  0.000000e+00      yes
Hmel.Hmel220003oG296 	                        CG9257	0.824  1.3587484 0.648    3.7160294   0.758   0.225  7.176185e-227    3.7160294   0.758   0.225 1.115107e-222      yes
Hmel.Hmel220003oG342 	Hmel.Hmel220003oG342	tnc	0.915  2.6248215 0.830    4.9762101   0.948   0.393   0.000000e+00    4.9762101   0.948   0.393  0.000000e+00      yes
Hmel.Hmel220003oG414 	Hmel.Hmel220003oG414	Rbfox1	0.173 -0.6916363 0.654   -1.1708811   0.870   0.986  8.881684e-212   -1.1708811   0.870   0.986 1.380125e-207      yes
Hmel.Hmel221001oG263 	Hmel.Hmel221001oG263	yin	0.113 -1.3738747 0.774   -3.2754589   0.273   0.8795  0.000000e+00   -3.2754589   0.273   0.879  0.000000e+00      yes
Hmel.Hmel221001oG385 	Hmel.Hmel221001oG385	trio	0.112 -1.4002963 0.776   -2.2217327   0.870   0.9731 3.029151e-297   -2.2217327   0.870   0.973 4.706997e-293      yes
Hmel.Hmel221001oG563 	Hmel.Hmel221001oG563	Glug_CG31100	0.816  1.6954509 0.632    7.1059888   0.646   0.0302 6.930462e-247    7.1059888   0.646   0.030 1.076925e-242      yes
          HMEL000025 		0.816  1.5786346 0.632    2.4593176   0.942   0.9233 5.883465e-197    2.4593176   0.942   0.923 9.142316e-193      yes
               Ivory 		0.822  2.6510070 0.644    4.1976853   0.920   0.7881 7.931174e-206    4.1976853   0.920   0.788 1.232425e-201      yes



markersMutTimedMidAUC <- FindMarkers(subset(Hm.sops.sct, orig.ident %in% c('iv_20','iv_55','Hmr_35','Hmr_40','Hmm_40')), ident.1='wt', ident.2='mut',recorrect_umi=FALSE, pseudocount.use=0,test.use='roc')
markersMutTimedMid <- FindMarkers(subset(Hm.sops.sct, orig.ident %in% c('iv_20','iv_55','Hmr_35','Hmr_40','Hmm_40')), ident.1='wt', ident.2='mut',recorrect_umi=FALSE, pseudocount.use=0)

markersMutTimedMidAUC$gene <- rownames(markersMutTimedMidAUC)
markersMutTimedMid$gene    <- rownames(markersMutTimedMid)

markersMidMerged <- merge(markersMutTimedMidAUC, markersMutTimedMid, by='gene')
markersMidMerged$powerSig <- ifelse(markersMidMerged$power > 0.6, 'yes','no')

                     myAUC  avg_diff power avg_log2FC pct.1 pct.2
Hmel.Hmel200135oG2   0.000 -1.594887 1.000  -2.440909 0.989 1.000
Hmel.Hmel215003oG93  0.004 -2.890817 0.992  -4.638117 0.911 1.000
Hmel.Hmel201001oG657 0.990  2.103713 0.980   3.689927 0.998 0.734
Hmel.Hmel207001oG212 0.011 -3.064369 0.978  -5.816750 0.181 0.994
Hmel.Hmel203003oG365 0.013 -1.558272 0.974  -5.568042 0.079 0.988
Hmel.Hmel207001oG179 0.018 -2.072458 0.964  -3.713930 0.359 0.998
Hmel.Hmel211001oG762 0.022 -1.871981 0.956  -2.724094 0.784 1.000
Hmel.Hmel200219oG4   0.976  1.501455 0.952   5.332693 0.966 0.123
Hmel.Hmel212001oG261 0.024 -2.765159 0.952  -7.257254 0.036 0.959
Hmel.Hmel221001oG563 0.975  1.905648 0.950   9.583168 0.952 0.011
Hmel.Hmel213001oG837 0.027 -2.004843 0.946  -4.569484 0.097 0.985
Hmel.Hmel218003oG32  0.028 -1.693778 0.944  -2.827662 0.701 0.998
Hmel.Hmel211001oG760 0.028 -2.130049 0.944  -3.931680 0.194 0.992
Hmel.Hmel215003oG92  0.029 -1.959294 0.942  -5.371618 0.218 0.962
Hmel.Hmel215003oG438 0.969  1.607628 0.938   5.657393 0.952 0.134
Hmel.Hmel218003oG716 0.031 -1.944430 0.938  -4.216828 0.309 0.992
Hmel.Hmel201001oG241 0.031 -2.736831 0.938  -4.307387 0.887 1.000
Hmel.Hmel215003oG104 0.036 -1.615214 0.928  -5.159224 0.183 0.950
Hmel.Hmel217001oG590 0.037 -2.329176 0.926  -5.674090 0.186 0.950
Hmel.Hmel211001oG590 0.038 -2.187904 0.924  -5.822595 0.141 0.944
Hmel.Hmel208001oG305 0.038 -1.749555 0.924  -2.577097 0.554 1.000
Hmel.Hmel201001oG566 0.038 -2.381852 0.924  -3.750678 0.617 0.997
Hmel.Hmel212001oG588 0.040 -1.905606 0.920  -4.009245 0.300 0.974
Hmel.Hmel210001oG949 0.040 -2.072479 0.920  -3.826589 0.484 0.980
Hmel.Hmel212001oG540 0.956  2.207295 0.912   3.559009 0.994 0.780
Hmel.Hmel213001oG322 0.045 -1.375923 0.910  -2.208767 0.834 1.000
Hmel.Hmel220003oG761 0.047 -2.169449 0.906  -5.053921 0.208 0.944
Hmel.Hmel220003oG328 0.047 -1.698047 0.906  -2.542903 0.737 1.000
Hmel.Hmel212001oG103 0.049 -1.672372 0.902  -2.982295 0.366 0.992
Hmel.Hmel211001oG769 0.049 -1.491034 0.902  -3.088967 0.282 0.985
Hmel.Hmel221001oG788 0.951  1.258532 0.902   2.113471 0.994 0.947
                 gene	
   Hmel.Hmel200135oG2	
  Hmel.Hmel215003oG93	
 Hmel.Hmel201001oG657	Osi21
 Hmel.Hmel207001oG212	CG30091
 Hmel.Hmel203003oG365	CG11836
 Hmel.Hmel207001oG179	CG8170
 Hmel.Hmel211001oG762	Eip74EF
   Hmel.Hmel200219oG4	
 Hmel.Hmel212001oG261	peb
 Hmel.Hmel221001oG563	Glug_CG31100
 Hmel.Hmel213001oG837	CG5731
  Hmel.Hmel218003oG32	
 Hmel.Hmel211001oG760	Raf
  Hmel.Hmel215003oG92	Rpn10
 Hmel.Hmel215003oG438	Ca-alpha1T
 Hmel.Hmel218003oG716	unannotated neighbours pan
 Hmel.Hmel201001oG241	sm
 Hmel.Hmel215003oG104	# (nested in LMTK, antiparallel. )
 Hmel.Hmel217001oG590	jing
 Hmel.Hmel211001oG590	Hr4
 Hmel.Hmel208001oG305	??
 Hmel.Hmel201001oG566	CG2225
 Hmel.Hmel212001oG588	Picot
 Hmel.Hmel210001oG949	Cda5
 Hmel.Hmel212001oG540	yellow-h
 Hmel.Hmel213001oG322	Blimp-1
 Hmel.Hmel220003oG761	next to fj
 Hmel.Hmel220003oG328	Eip75B
 Hmel.Hmel212001oG103	?? next to wls
 Hmel.Hmel211001oG769	Lsm11
 Hmel.Hmel221001oG788	fax


#names <- read.table("gene_names_for_volcplot.txt",h=F)
colnames(names) <- c("gene", "gene_name")
test <- merge(test,names, by='gene',all.x=T)
test$gene_name[test$gene_name==NA] <- ""
a <- ggplot(markersMidMerged, aes(-log10(p_val_adj),avg_log2FC.y, colour=powerSig))+
   geom_point(alpha=0.4)+
   scale_colour_manual(values = c('black','blue'))+
 #  geom_text_repel(data = filter(test, !is.na(gene_name)), aes(label = gene_name),size = 3,max.overlaps = Inf)+
   theme_classic()+ylim(-10,10)

b <- ggplot(markersEarlyMerged, aes(-log10(p_val_adj),avg_log2FC.y, colour=powerSig))+
   geom_point(alpha=0.4)+
   scale_colour_manual(values = c('black','blue'))+
 #  geom_text_repel(data = filter(test, !is.na(gene_name)), aes(label = gene_name),size = 3,max.overlaps = Inf)+
   theme_classic()+ylim(-10,10)


Idents(Hm_merged) <- Hm_merged$orig.ident
markersMutLateAUC <- FindMarkers(subset(Hm_merged), ident.1=c('Hmr_45','Hmr_50'), ident.2='iv_38',recorrect_umi=FALSE, pseudocount.use=0,test.use='roc')
markersMutLate <- FindMarkers(subset(Hm_merged), ident.1=c('Hmr_45','Hmr_50'), ident.2='iv_38',recorrect_umi=FALSE, pseudocount.use=0)

markersMutLateAUC$gene <- rownames(markersMutLateAUC)
markersMutLate$gene    <- rownames(markersMutLate)

markersLateMerged <- merge(markersMutLateAUC, markersMutLate, by='gene')
markersLateMerged$powerSig <- ifelse(markersLateMerged$power > 0.5, 'yes','no')
                      myAUC   avg_diff power avg_log2FC pct.1 pct.2
Hmel.Hmel200012oG2	0.959	1.394721	0.918	2.2725	0.996	0.938		
Hmel.Hmel200012oG3	0.979	1.7049443	0.958	5.070076	0.975	0.149		
Hmel.Hmel200012oG4	0.901	1.2073153	0.802	4.600495	0.84	0.121		
Hmel.Hmel200012oG5	0.943	1.5364439	0.886	3.980367	0.953	0.27		
Hmel.Hmel200075oG1	0.087	-1.1736048	0.826	#NAME?	0	0.827		
Hmel.Hmel200135oG2	0	-1.7486599	1	-2.687472	0.967	1		nxt_to_Dnai2
Hmel.Hmel201001oG1075	0.905	1.1460759	0.81	5.756179	0.825	0.042		
Hmel.Hmel201001oG163	0.905	1.1796339	0.81	3.170659	0.902	0.173	Hmel.Hmel201001oG163	RyR
Hmel.Hmel202001oG239	0.988	1.8300161	0.976	5.783603	0.985	0.08		
Hmel.Hmel202001oG4	0.083	-1.3402339	0.834	-5.015678	0.084	0.862		
Hmel.Hmel202001oG461	0.904	1.0667205	0.808	5.292419	0.833	0.059		
Hmel.Hmel202001oG567	0.927	1.1131478	0.854	2.124149	0.989	0.72	Hmel.Hmel202001oG567	cup
Hmel.Hmel203003oG410	0.062	-0.9224505	0.876	-1.419944	1	1	Hmel.Hmel203003oG410	mod(mdg4)
Hmel.Hmel204001oG307	0.078	-1.1954706	0.844	-5.495291	0.051	0.865		
Hmel.Hmel205001oG359	0.932	1.3017026	0.864	3.98385	0.916	0.187		
Hmel.Hmel205001oG400	0.099	-2.8902303	0.802	-4.590474	0.295	0.92	Hmel.Hmel205001oG400	CG10348
Hmel.Hmel205001oG508	0.084	-1.1779162	0.832	-2.840496	0.378	0.938		
Hmel.Hmel205001oG581	0.914	1.0044731	0.828	1.793487	0.993	0.81	Hmel.Hmel205001oG581	Npc1a
Hmel.Hmel207001oG772	0.973	1.3356199	0.946	2.235331	1	0.896	Hmel.Hmel207001oG772	mahe
Hmel.Hmel207001oG946	0.003	-3.0374094	0.994	-5.371502	0.4	1	Hmel.Hmel207001oG946	Ac76E
Hmel.Hmel208001oG554	0.045	-1.4171786	0.91	#NAME?	0	0.91		
Hmel.Hmel209001oG244	0.911	0.9876925	0.822	1.731742	0.993	0.869	Hmel.Hmel209001oG244	Saf-B
Hmel.Hmel209001oG530	0.948	1.1755896	0.896	2.182911	0.996	0.758	Hmel.Hmel209001oG530	Pus1
Hmel.Hmel209001oG564	0	-3.1661399	1	-6.670378	0.233	1		
Hmel.Hmel210001oG101	0.934	2.3012026	0.868	3.92714	0.945	0.287	Hmel.Hmel210001oG101	CG7110
Hmel.Hmel210001oG285	0.91	1.0477984	0.82	5.828609	0.836	0.035		
Hmel.Hmel210001oG668	0.088	-1.1502987	0.824	#NAME?	0	0.824		
Hmel.Hmel211001oG705	0.086	-1.5853236	0.828	-3.120073	0.473	0.962	Hmel.Hmel211001oG705	pot
Hmel.Hmel212001oG177	0.084	-2.8436209	0.832	-4.422218	0.775	0.993	Hmel.Hmel212001oG177	Fak
Hmel.Hmel212001oG707	0.083	-1.2137492	0.834	-5.866113	0.036	0.851		
Hmel.Hmel213001oG222	0.954	4.1161394	0.908	8.270536	0.927	0.093		
Hmel.Hmel213001oG460	0.91	1.2191246	0.82	2.062541	1	0.872	Hmel.Hmel213001oG460	Nuak
Hmel.Hmel213001oG573	0.087	-2.0242723	0.826	-3.702186	0.258	0.92	Hmel.Hmel213001oG573	CG33110
Hmel.Hmel214004oG273	0.962	1.3914729	0.924	6.601029	0.931	0.031		
Hmel.Hmel215003oG172	0.903	1.0970469	0.806	4.250237	0.855	0.08		
Hmel.Hmel215003oG20	0.094	-1.3164815	0.812	-2.532512	0.436	0.969		
Hmel.Hmel215003oG428	0.935	1.3216086	0.87	6.874205	0.876	0.031		
Hmel.Hmel216002oG203	0.035	-1.5095531	0.93	-5.016177	0.116	0.948		
Hmel.Hmel216002oG468	0.083	-1.200533	0.834	-1.79936	0.876	1		
Hmel.Hmel216002oG538	0.055	-1.5976343	0.89	-2.360826	0.887	1		
Hmel.Hmel216002oG564	0.047	-2.1351878	0.906	-3.590249	0.422	0.997		
Hmel.Hmel217001oG10	0.939	1.1672262	0.878	2.452894	0.975	0.606	Hmel.Hmel217001oG10	Ptp4E
Hmel.Hmel217001oG128	0.056	-1.3631758	0.888	-6.404753	0.036	0.9		
Hmel.Hmel218003oG285	0.913	1.6071798	0.826	4.644738	0.876	0.128	Hmel.Hmel218003oG285	PLCXD
Hmel.Hmel218003oG330	0.073	-1.7475361	0.854	-3.317464	0.309	0.965		
Hmel.Hmel218003oG551	0.954	1.4811211	0.908	7.760008	0.913	0.017		
Hmel.Hmel218003oG779	0.914	1.8894543	0.828	3.085846	0.953	0.232	Hmel.Hmel218003oG779	sv
Hmel.Hmel219001oG225	0.935	1.1642979	0.87	2.469561	0.978	0.571	Hmel.Hmel219001oG225	polybromo
Hmel.Hmel219001oG490	0.93	1.2937833	0.86	9.653132	0.862	0.007		
Hmel.Hmel219001oG548	0.935	1.2734844	0.87	Inf	0.869	0		
Hmel.Hmel220003oG717	0.923	1.1799592	0.846	2.365313	0.982	0.63		
Hmel.Hmel220003oG718	0.901	1.3648486	0.802	4.419596	0.855	0.17		
Hmel.Hmel220003oG746	0.015	-1.7912901	0.97	-3.834154	0.418	0.997		
Hmel.Hmel220003oG747	0.088	-1.4278913	0.824	-2.482992	0.822	0.965	Hmel.Hmel220003oG747	emp
Hmel.Hmel221001oG289	0.953	1.2543227	0.906	1.954869	1	0.965	Hmel.Hmel221001oG289	per
Hmel.Hmel221001oG46	0.956	1.4638259	0.912	3.587488	0.96	0.346	Hmel.Hmel221001oG46	Snmp2
Hmel.Hmel221001oG554	0.021	-1.515494	0.958	-3.002297	0.604	0.997	Hmel.Hmel221001oG554	Neurl4
Hmel.Hmel221001oG591	0.004	-2.0521049	0.992	-5.938904	0.124	0.997	Hmel.Hmel221001oG591	Gk1
Hmel.Hmel221001oG684	0.074	-1.3530059	0.852	-3.793855	0.255	0.907	Hmel.Hmel221001oG684	alpha-Man-Ib
Hmel.Hmel221001oG797	0.987	2.1572645	0.974	8.832505	0.975	0.014		
Hmel.Hmel221001oG86	0.043	-1.4642761	0.914	#NAME?	0	0.913		
Ivory	0.946	2.3261588	0.892	5.083668	0.927	0.27		
								
								
								
Hmel.Hmel201001oG163	RyR
Hmel.Hmel202001oG567	cup
Hmel.Hmel203003oG410	mod(mdg4)
Hmel.Hmel205001oG400	CG10348
Hmel.Hmel205001oG581	Npc1a
Hmel.Hmel207001oG772	mahe
Hmel.Hmel207001oG946	Ac76E
Hmel.Hmel209001oG244	Saf-B
Hmel.Hmel209001oG530	Pus1
Hmel.Hmel210001oG101	CG7110
Hmel.Hmel211001oG705	pot
Hmel.Hmel212001oG177	Fak
Hmel.Hmel213001oG460	Nuak
Hmel.Hmel213001oG573	CG33110
Hmel.Hmel217001oG10	Ptp4E
Hmel.Hmel218003oG285	PLCXD
Hmel.Hmel218003oG779	sv
Hmel.Hmel219001oG225	polybromo
Hmel.Hmel220003oG747	emp
Hmel.Hmel221001oG289	per
Hmel.Hmel221001oG46	Snmp2
Hmel.Hmel221001oG554	Neurl4
Hmel.Hmel221001oG591	Gk1
Hmel.Hmel221001oG684	alpha-Man-Ib


c <- ggplot(markersLateMerged, aes(-log10(p_val_adj),avg_log2FC.y, colour=powerSig))+
   geom_point(alpha=0.4)+
   scale_colour_manual(values = c('black','blue'))+
 #  geom_text_repel(data = filter(test, !is.na(gene_name)), aes(label = gene_name),size = 3,max.overlaps = Inf)+
   theme_classic()+ylim(-10,10)



Idents(Hm.sops.sct) <- 'seurat_clusters'
markersSOPsAll <- FindAllMarkers(Hm.sops.sct,recorrect_umi=FALSE, only.pos=TRUE)

names(new.cluster.ids) <- levels(Hm.sops.sct)
Hm.sops.sct <- RenameIdents(Hm.sops.sct, new.cluster.ids)
Hm.sops.sct$celltypes <- Idents(Hm.sops.sct)


g <- VlnPlot(subset(Hm.sops.sct,race !='iv' & celltypes=='scale'),'PC_1', group.by='time',split.by='race')+ylim(-60,80)
h <- VlnPlot(subset(Hm.sops.sct,race =='iv' & celltypes=='scale'),'PC_1', group.by='time',split.by='mut')+ylim(-60,80)
plot_grid(g,h,nrow=1, rel_widths=c(9,5))




a <- DimPlot(subset(Hm.sops.sct, celltypes != 'sockets'),group.by='seurat_clusters',label=T,repel=T)
b <- FeaturePlot(subset(Hm.sops.sct, celltypes != 'sockets'),'Ivory',min.cutoff=0,order=T)
c <- FeaturePlot(subset(Hm.sops.sct, celltypes != 'sockets'),'Hmel.Hmel218003oG779',min.cutoff=0,order=T)
d <- FeaturePlot(subset(Hm.sops.sct, celltypes != 'sockets'),'Hmel.Hmel212001oG165',min.cutoff=0,order=T)
e <- DimPlot(subset(Hm.sops.sct, celltypes != 'sockets'),group.by='time',shuffle=T)
f <- DimPlot(subset(Hm.sops.sct, celltypes != 'sockets'), group.by='Phase',shuffle=T)
a+b+c+d+e+f

Idents(Hm.sops.sct) <- 'orig.ident'
iv61_vs_Hm25 <- FindMarkers(subset(Hm.sops.sct, celltypes != 'sockets'), ident.1="iv_61", ident.2=c('Hmm_25','Hmr_25'),recorrect_umi=FALSE)
iv61_vs_Hm25$pct.d <- abs(iv61_vs_Hm25$pct.1 - iv61_vs_Hm25$pct.2)
iv61_vs_Hm25[iv61_vs_Hm25$pct.d > 0.75,]

> iv61_vs_Hm25[iv61_vs_Hm25$pct.d > 0.75,]
                              p_val avg_log2FC pct.1 pct.2     p_val_adj pct.d
Hmel.Hmel210001oG964   0.000000e+00  -8.545601 0.189 1.000  0.000000e+00 0.811
Hmel.Hmel220003oG296  1.270779e-298  -4.558170 0.205 0.994 1.974664e-294 0.789
Hmel.Hmel217001oG590  1.695188e-292  -5.480643 0.150 0.967 2.634152e-288 0.817
Hmel.Hmel213001oG738  6.172469e-289  -4.668037 0.225 0.993 9.591399e-285 0.768
Hmel.Hmel201001oG790  1.027586e-286   7.503869 1.000 0.112 1.596766e-282 0.888
Hmel.Hmel210001oG642  1.207435e-285   7.237853 1.000 0.138 1.876233e-281 0.862
Hmel.Hmel215003oG55   4.250098e-285   7.623380 0.998 0.124 6.604228e-281 0.874
Hmel.Hmel219001oG135  2.025600e-284   6.185687 1.000 0.160 3.147580e-280 0.840
Hmel.Hmel201001oG657  1.225782e-283   6.561017 1.000 0.192 1.904743e-279 0.808
Hmel.Hmel211001oG23   9.850796e-283   6.098769 1.000 0.219 1.530715e-278 0.781
Hmel.Hmel201001oG939  3.628821e-282   6.458401 0.995 0.137 5.638825e-278 0.858
Hmel.Hmel206001oG841  4.062485e-282   7.300004 1.000 0.240 6.312696e-278 0.760
Hmel.Hmel201001oG993  2.115043e-281   7.784053 0.983 0.035 3.286566e-277 0.948
Hmel.Hmel219001oG144  1.137231e-280   5.817761 0.993 0.119 1.767144e-276 0.874
Hmel.Hmel203003oG148  4.241711e-280   7.403490 0.984 0.044 6.591195e-276 0.940
Hmel.Hmel207001oG31   2.027640e-279  12.589592 0.974 0.000 3.150750e-275 0.974
Hmel.Hmel219001oG657  5.002236e-278   6.871346 0.980 0.049 7.772974e-274 0.931
Hmel.Hmel215003oG487  1.219074e-277   5.204163 0.995 0.206 1.894320e-273 0.789
Hmel.Hmel217001oG591  9.812864e-273  -5.004320 0.105 0.922 1.524821e-268 0.817
Hmel.Hmel202001oG522  1.072669e-269   4.470255 0.995 0.233 1.666821e-265 0.762
Hmel.Hmel207001oG193  1.568349e-267   7.306725 0.955 0.028 2.437058e-263 0.927
Hmel.Hmel204001oG251  3.943312e-267   8.303025 0.951 0.012 6.127513e-263 0.939
Hmel.Hmel210001oG474  1.058069e-260   6.791056 0.957 0.109 1.644133e-256 0.848
Hmel.Hmel218003oG735  1.229502e-258   4.222376 0.979 0.219 1.910523e-254 0.760
Hmel.Hmel214004oG525  7.431598e-258  -4.582527 0.172 0.924 1.154796e-253 0.752
Hmel.Hmel201001oG802  1.799327e-257   6.124813 0.944 0.067 2.795975e-253 0.877
Hmel.Hmel209001oG79   1.443201e-256   6.751556 0.935 0.033 2.242590e-252 0.902
Hmel.Hmel201001oG1041 2.168439e-256  -5.407656 0.092 0.884 3.369537e-252 0.792
Hmel.Hmel212001oG143  5.579312e-255   6.282034 0.935 0.048 8.669693e-251 0.887
Hmel.Hmel216002oG578  1.034826e-253   5.777823 0.938 0.080 1.608017e-249 0.858
Hmel.Hmel205001oG273  3.031901e-251   5.250106 0.942 0.099 4.711271e-247 0.843
Hmel.Hmel221001oG263  4.512669e-250   6.324729 0.930 0.065 7.012237e-246 0.865
Hmel.Hmel216002oG536  1.031679e-249   5.863177 0.930 0.067 1.603126e-245 0.863
Hmel.Hmel204001oG396  1.353511e-245  -4.313906 0.112 0.887 2.103221e-241 0.775
Hmel.Hmel206001oG635  1.673633e-244   5.680377 0.925 0.062 2.600658e-240 0.863
Hmel.Hmel220003oG241  1.802719e-244   5.874542 0.921 0.057 2.801244e-240 0.864
Hmel.Hmel208001oG546  2.110597e-240   4.515016 0.935 0.176 3.279657e-236 0.759
Hmel.Hmel207001oG42   4.582662e-240   4.097806 0.949 0.192 7.120998e-236 0.757
Hmel.Hmel210001oG940  5.639434e-238  -5.098008 0.073 0.844 8.763116e-234 0.771
Hmel.Hmel206001oG345  1.153090e-233  -5.408698 0.082 0.834 1.791787e-229 0.752
Hmel.Hmel215003oG50   2.982830e-230   5.431761 0.899 0.084 4.635019e-226 0.815
Hmel.Hmel219001oG137  4.863873e-230   5.889875 0.889 0.055 7.557972e-226 0.834
Hmel.Hmel219001oG136  1.639166e-212   6.554712 0.835 0.025 2.547100e-208 0.810
Hmel.Hmel216002oG22   6.054815e-210   4.940375 0.865 0.109 9.408577e-206 0.756
Hmel.Hmel219001oG127  1.210922e-202   7.461738 0.803 0.010 1.881651e-198 0.793
Hmel.Hmel202001oG423  1.291582e-201   6.935612 0.811 0.033 2.006989e-197 0.778
Hmel.Hmel219001oG230  2.892187e-196   5.251084 0.815 0.049 4.494170e-192 0.766



iv61_vs_25_30 <- FindMarkers(subset(Hm.sops.sct, celltypes != 'sockets'), ident.1="iv_61", ident.2=c('Hmm_25','Hmr_25',"Hmr_30","iv_20"),recorrect_umi=FALSE)
iv61_vs_25_30$pct.d <- abs(iv61_vs_25_30$pct.1-iv61_vs_25_30$pct.2)
> head(iv61_vs_25_30,20)
                             p_val avg_log2FC pct.1 pct.2     p_val_adj pct.d#
Hmel.Hmel202001oG423  0.000000e+00  6.7061055 0.811 0.040  0.000000e+00 0.771	# CG18586
Hmel.Hmel221001oG263  0.000000e+00  3.1455092 0.930 0.324  0.000000e+00 0.606	# yin
Hmel.Hmel210001oG572  0.000000e+00 -4.0400805 0.613 1.000  0.000000e+00 0.387	# 
Hmel.Hmel218003oG722  0.000000e+00  2.4799711 0.995 0.609  0.000000e+00 0.386	# CG33978
Hmel.Hmel210001oG623  0.000000e+00 -4.1684238 0.749 0.999  0.000000e+00 0.250	# 
Hmel.Hmel215003oG93   0.000000e+00  3.1511793 0.995 0.833  0.000000e+00 0.162	# 
Hmel.Hmel221001oG385  0.000000e+00  2.1371232 0.999 0.942  0.000000e+00 0.057	# trio
Hmel.Hmel203003oG386  0.000000e+00  1.3681748 1.000 0.990  0.000000e+00 0.010	# CG7800
Hmel.Hmel200185oG1    0.000000e+00  0.9917599 1.000 1.000  0.000000e+00 0.000	# 
Hmel.Hmel219001oG17  1.732284e-303  3.5617531 0.853 0.193 2.691796e-299 0.660	# CG11498
Hmel.Hmel206001oG841 8.516100e-298  1.9171998 1.000 0.572 1.323317e-293 0.428	# 
Hmel.Hmel219001oG127 1.387320e-297  3.8084505 0.803 0.124 2.155756e-293 0.679	# Osi7
Hmel.Hmel210001oG474 4.432012e-291  2.9313528 0.957 0.406 6.886903e-287 0.551	# exp
Hmel.Hmel220003oG342 1.121602e-289 -4.4008352 0.357 0.936 1.742857e-285 0.579	# tnc
Hmel.Hmel207001oG518 6.503238e-289  1.5936717 1.000 0.980 1.010538e-284 0.020	# Svil
Hmel.Hmel217001oG182 2.622883e-282  2.0605487 0.998 0.623 4.075698e-278 0.375	# CG16786
Hmel.Hmel210001oG711 5.567680e-282 -2.3436157 0.701 0.987 8.651618e-278 0.286	# 
Hmel.Hmel219001oG15  9.801249e-281  2.9892729 0.874 0.300 1.523016e-276 0.574	# Rsod
Hmel.Hmel215003oG101 4.310719e-267  3.9433540 0.770 0.162 6.698426e-263 0.608	# Dnaaf5
Hmel.Hmel221001oG563 7.581751e-265 -7.3334872 0.031 0.738 1.178128e-260 0.707	# yin

VlnPlot(subset(Hm.sops.sct,orig.ident %in% c('Hmm_25','Hmr_25',"Hmr_30","iv_20", 'iv_61')), a,group.by='mut')

Idents(Hm.sops.sct) <- 'race'
iv_vs_malros <- FindMarkers(subset(Hm.sops.sct, celltypes != 'sockets'), ident.1="iv", ident.2=c('mal','ros'),recorrect_umi=FALSE)
iv_vs_malros$pct.d <- abs(iv_vs_malros$pct.1-iv_vs_malros$pct.2)
iv_vs_malros[iv_vs_malros$pct.d > 0.6,]
                     p_val avg_log2FC pct.1 pct.2 p_val_adj pct.d
Hmel.Hmel201001oG993     0   2.063879 0.967 0.214         0 0.753
Hmel.Hmel204001oG251     0   1.105943 0.932 0.208         0 0.724
Hmel.Hmel207001oG31      0   1.734725 0.953 0.254         0 0.699
Hmel.Hmel221001oG263     0   4.173591 0.896 0.201         0 0.695
Hmel.Hmel219001oG135     0   2.442260 0.981 0.297         0 0.684
Hmel.Hmel210001oG642     0   2.492070 0.977 0.306         0 0.671
Hmel.Hmel219001oG230     0   2.496129 0.792 0.128         0 0.664
Hmel.Hmel215003oG487     0   2.006210 0.984 0.331         0 0.653
Hmel.Hmel220003oG241     0   1.051925 0.962 0.317         0 0.645
Hmel.Hmel216002oG312     0   2.457429 0.986 0.348         0 0.638
Hmel.Hmel202001oG522     0   1.540297 0.990 0.353         0 0.637
Hmel.Hmel209001oG79      0   2.303082 0.851 0.215         0 0.636
Hmel.Hmel201001oG657     0   1.837370 0.949 0.320         0 0.629
Hmel.Hmel212001oG819     0   2.845365 0.840 0.213         0 0.627
Hmel.Hmel211001oG23      0   2.924945 0.995 0.372         0 0.623
Hmel.Hmel219001oG144     0   1.378070 0.980 0.361         0 0.619
Hmel.Hmel219001oG136     0   3.145065 0.797 0.181         0 0.616
Hmel.Hmel212001oG143     0   2.454984 0.871 0.261         0 0.610
Hmel.Hmel215003oG50      0   3.688286 0.731 0.122         0 0.609


Idents(Hm.sops.sct) <- 'orig.ident'
hm30_vs_iv20 <- FindMarkers(subset(Hm.sops.sct, celltypes != 'sockets'), ident.1="Hmr_30", ident.2='iv_20',recorrect_umi=FALSE)

Idents(Hm.sops.sct) <- 'mut'
scales_wtVmut <- FindMarkers(subset(Hm.sops.sct, celltypes == 'scale' & orig.ident !='iv_33' & time > 21), ident.1 ='wt', ident.2='mut', recorrect_umi=FALSE)
scales_wtVmut$pct.d <- abs(scales_wtVmut$pct.1-scales_wtVmut$pct.2)
scales_wtVmut[scales_wtVmut$pct.d > 0.5,]
                             p_val avg_log2FC pct.1 pct.2     p_val_adj pct.d
Hmel.Hmel221001oG563  0.000000e+00   7.485974 0.698 0.024  0.000000e+00 0.674
Hmel.Hmel215003oG101  0.000000e+00  -4.270087 0.156 0.787  0.000000e+00 0.631
Hmel.Hmel215003oG104  0.000000e+00  -4.111967 0.172 0.751  0.000000e+00 0.579
Hmel.Hmel215003oG92   0.000000e+00  -4.070336 0.231 0.783  0.000000e+00 0.552
Hmel.Hmel200219oG4    0.000000e+00   2.717029 0.855 0.308  0.000000e+00 0.547
Hmel.Hmel219001oG17   0.000000e+00  -2.656100 0.270 0.817  0.000000e+00 0.547
Hmel.Hmel202001oG423  0.000000e+00  -6.030705 0.046 0.581  0.000000e+00 0.535
Hmel.Hmel215003oG103  0.000000e+00  -3.373463 0.318 0.850  0.000000e+00 0.532
Hmel.Hmel221001oG263  0.000000e+00  -3.102772 0.424 0.947  0.000000e+00 0.523
Hmel.Hmel215003oG90   0.000000e+00  -2.181176 0.679 0.926  0.000000e+00 0.247
Hmel.Hmel213001oG244  0.000000e+00  -1.713329 0.756 0.990  0.000000e+00 0.234
Hmel.Hmel215003oG93   0.000000e+00  -3.880415 0.866 0.995  0.000000e+00 0.129
Hmel.Hmel215003oG438 9.536345e-280   5.353596 0.621 0.082 1.481853e-275 0.539
Hmel.Hmel208001oG189 4.698173e-274   1.371450 0.993 0.958 7.300491e-270 0.035
Hmel.Hmel217001oG932 3.284321e-270  -1.879845 0.491 0.864 5.103506e-266 0.373
Hmel.Hmel204001oG224 4.979199e-269  -1.153734 0.900 0.996 7.737178e-265 0.096
Ivory                1.157757e-268   4.218565 0.946 0.894 1.799038e-264 0.052
Hmel.Hmel219001oG15  3.953691e-266  -2.181440 0.409 0.835 6.143640e-262 0.426
HMEL000025           3.922545e-262   2.380173 0.960 0.948 6.095242e-258 0.012
Hmel.Hmel219001oG138 3.620328e-248  -1.907531 0.438 0.876 5.625628e-244 0.438
VlnPlot(subset(Hm.sops.sct, orig.ident != 'iv_33'), a,group.by='mut')

iv55_vs_Hm40 <- FindMarkers(subset(Hm.sops.sct, celltypes != 'socket'), ident.1="iv_55", ident.2=c('Hmm_40','Hmr_40'),recorrect_umi=FALSE)


####################################################################
##### realnalysing all non sops, to see if time slipping holds #####
####################################################################

Hm_merged$celltypes <- Hm_merged$sop
Idents(Hm_merged) <- Hm_merged$celltypes

'%ni%' <- function(x,y)!('%in%'(x,y))
Hm_Notsops <- subset(Hm_merged, celltypes %ni% c('sop','socket','scale'))
DimPlot(Hm_Notsops)


DefaultAssay(Hm_Notsops) <- 'integrated'

Hm_Notsops <- FindVariableFeatures(Hm_Notsops)
#features <- FindVariableFeatures(Hm_Notsops)
Hm_Notsops <- ScaleData(Hm_Notsops, verbose = FALSE)

Hm_Notsops <- RunPCA(Hm_Notsops, npcs=50, verbose = FALSE)

##based on elbow plot i'm goign to go for 20 dims
DefaultAssay(Hm_Notsops) <- 'integrated'
Hm_Notsops <- RunUMAP(Hm_Notsops, reduction = "pca", dims = 1:8)
Hm_Notsops <- FindNeighbors(Hm_Notsops, reduction = "pca", dims = 1:8)
Hm_Notsops <- FindClusters(Hm_Notsops, resolution = 2)

DefaultAssay(Hm_Notsops) <- 'SCT'
a <- DimPlot(Hm_Notsops,group.by='time')
b <- FeaturePlot(Hm_Notsops,'Ivory')
c <- FeaturePlot(Hm_Notsops,'Hmel.Hmel218003oG779')
d <- FeaturePlot(Hm_Notsops,'Hmel.Hmel207001oG974')
a+b+c+d

### need to reupload to the cluster to run this bit
##SeuratDisk::SaveH5Seurat(Hm_sops, 'Hm_sops', overwrite = TRUE, verbose = TRUE)
##Hm_sops <- LoadH5Seurat('Hm_sops*')


# it worked really well with 5 dims and then running 1:3 locally, but trying again with 1:4 here
DefaultAssay(Hm_sops) <- 'RNA'
Hm_Notsops <- SCTransform(Hm_Notsops, verbose = T, variable.features.n = 5000)
Hm_Notsops <- RunPCA(Hm_Notsops,  npcs = 50, verbose = T)
Hm_Notsops <- RunUMAP(Hm_Notsops,  dims = 1:5)
Hm_Notsops <- FindNeighbors(Hm_Notsops, dims = 1:5)
Hm_Notsops <- FindClusters(Hm_Notsops, resolution = 2)

#Hm_sops$race <- ifelse(Hm_sops$orig.ident == "Hmm_10" | Hm_merged$orig.ident == "Hmm_25" | Hm_merged$orig.ident == 'Hmm_40', "mal", "ros")


Hm_list <- SplitObject(Hm_Notsops, split.by = "race")
#######Hm_list <- lapply(X = Hm_list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = Hm_list, nfeatures = 7000)
Hm_list <- PrepSCTIntegration(object.list = Hm_list, anchor.features = features)

Hm.anchors <- FindIntegrationAnchors(object.list = Hm_list, normalization.method = "SCT",
    anchor.features = features)
Hm.Notsops.sct <- IntegrateData(anchorset = Hm.anchors, normalization.method = "SCT")

#SaveH5Seurat(Hm.sops.sct, 'Hm.sops.sct', overwrite = TRUE, verbose = TRUE)
#Hm.sops.sct <- LoadH5Seurat('Hm.sops.sct.h5seurat')

### so far, 1:3 works best
DefaultAssay(Hm_sops) <- 'integrated'
Hm.Notsops.sct <- RunPCA(Hm.Notsops.sct,  npcs = 50, verbose = F)
Hm.Notsops.sct <- RunUMAP(Hm.Notsops.sct,min.dist= 0.5,  dims = c(1:3))
Hm.Notsops.sct <- FindNeighbors(Hm.Notsops.sct, dims =  c(1:3))
Hm.Notsops.sct <- FindClusters(Hm.Notsops.sct, resolution =1)
#Hm.Notsops.sct  <- NormalizeData(Hm.Notsops.sct , normalization.method = "LogNormalize", assay = 'RNA')

DefaultAssay(Hm.Notsops.sct) <- 'RNA'
a <- DimPlot(Hm.Notsops.sct,group.by='celltypes',label=T,repel=T)
b <- DimPlot(Hm.Notsops.sct,group.by='race',order=T)
c <- FeaturePlot(Hm.Notsops.sct,'Hmel.Hmel218003oG779',min.cutoff=0,order=T) # sv
d <- FeaturePlot(Hm.Notsops.sct,'Hmel.Hmel212001oG165',min.cutoff=0,order=T) # sox15
e <- DimPlot(Hm.Notsops.sct,group.by='time',shuffle=T)
f <- DimPlot(Hm.Notsops.sct, group.by='orig.ident',shuffle=T,label=T)
a+b+c+d+e+f

PCAPlot(Hm.Notsops.sct,group.by='time')

### label wts and muts
barcodes <- subset(x = Hm.Notsops.sct, orig.ident %in% c('iv_61', 'iv_38','iv_55'))
Hm.Notsops.sct$mut <- ""
Hm.Notsops.sct$mut[WhichCells(Hm.Notsops.sct) %in%  WhichCells(barcodes)] <- "mut"
Hm.Notsops.sct$mut[!(WhichCells(Hm.Notsops.sct) %in%  WhichCells(barcodes))] <- "wt"

### make PC1 'pseudotime' comparison plot
i <- VlnPlot(subset(Hm.Notsops.sct,race !='iv'),'PC_1', group.by='time',split.by='mut')+ylim(-30,90)
j <- VlnPlot(subset(Hm.Notsops.sct,race =='iv' & time < 50),'PC_1', group.by='time',split.by='mut')+ylim(-30,90)
plot_grid(i,j,nrow=1, rel_widths=c(7,4))

### do wt vs mut in non-sops
Idents(Hm.Notsops.sct) <- Hm.Notsops.sct$mut
notSop.markers <- FindMarkers(subset(Hm.Notsops.sct, time > 20 & time < 40), ident.1="wt", ident.2="mut",recorrect_umi=FALSE)
notSop.markers$pct.d <- abs(notSop.markers$pct.1 - notSop.markers$pct.2)

Idents(Hm.sops.sct) <- Hm.sops.sct$mut
Sop.markers <- FindMarkers(subset(Hm.sops.sct,time > 20 & time < 40), ident.1="wt", ident.2="mut",recorrect_umi=FALSE)
Sop.markers$pct.d <- abs(Sop.markers$pct.1 - Sop.markers$pct.2)

### trying to compare this to the differences between ros and mel
test <- subset(Hm.sops.sct, orig.ident %in% c('iv_20','iv_61','Hmm_25','Hmr_25','Hmr_30') & celltypes=='scale')
barcodes <- subset(x = test, Ivory > 1 & race=='iv')
barcodes2<- subset(x = test, Ivory > 4 & race!='iv')
test$Iv <- ""
test$Iv[WhichCells(test) %in%  WhichCells(barcodes)] <- "hi"
test$Iv[WhichCells(test) %in%  WhichCells(barcodes2)] <- 'hi'
test$Iv[test$Iv == ''] <- 'lo'
#test$Iv[WhichCells(test) %ni%  WhichCells(barcodes)] <- "lo"
#test$Iv[WhichCells(test) %ni%  WhichCells(barcodes2)] <- 'lo'



Idents(test) <- test$mut
WtVMutInHi.Markers <- FindMarkers(test, ident.1='wt',recorrect_umi=FALSE)
genes <- c("Hmel.Hmel201001oG551","Hmel.Hmel201001oG657","Hmel.Hmel201001oG790","Hmel.Hmel201001oG802","Hmel.Hmel201001oG993","Hmel.Hmel202001oG423","Hmel.Hmel203003oG153","Hmel.Hmel203003oG386","Hmel.Hmel205001oG273","Hmel.Hmel206001oG635","Hmel.Hmel206001oG777","Hmel.Hmel206001oG841","Hmel.Hmel207001oG518","Hmel.Hmel210001oG474","Hmel.Hmel210001oG572","Hmel.Hmel210001oG623","Hmel.Hmel210001oG711","Hmel.Hmel210001oG964","Hmel.Hmel211001oG624","Hmel.Hmel213001oG244","Hmel.Hmel213001oG68","Hmel.Hmel214004oG555","Hmel.Hmel215003oG101","Hmel.Hmel215003oG50","Hmel.Hmel215003oG55","Hmel.Hmel215003oG93","Hmel.Hmel216002oG269","Hmel.Hmel217001oG182","Hmel.Hmel218003oG36","Hmel.Hmel218003oG722","Hmel.Hmel219001oG125","Hmel.Hmel219001oG127","Hmel.Hmel219001oG15","Hmel.Hmel219001oG17","Hmel.Hmel219001oG657","Hmel.Hmel220003oG342","Hmel.Hmel221001oG263","Hmel.Hmel221001oG385","Hmel.Hmel221001oG563")
names <- c("JhI-21","Osi21","CG6055","Pka-C3","mfas","CG18586","m","CG7800","CG3961","LManI","Invadolysin","Hmel.Hmel206001oG841","Svil","exp","nw","CG1273","Hmel.Hmel210001oG711","NLaz","Gasp","yellow-c","Lmpt","dpr21","Dnaaf5","Hmel.Hmel215003oG50 ","Cht7","Hmel.Hmel215003oG93 ","CG13377","CG16786","CG42788","CG33978","Osi4","Osi7","Rsod","CG11498","Hmel.Hmel219001oG657","tnc","yin","trio","Hmel.Hmel221001oG563")
a <- DotPlot(test,features=genes,scale=F)+coord_flip()+ theme(axis.text.x = element_text(angle = 30, hjust=1))+scale_color_batlow(reverse=T)+scale_x_discrete(labels=names)
tab <- data.frame(genes,names)
WtVMutInHi.Markers$genes <- row.names(WtVMutInHi.Markers)
WtVMutInHi.Markers <- merge(WtVMutInHi.Markers,tab,by='genes',all=TRUE)
WtVMutInHi.Markers <- WtVMutInHi.Markers[order(WtVMutInHi.Markers$p_val_adj), ]

a2 <- ggplot(WtVMutInHi.Markers, aes(avg_log2FC, -log10(p_val_adj)))+
             geom_point()+
             geom_text_repel(data=subset(WtVMutInHi.Markers, avg_log2FC < 0),aes(label=names),force=10,max.overlaps=Inf,xlim=c(-8,-2))+
             geom_text_repel(data=subset(WtVMutInHi.Markers, avg_log2FC > 0),aes(label=names),force=10,max.overlaps=Inf,xlim=c(2,10))+
             theme_classic()+coord_flip()


Idents(test) <- test$Iv
HiVsLo.markers <- FindMarkers(test, ident.1="hi", ident.2="lo",recorrect_umi=FALSE)
genes <- c("Hmel.Hmel201001oG540","Hmel.Hmel201001oG768","Hmel.Hmel201001oG923","Hmel.Hmel201001oG924","Hmel.Hmel203003oG329","Hmel.Hmel203003oG331","Hmel.Hmel205001oG243","Hmel.Hmel206001oG108","Hmel.Hmel206001oG345","Hmel.Hmel206001oG497","Hmel.Hmel206001oG949","Hmel.Hmel207001oG801","Hmel.Hmel208001oG159","Hmel.Hmel209001oG214","Hmel.Hmel209001oG285","Hmel.Hmel209001oG81","Hmel.Hmel210001oG1051","Hmel.Hmel210001oG276","Hmel.Hmel211001oG153","Hmel.Hmel211001oG557","Hmel.Hmel212001oG989","Hmel.Hmel213001oG153","Hmel.Hmel214004oG334","Hmel.Hmel215003oG461","Hmel.Hmel215003oG87","Hmel.Hmel215003oG90","Hmel.Hmel216002oG72","Hmel.Hmel217001oG40","Hmel.Hmel217001oG594","Hmel.Hmel218003oG49","Hmel.Hmel218003oG749","Hmel.Hmel218003oG779","Hmel.Hmel219001oG214","Hmel.Hmel219001oG215","Hmel.Hmel219001oG281","Hmel.Hmel220003oG799","Hmel.Hmel220003oG90","Hmel.Hmel221001oG246","HMEL000025","Ivory")
names <- c("mid","Nha1","UNK_next_to_al1","al1","UNK_next_to_Antp","Antp","Hmel.Hmel205001oG243","mamo","CG13954","Hmel.Hmel206001oG497","Hmel.Hmel206001oG949","Src64B","dsx","Spn100A","shd","ss","conu","tup","dally","B-H1","jhamt","grn","sha","ckn","Hmel.Hmel215003oG87","Hmel.Hmel215003oG90","Hmel.Hmel216002oG72","CG43980","Hmel.Hmel217001oG594","Optix","rn","sv","Dr","Dr_split_annot)","CG11147","CG42741","CG7991","Hmel.Hmel221001oG246","HMEL000025","Ivory")
tab <- data.frame(genes,names)
HiVsLo.markers$genes <- row.names(HiVsLo.markers)
HiVsLo.markers <- merge(HiVsLo.markers,tab,by='genes',all=TRUE)
HiVsLo.markers <- HiVsLo.markers[order(HiVsLo.markers$p_val_adj), ]
b <- DotPlot(test,features=genes,scale=FALSE,group.by='Iv')+coord_flip()+ theme(axis.text.x = element_text(angle = 30, hjust=1))+scale_color_batlow(reverse=T)+scale_x_discrete(labels=names)
b1 <- ggplot(HiVsLo.markers, aes(avg_log2FC, -log10(p_val_adj),label=names))+
        #   scale_y_break(c(50,100))+
           geom_point()+
           geom_text_repel(data=subset(HiVsLo.markers, avg_log2FC < 0),aes(label = names),force=10,max.overlaps=Inf,xlim=c(-4,-2.5),ylim=c(1,Inf),nudge_x=-1)+
           geom_text_repel(data=subset(HiVsLo.markers, avg_log2FC > 0),aes(label = names),force=10,max.overlaps=Inf,xlim=c(.5,4),ylim=c(1,Inf),nudge_x=1)+
           theme_classic()+coord_flip()

### running go term analysis
library("BiocManager")

library('clusterProfiler')
library('AnnotationDbi')
library('enrichplot')
library('ggnewscale')
library('cowplot')
library('fly.db0')
library('org.Dm.eg.db')

#following this code put together by ling
x <- org.Dm.egFLYBASEPROT # Get the entrez gene IDs that are mapped to a Flybase prot ID 
mapped_genes <- mappedkeys(x) # Convert to a list 
# Get the FlyBase protein IDs for the first five genes xx[1:5] # Get the first one 
xx <- as.list(x[mapped_genes])
 if(length(xx) > 0){
    xx[[1]] 
    }







Idents(Hm.sops.sct) <- Hm.sops.sct$Iv
HiVsLoHmm25.markers <- FindMarkers(subset(Hm.sops.sct, orig.ident %in% c("Hmm_25")), ident.1="hi", ident.2="lo",recorrect_umi=FALSE)
HiVsLoHmm40.markers <- FindMarkers(subset(Hm.sops.sct, orig.ident %in% c("Hmm_40")), ident.1="hi", ident.2="lo",recorrect_umi=FALSE)
HiVsLoHmm.markers <- FindMarkers(subset(Hm.sops.sct, orig.ident %in% c("Hmm_40","Hmm_25")), ident.1="hi", ident.2="lo",recorrect_umi=FALSE)


Idents(Hm.sops.sct) <- Hm.sops.sct$orig.ident
Hmr25_vs_iv63.markers <- FindMarkers(Hm.sops.sct, ident.1="Hmr_25", ident.2="iv_61",recorrect_umi=FALSE)

#### use ivory and optix to define things into four classes
Hm.sops.sct$scaletype <- ""
barcodes <- subset(x = Hm.sops.sct, Ivory > 0.5 & Hmel.Hmel218003oG49 > 0.1)
Hm.sops.sct$scaletype[WhichCells(Hm.sops.sct) %in%  WhichCells(barcodes)] <- "TypeIII"
barcodes <- subset(x = Hm.sops.sct, Ivory < 0.5 & Hmel.Hmel218003oG49 < 0.1)
Hm.sops.sct$scaletype[WhichCells(Hm.sops.sct) %in%  WhichCells(barcodes)] <- "TypeI"
barcodes <- subset(x = Hm.sops.sct, Ivory > 0.5 & Hmel.Hmel218003oG49 < 0.1)
Hm.sops.sct$scaletype[WhichCells(Hm.sops.sct) %in%  WhichCells(barcodes)] <- "TypeII"
barcodes <- subset(x = Hm.sops.sct, Ivory < 0.5 & Hmel.Hmel218003oG49 > 0.1)
Hm.sops.sct$scaletype[WhichCells(Hm.sops.sct) %in%  WhichCells(barcodes)] <- "Overlap"
## force all IV scales to be type 1
barcodes <- subset(x = Hm.sops.sct, mut == "mut" & scaletype %in% c("TypeII","TypeIII"))
Hm.sops.sct$scaletype[WhichCells(Hm.sops.sct) %in%  WhichCells(barcodes)] <- "TypeI"
## force all rosina scales that get called type I above to be called type II
barcodes <- subset(x = Hm.sops.sct, race == 'ros' & scaletype == 'TypeI')
Hm.sops.sct$scaletype[WhichCells(Hm.sops.sct) %in%  WhichCells(barcodes)] <- "TypeII"


Idents(Hm.sops.sct) <- Hm.sops.sct$scaletype
TypeIII.Early.markers <- FindMarkers(subset(Hm.sops.sct,celltypes != 'sockets' & time < 25), ident.1="TypeIII",recorrect_umi=FALSE)

Type1vs2 <- FindMarkers(subset(Hm.sops.sct,celltypes != 'sockets'), ident.1="TypeI", ident.2='TypeII',recorrect_umi=FALSE)
Type1vs2$pct.d <- abs(Type1vs2$pct.1 - Type1vs2$pct.2)

VlnPlot(subset(Hm.sops.sct, race=='ros'), 'Ivory', group.by='scaletypes')

#### use drop to label D vs V
Hm.sops.sct$DV <- ""
barcodes <- subset(x = Hm.sops.sct, Hmel.Hmel219001oG214 > 0.5)
Hm.sops.sct$DV[WhichCells(Hm.sops.sct) %in%  WhichCells(barcodes)] <- "Dorsal"
barcodes <- subset(x = Hm.sops.sct, Hmel.Hmel219001oG214 < 0.5 )
Hm.sops.sct$DV[WhichCells(Hm.sops.sct) %in%  WhichCells(barcodes)] <- "Ventral"
## force all overlap scales to be ventral
barcodes <- subset(x = Hm.sops.sct, scaletype == 'Overlap' )
Hm.sops.sct$DV[WhichCells(Hm.sops.sct) %in%  WhichCells(barcodes)] <- "Ventral"

Idents(Hm.sops.sct) <- Hm.sops.sct$DV
DV.markers <- FindMarkers(Hm.sops.sct, ident.1="Dorsal",recorrect_umi=FALSE)
Hmel.Hmel219001oG214  0.000000e+00 14.1206353 1.000 0.000  0.000000e+00	Dr
Hmel.Hmel219001oG215  0.000000e+00  2.2622439 0.583 0.152  0.000000e+00	
Hmel.Hmel206001oG108 1.197673e-289  2.6621943 0.619 0.262 1.788485e-285	mamo
Hmel.Hmel220003oG799 5.023422e-160  1.4430844 0.624 0.376 7.501476e-156	CG42741
Hmel.Hmel202001oG493 1.046733e-148 -0.9150466 0.862 0.937 1.563087e-144	bi
HMEL000025           5.827169e-143  1.0145828 0.933 0.811 8.701712e-139	cortex
Ivory                1.130451e-130  1.1086654 0.908 0.757 1.688102e-126	ivory
Hmel.Hmel201001oG873 8.265083e-110  1.2142193 0.840 0.729 1.234225e-105	CG5708
Hmel.Hmel211001oG557 9.136663e-107 -1.5606486 0.193 0.413 1.364378e-102	BH-1
Hmel.Hmel213001oG451  2.900928e-77  0.5108780 0.999 0.995  4.331956e-73	heph

DV_Topmarkers <- c('Hmel.Hmel219001oG214',
'Hmel.Hmel206001oG108',
'Hmel.Hmel220003oG799',
'Hmel.Hmel202001oG493',
'HMEL000025',
'Ivory',
'Hmel.Hmel201001oG873',
'Hmel.Hmel211001oG557',
'Hmel.Hmel213001oG451')

labels <- c('Dr'
'mamo'
'CG4274'
'bi'
'cortex'
'ivory'
'CG5708'
'BH-1'
'heph')

#### use dve to label ground vs cover
Hm.sops.sct$ground <- ""
barcodes <- subset(x = Hm.sops.sct, Hmel.Hmel210001oG901 < 0.1)
Hm.sops.sct$ground[WhichCells(Hm.sops.sct) %in%  WhichCells(barcodes)] <- "cover"
Hm.sops.sct$ground[WhichCells(Hm.sops.sct) %ni% WhichCells(barcodes)] <- "ground"
Idents(Hm.sops.sct) <- Hm.sops.sct$ground
ground.markers <- FindMarkers(Hm.sops.sct, ident.1="ground",recorrect_umi=FALSE)

groundHet.markers <- FindMarkers(subset(Hm.sops.sct,orig.ident=='iv_33'), ident.1="ground",recorrect_umi=FALSE)



#### use cut to label hairs
Hm.sops.sct$hair <- ""
barcodes <- subset(x = Hm.sops.sct, Hmel.Hmel202001oG51 > 1)
Hm.sops.sct$hair[WhichCells(Hm.sops.sct) %in%  WhichCells(barcodes)] <- "hair"
Hm.sops.sct$hair[WhichCells(Hm.sops.sct) %ni% WhichCells(barcodes)] <- "NotHair"
Idents(Hm.sops.sct) <- Hm.sops.sct$hair
hair.markers <- FindMarkers(Hm.sops.sct, ident.1="hair",recorrect_umi=FALSE)

#### trying mut v wt, but ignoring wt type i
Idents(Hm.sops.sct) <- Hm.sops.sct$mut
Hm.sops.sct$wtTI <- ""
barcodes <- subset(x = Hm.sops.sct, scaletype == 'TypeI' & mut == 'wt')
Hm.sops.sct$wtTI[WhichCells(Hm.sops.sct) %in%  WhichCells(barcodes)] <- "wt"
Hm.sops.sct$wtTI[WhichCells(Hm.sops.sct) %ni% WhichCells(barcodes)] <- "other"

mut.markers <- FindMarkers(subset(Hm.sops.sct, wtTI != 'wt'), ident.1="mut",recorrect_umi=FALSE)


#### mut v wt but just in notsops
Idents(Hm.Notsops.sct) <- Hm.Notsops.sct$mut

mut.markers <- FindMarkers(subset(Hm.sops.sct), ident.1="mut",recorrect_umi=FALSE)


#### trying "iv_63" vs c("Hmm_25","Hmr_30","iv_20"), followed by iv_20 vs het
Idents(Hm.sops.sct) <- Hm.sops.sct$orig.ident

mut.markers <- FindMarkers(subset(Hm.sops.sct), ident.1="iv_61", ident.2=c("Hmm_25","Hmr_30","iv_20"),recorrect_umi=FALSE)
mut.markers$pct.d.Wt <- abs(mut.markers$pct.1-mut.markers$pct.2)
mutVhet.markers <- FindMarkers(subset(Hm.sops.sct), ident.1="iv_61", ident.2=c("iv_33"),recorrect_umi=FALSE)
mutVhet.markers$pct.d.het <- abs(mutVhet.markers$pct.1-mutVhet.markers$pct.2)

mut.markers$gene <- rownames(mut.markers)
mutVhet.markers$gene <- rownames(mutVhet.markers)

a <- merge(mut.markers,mutVhet.markers,by="gene")

####scaletypes 
Idents(Hm.sops.sct) <- Hm.sops.sct$scaletype
T1.markers <- FindMarkers(subset(Hm.sops.sct, orig.ident %in% c('Hmm_25','Hmm_40','iv_20')), ident.1="TypeI",recorrect_umi=FALSE)

DefaultAssay(Hm_merged) <- 'SCT'
## looking for ivory correlation in the het

# List of sample IDs to iterate over
sample_ids <- c("iv_20","iv_61","iv_33", "iv_55", "iv_38", "Hmm_10","Hmr_10","Hmr_15","Hmr_20","Hmr_25","Hmm_25","Hmr_30","Hmr_35","Hmr_40","Hmm_40","Hmr_45","Hmr_50") 
# Initialize a named list to hold the results
IvoryCorrList <- list()
# Loop through each sample
for (sample_id in sample_ids) {
  
  # Subset Seurat object and get data
  matrix <- GetAssayData(object = subset(Hm_merged, orig.ident %in% sample_id),
                         assay = "SCT",
                         slot = "data")
  matrix_mod <- as.matrix(matrix)
  
  # Check if "Ivory" exists in filtered matrix
  if (!"Ivory" %in% rownames(matrix_mod)) {
    warning(paste("Ivory gene not found in", sample_id, "- skipping"))
    next
  }
  
  # Filter genes by expression
  cell_count <- ncol(matrix_mod)
  min_cells <- ceiling(0.001 * cell_count)
  min_expr <- 0
  gene_filter <- rowSums(matrix_mod > 0) >= min_cells & rowMeans(matrix_mod) > min_expr
  matrix_mod <- matrix_mod[gene_filter, ]
  
  # Double-check again after filtering
  if (!"Ivory" %in% rownames(matrix_mod)) {
    warning(paste("Ivory gene filtered out in", sample_id, "- skipping"))
    next
  }
  
  # Correlation and expression
  gene <- as.numeric(matrix_mod["Ivory", ])
  correlations <- apply(matrix_mod, 1, function(x) cor(gene, x))
  mean_expression <- rowMeans(matrix_mod)
  
  # Assemble result
  result_df <- data.frame(
    gene = names(correlations),
    correlation = correlations,
    mean_expression = mean_expression[names(correlations)]
  )
  result_df <- result_df[order(-result_df$correlation), ]
  
  # Store in named list
  IvoryCorrList[[sample_id]] <- result_df
}

#IvoryCorrList_filtered <- IvoryCorrList[!sapply(IvoryCorrList, is.null)]

# Now merge
IvoryCorrWide <- Reduce(function(x, y) full_join(x, y, by = "gene"), 
  lapply(names(IvoryCorrList), function(sample_id) {
    df <- IvoryCorrList[[sample_id]]
    df <- df %>%
      select(gene, correlation, mean_expression) %>%
      rename(!!paste0("cor_", sample_id) := correlation,
             !!paste0("expr_", sample_id) := mean_expression)
    return(df)
  })
)


plot(IvoryCorrWide$expr_iv_61, IvoryCorrWide$expr_Hmr_30)
test <- IvoryCorrWide[IvoryCorrWide$expr_iv_61 > 2 & IvoryCorrWide$expr_Hmr_30 < 1,]
test <- test[,c(1,4,5,24,25)]
test <- na.omit(test)


### MERGING ALL THESE CORRELATION MATRICES
# Get all unique gene names
all_genes <- IvoryCorrWide$gene
# Build data frame
df <- data.frame(
  Gene = all_genes,
  iv61  = as.numeric(IvoryCorrWide$cor_iv_61),
  Hm20  = as.numeric(IvoryCorrWide$cor_Hmr_20),
  Hm35 = as.numeric(IvoryCorrWide$cor_Hmr_35),
  iv20  = as.numeric(IvoryCorrWide$cor_iv_20),
  Hm10  = as.numeric(IvoryCorrWide$cor_Hmr_10),
  Hm25  = as.numeric(IvoryCorrWide$cor_Hmr_25),
  Hm40  = as.numeric(IvoryCorrWide$cor_Hmr_40),
  iv33  = as.numeric(IvoryCorrWide$cor_iv_33),
  Hm15  = as.numeric(IvoryCorrWide$cor_Hmr_15),
  Hm30  = as.numeric(IvoryCorrWide$cor_Hmr_30),
  iv38  = as.numeric(IvoryCorrWide$cor_iv_38),
  iv55  = as.numeric(IvoryCorrWide$cor_iv_55),
  Hm45  = as.numeric(IvoryCorrWide$cor_Hmr_45),
  Hmm10  = as.numeric(IvoryCorrWide$cor_Hmm_10),
  Hmm25  = as.numeric(IvoryCorrWide$cor_Hmm_25),
  Hmm40  = as.numeric(IvoryCorrWide$cor_Hmm_40),
  Hm50  = as.numeric(IvoryCorrWide$cor_Hmr_50)
)

a <- ggplot(IvoryCorrWide, aes(cor_Hmr_30,expr_Hmr_30))+
   geom_point()+
   xlim(-1,1)+
   theme_classic()
b <- ggplot(IvoryCorrWide, aes(cor_iv_61,expr_iv_61))+
   geom_point()+
   xlim(-1,1)+
   theme_classic()
plot_grid(a,b,ncol=1)

a <- ggplot(df,aes(iv61,iv20))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()
    
b <- ggplot(df,aes(iv61,Hm25))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()

c <- ggplot(df,aes(iv61,Hm30))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()

d <- ggplot(df,aes(iv61,Hm35))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()
    
e <- ggplot(df,aes(Hm25,iv61))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()
    
f <- ggplot(df,aes(Hm35,iv55))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()

g <- ggplot(df,aes(Hm40,iv55))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()

plot_grid(a,b,c,d,e,f,g,nrow=2)

e <- ggplot(df,aes(Hm40,iv55))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_abline(slope=1,intercept=0)+
    theme_classic()
f <- ggplot(df,aes(Hm45,iv55))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_abline(slope=1,intercept=0)+
    theme_classic()
g <- ggplot(df,aes(Hm45,iv38))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_abline(slope=1,intercept=0)+
    theme_classic()
h <- ggplot(df,aes(iv20,iv61))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()


a <- ggplot(df,aes(Hm25,iv20))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()
b <- ggplot(df,aes(Hm30,iv20))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()
c <- ggplot(df,aes(Hm35,iv20))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()
d <- ggplot(df,aes(Hm40,iv20))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()


e <- ggplot(df,aes(Hm25,iv61))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()
f <- ggplot(df,aes(Hm30,iv61))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()
g <- ggplot(df,aes(Hm35,iv61))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()
h <- ggplot(df,aes(Hm40,iv61))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()

iv_20grid <- plot_grid(a,b,c,d,nrow=1)
iv_61grid <- plot_grid(e,f,g,h,nrow=1)
plot_grid(iv_20grid,iv_61grid, ncol=1)


df$MeanMut <- rowMeans(df[, c("iv61", "iv55", "iv38")], na.rm = TRUE)
df$MeanWt <- rowMeans(df[, c("Hm30","iv20","iv33")], na.rm = TRUE)

i <- ggplot(df,aes(MeanWt,iv61))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()
    
j <- ggplot(df,aes(Hm35,iv38))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()
k <- ggplot(df,aes(Hm40,iv38))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()
l <- ggplot(df,aes(Hm45,iv38))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()
m <- ggplot(df,aes(Hm50,iv38))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()
    
    
    
### checking malleti vs rosina
matrix <- GetAssayData(object = subset(Hm_merged,orig.ident %in% c('Hmr_10')))
matrix_mod<-as.matrix(matrix)
cell_count <- ncol(matrix_mod)
min_cells <- ceiling(0.1 * cell_count) # set minimum % of cells gene is expressed in
min_expr <- 0.5  # set minimum expression value
gene_filter <- rowSums(matrix_mod > 0) >= min_cells & rowMeans(matrix_mod) > min_expr
matrix_mod <- matrix_mod[gene_filter, ]
gene<-as.numeric(matrix_mod["Ivory",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
IvoryCorrHmr10 <- sort(correlations)

matrix <- GetAssayData(object = subset(Hm_merged,orig.ident %in% c('Hmm_10')))
matrix_mod<-as.matrix(matrix)
cell_count <- ncol(matrix_mod)
min_cells <- ceiling(0.1 * cell_count) # set minimum % of cells gene is expressed in
min_expr <- 0.5  # set minimum expression value
gene_filter <- rowSums(matrix_mod > 0) >= min_cells & rowMeans(matrix_mod) > min_expr
matrix_mod <- matrix_mod[gene_filter, ]
gene<-as.numeric(matrix_mod["Ivory",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
IvoryCorrHmm10 <- sort(correlations)

matrix <- GetAssayData(object = subset(Hm_merged,orig.ident %in% c('Hmr_25')))
matrix_mod<-as.matrix(matrix)
cell_count <- ncol(matrix_mod)
min_cells <- ceiling(0.1 * cell_count) # set minimum % of cells gene is expressed in
min_expr <- 0.5  # set minimum expression value
gene_filter <- rowSums(matrix_mod > 0) >= min_cells & rowMeans(matrix_mod) > min_expr
matrix_mod <- matrix_mod[gene_filter, ]
gene<-as.numeric(matrix_mod["Ivory",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
IvoryCorrHmr25 <- sort(correlations)

matrix <- GetAssayData(object = subset(Hm_merged,orig.ident %in% c('Hmm_25')))
matrix_mod<-as.matrix(matrix)
cell_count <- ncol(matrix_mod)
min_cells <- ceiling(0.1 * cell_count) # set minimum % of cells gene is expressed in
min_expr <- 0.5  # set minimum expression value
gene_filter <- rowSums(matrix_mod > 0) >= min_cells & rowMeans(matrix_mod) > min_expr
matrix_mod <- matrix_mod[gene_filter, ]
gene<-as.numeric(matrix_mod["Ivory",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
IvoryCorrHmm25 <- sort(correlations)

matrix <- GetAssayData(object = subset(Hm_merged,orig.ident %in% c('Hmr_40')))
matrix_mod<-as.matrix(matrix)
cell_count <- ncol(matrix_mod)
min_cells <- ceiling(0.1 * cell_count) # set minimum % of cells gene is expressed in
min_expr <- 0.5  # set minimum expression value
gene_filter <- rowSums(matrix_mod > 0) >= min_cells & rowMeans(matrix_mod) > min_expr
matrix_mod <- matrix_mod[gene_filter, ]
gene<-as.numeric(matrix_mod["Ivory",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
IvoryCorrHmr40 <- sort(correlations)

matrix <- GetAssayData(object = subset(Hm_merged,orig.ident %in% c('Hmm_40')))
matrix_mod<-as.matrix(matrix)
cell_count <- ncol(matrix_mod)
min_cells <- ceiling(0.1 * cell_count) # set minimum % of cells gene is expressed in
min_expr <- 0.5  # set minimum expression value
gene_filter <- rowSums(matrix_mod > 0) >= min_cells & rowMeans(matrix_mod) > min_expr
matrix_mod <- matrix_mod[gene_filter, ]
gene<-as.numeric(matrix_mod["Ivory",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
IvoryCorrHmm40 <- sort(correlations)


all_genes <- unique(c(
  names(IvoryCorrHmr10),
  names(IvoryCorrHmm10),
  names(IvoryCorrHmr25),
  names(IvoryCorrHmm25),
  names(IvoryCorrHmr40),
  names(IvoryCorrHmm40)
))
# Build data frame
df <- data.frame(
  Gene = all_genes,
  Hmr10 = as.numeric(IvoryCorrHmr10[all_genes]),
  Hmm10 = as.numeric(IvoryCorrHmm10[all_genes]),
  Hmr25 = as.numeric(IvoryCorrHmr25[all_genes]),
  Hmm25 = as.numeric(IvoryCorrHmm25[all_genes]),
  Hmr40 = as.numeric(IvoryCorrHmr40[all_genes]),
  Hmm40 = as.numeric(IvoryCorrHmm40[all_genes])
)

a1 <- ggplot(df,aes(Hmr10,Hmm10))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()

b1 <- ggplot(df,aes(Hmr25,Hmm25))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()
    
c1 <- ggplot(df,aes(Hmr40,Hmm40))+
    geom_point(alpha=0.2)+
    xlim(-1,1)+ylim(-1,1)+
    geom_hline(yintercept=c(-0.25,0,0.25))+geom_vline(xintercept=c(-0.25,0,0.25))+geom_abline(slope=1,intercept=0)+
    theme_classic()

a1+b1+c1

## interesting gene IDs
Hmel.Hmel210001oG184 = WntA
Hmel.Hmel207001oG974 = Fz2
Hmel.Hmel205001oG103 = Fz3
Hmel.Hmel211001oG677 = Notch
Hmel.Hmel209001oG161 = delta
Hmel.Hmel201001oG922 = aristaless 2
Hmel.Hmel201001oG951 = senseless
Hmel.Hmel212001oG165 = sox15
Hmel.Hmel201001oG914 = Wg
Hmel.Hmel208001oG566 = Hth
Hmel.Hmel201001oG412 = Su(H)
Hmel.Hmel203003oG331 = antennapedia
Hmel.Hmel218003oG779 = shaven ; sv
Hmel.Hmel209001oG81 = spineless; ss
Hmel.Hmel218003oG49 = optix
Hmel.Hmel220003oG638 = melt
Hmel.Hmel208001oG166 = Amon
Hmel.Hmel201001oG204 = ocelliless
Hmel.Hmel210001oG901 = dve
Hmel.Hmel211001oG401 = repo; glia marker
Hmel.Hmel211001oG557 = BarH1
Hmel.Hmel217001oG226 = ac (ac-sc, ASH1?)
Hmel.Hmel217001oG232 = sc (ac-sc, ASH1?)
Hmel.Hmel221001oG429 = ap (one homolog)
Hmel.Hmel221001oG423 = ap (one homolog)
Hmel.Hmel207001oG279 = spdo (poor hit)
Hmel.Hmel208001oG149 = pros
Hmel.Hmel202001oG51 = ct
Hmel.Hmel214004oG20 = sox14
Hmel.Hmel218003oG472 = bs ; blistered
Hmel.Hmel210001oG272 = fne ; neuron makrer
Hmel.Hmel206001oG530 = NimC4 ; haemocyte marker
Hmel.Hmel217001oG386 = toe ; canpaniform sensilla marker
Hmel.Hmel202001oG493 = bi / omb (bifid / optomotor blind)
Hmel.Hmel203003oG325 = Ubx
Hmel.Hmel201001oG758 = Fz1
Hmel.Hmel205001oG198 = Fz4
Hmel.Hmel207001oG300 = en
Hmel.Hmel218003oG671 = ci


## checking id of weird iv-only lineage: 

21,26,18,16,20

Idents(Hm_merged) <- 'seurat_clusters'
markersMystery <- FindMarkers(Hm_merged, ident.1=c('21','26','18','16','20'),recorrect_umi=FALSE)
markersMystery$pct.d <- abs(markersMystery$pct.1 - markersMystery$pct.2)



















