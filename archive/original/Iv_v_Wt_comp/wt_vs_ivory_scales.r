library(Seurat)
library(ggplot2)
library(sctransform)

# read data
iv_20.data <- Read10X(data.dir="iv_20.filtered")
iv_33.data <- Read10X(data.dir="iv_33.filtered")
iv_38.data <- Read10X(data.dir="iv_38.filtered")
iv_55.data <- Read10X(data.dir="iv_55.filtered")
iv_61.data <- Read10X(data.dir="iv_61.filtered")
Hmm_40.data <- Read10X(data.dir='filtered.Hmm_40')
Hmr_35.data <- Read10X(data.dir='filtered.Hmr_35')
Hmr_40.data <- Read10X(data.dir='filtered.Hmr_40')
Hmr_45.data <- Read10X(data.dir='filtered.Hmr_45')

#use sctransform
iv_20 <- CreateSeuratObject(counts = iv_20.data, project = "iv_20", min.cells = 30, min.features = 1000)
iv_33 <- CreateSeuratObject(counts = iv_33.data, project = "iv_33", min.cells = 30, min.features = 1000)
iv_38 <- CreateSeuratObject(counts = iv_38.data, project = "iv_38", min.cells = 30, min.features = 1000)
iv_55 <- CreateSeuratObject(counts = iv_55.data, project = "iv_55", min.cells = 30, min.features = 1000)
iv_61 <- CreateSeuratObject(counts = iv_61.data, project = "iv_61", min.cells = 30, min.features = 1000)
Hmm_40 <- CreateSeuratObject(counts= Hmm_40.data, project= 'Hmm_40',min.cells= 30, min.features= 1000)
Hmr_35 <- CreateSeuratObject(counts= Hmr_35.data, project= 'Hmr_35',min.cells= 30, min.features= 1000)
Hmr_40 <- CreateSeuratObject(counts= Hmr_40.data, project= 'Hmr_40',min.cells= 30, min.features= 1000)
Hmr_45 <- CreateSeuratObject(counts= Hmr_45.data, project= 'Hmr_45',min.cells= 30, min.features= 1000)

# use sctransform, regressing on count
iv_20  <- SCTransform(iv_20 , method = "glmGamPoi", vars.to.regress = "nCount_RNA", verbose = FALSE)
iv_33  <- SCTransform(iv_33 , method = "glmGamPoi", vars.to.regress = "nCount_RNA", verbose = FALSE)
iv_38  <- SCTransform(iv_38 , method = "glmGamPoi", vars.to.regress = "nCount_RNA", verbose = FALSE)
iv_55  <- SCTransform(iv_55 , method = "glmGamPoi", vars.to.regress = "nCount_RNA", verbose = FALSE)
iv_61  <- SCTransform(iv_61 , method = "glmGamPoi", vars.to.regress = "nCount_RNA", verbose = FALSE)
Hmm_40 <- SCTransform(Hmm_40, method = "glmGamPoi", vars.to.regress = "nCount_RNA", verbose = FALSE)
Hmr_35 <- SCTransform(Hmr_35, method = "glmGamPoi", vars.to.regress = "nCount_RNA", verbose = FALSE)
Hmr_40 <- SCTransform(Hmr_40, method = "glmGamPoi", vars.to.regress = "nCount_RNA", verbose = FALSE)
Hmr_45 <- SCTransform(Hmr_45, method = "glmGamPoi", vars.to.regress = "nCount_RNA", verbose = FALSE)

# make object list
object.list = c(iv_20,
                iv_33,
                iv_38,
                iv_55,
                iv_61,
                Hmm_40,
                Hmr_35,
                Hmr_40,
                Hmr_45)

# run PrepSCTIntegration (does some necessary data reformatting)
features <- SelectIntegrationFeatures(object.list = wing.list, nfeatures = 3000)
wing.list <- PrepSCTIntegration(
  object.list,
  assay = NULL,
  anchor.features = features,
  sct.clip.range = NULL,
  verbose = TRUE
)

# find anchors and integrate 
wing.anchors <- FindIntegrationAnchors(object.list = wing.list, normalization.method = "SCT",anchor.features = features)
#(had to do this and the PCA step on pegasus)
wing.sct <- IntegrateData(anchorset = wing.anchors, normalization.method = "SCT")
# run PCA
wing.sct <- RunPCA(wing.sct, verbose = FALSE)

# select dims
ElbowPlot(wing.sct)
# selected 10

wing.sct <- RunUMAP(wing.sct, reduction = "pca", dims = 1:10)
wing.sct <- FindNeighbors(wing.sct, reduction = "pca", dims = 1:10)
wing.sct <- FindClusters(wing.sct, resolution = 0.5)


p1 <- DimPlot(wing.sct, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(wing.sct, reduction = "umap", group.by = "seurat_clusters", label = TRUE,repel = TRUE)
p1 + p2


### cell cycle 

### import list of canonical cell cycle markers
ccgenes <- read.table('../Ivory_snuc_reannot/CellCycleGenes.txt',h=T)
### make objects for each class
s.genes <- ccgenes$HmID[ccgenes$stage == 'S']
g2m.genes <- ccgenes$HmID[ccgenes$stage == 'G2/M']

### make copy of fre40 for tests
test <- wing.sct

### assign each cell a cell cycle score
test <- CellCycleScoring(wing.sct, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# run a pca on cell cycle genes to see if cells split out this way
test <- RunPCA(test, features = c(s.genes, g2m.genes))
PCAPlot(test)
DimPlot(test)


plt_dat <- FetchData(test, c("PC_1", "Phase","orig.ident"))
ggplot(plt_dat, aes(orig.ident, PC_1)) +
  geom_jitter(aes(color = Phase)) +
  labs(y = "PC1 (aka pseudotime)")

####trying using an alt protocl where you look at the difference rather than the phase. 
test$CC.Difference <- test$S.Score - test$G2M.Score

test <- SCTransform(test, method = "glmGamPoi", vars.to.regress = c("Phase"), verbose = FALSE)

test <- RunPCA(test)
PCAPlot(test)

# repeat the pc1 vs time plot:
plt_dat <- FetchData(test, c("PC_1", "PC_2","Phase","orig.ident","seurat_clusters"))
ggplot(plt_dat, aes(orig.ident, PC_1,col=Phase)) +
  geom_jitter(aes(color = seurat_clusters)) +
  labs(y = "PC1 (aka pseudotime)")

##rerun umap and clsutering
test <- RunUMAP(test, reduction = "pca", dims = 1:15)
test <- FindNeighbors(test, reduction = "pca", dims = 1:15)
test <- FindClusters(test, resolution = 0.5)

DimPlot(test, reduction = "umap", group.by = "orig.ident",shuffle=TRUE)
DimPlot(test, reduction = "umap", label = TRUE, repel = TRUE)

# cell cycle effects strongly mitigated in PCA
test <- RunPCA(test, features = VariableFeatures(test), nfeatures.print = 10)

### using the new (un-cell-cycle-sorted) object, look for differences between iv_61 and Hmr_35


Idents(wing.sct) <- "orig.ident"
avg.wing.sct <- as.data.frame(log1p(AverageExpression(wing.sct, verbose = FALSE)$SCT))
avg.wing.sct$gene <- rownames(avg.wing.sct)
ggplot(avg.wing.sct, aes(Hmr_35,iv_61)) + geom_point() 
ggplot(avg.wing.sct, aes(Hmr_35,iv_20)) + geom_point()


# then looking at DE
Hmr_35viv_61 <- FindMarkers(wing.sct, ident.1 = "Hmr_35", ident.2 = "iv_61", verbose = FALSE)
# volcoano plot
ggplot(Hmr_35viv_61,aes(avg_log2FC, -log10(p_val)))+geom_point()

# then looking at DE
iv20_35viv_61 <- FindMarkers(wing.sct, ident.1 = "iv_20", ident.2 = "iv_61", verbose = FALSE)
# volcoano plot
ggplot(iv20_35viv_61,aes(avg_log2FC, -log10(p_val)))+geom_point()


#### checking ID of cluster 13 from wing.sct
Idents(wing.sct) <- "seurat_clusters"
wingSCT.cluster13 <- FindMarkers(wing.sct, ident.1 = "13", verbose = FALSE)

Idents(wing.sct) <- "Phase"
wingSCT.sphase <- FindMarkers(wing.sct, ident.1 = "S", verbose = FALSE)
> wingSCT.sphase <- FindMarkers(wing.sct, ident.1 = "S", verbose = FALSE)
> head(wingSCT.sphase)
                              p_val avg_log2FC pct.1 pct.2     p_val_adj
Hmel.Hmel208001oG270   0.000000e+00  5.9298499 0.460 0.195  0.000000e+00	# claspin ATR-Chk1 checkpoint pathway.
Hmel.Hmel210001oG1026  0.000000e+00  4.9942209 0.441 0.235  0.000000e+00	# cyce control of the cell cycle at the G1/S (start) transition
Hmel.Hmel211001oG151   0.000000e+00  3.8460396 0.433 0.195  0.000000e+00	# pola1 DNA polymerase alpha subunit 1
Hmel.Hmel219001oG202   0.000000e+00  3.0230275 0.410 0.186  0.000000e+00	# mcm6 component of eukaryotic replicative DNA helicase
Hmel.Hmel218003oG276  4.834153e-273 -0.6557956 0.409 0.199 1.450246e-269	# spel1 post replication mismatch repair
Hmel.Hmel217001oG186  5.592355e-270  7.0895379 0.402 0.194 1.677707e-266	#



> ggplot(iv20_35viv_61,aes(avg_log2FC, -log10(p_val)))+geom_point()
> head(iv20_35viv_61)
                     p_val avg_log2FC pct.1 pct.2 p_val_adj
Hmel.Hmel221001oG580     0  -1.102666 0.185 0.108         0	# aracaun
Hmel.Hmel207001oG406     0  -6.483671 0.348 0.235         0	# NetA netrin A 
Hmel.Hmel210001oG318     0  -4.672710 0.191 0.152         0	# Nox NADPH oxidase 
Hmel.Hmel221001oG586     0  -4.673380 0.310 0.074         0	# mirr
Hmel.Hmel213001oG853     0  -6.836317 0.313 0.019         0	# cg15020
Hmel.Hmel219001oG506     0  -5.929883 0.289 0.023         0	# gl (glass, cell fate spec in retina)

### trying to pull out just the scale cells
Idents(wing.sct) <- "seurat_clusters"
cluster_id <- c(15,8,2)
cells_of_interest <- which(Idents(wing.sct) %in% cluster_id)

# Subset the Seurat object to include only the cells in the cluster of interest
scales <- subset(wing.sct, cells = cells_of_interest)


#t.cells <- subset(immune.combined, idents = "CD4 Naive T")
Idents(scales) <- "orig.ident"
avg.scales <- as.data.frame(log1p(AverageExpression(scales, verbose = FALSE)$RNA))
avg.scales$gene <- rownames(avg.scales)


ggplot(avg.scales, aes(Hmr_30,iv_61)) + geom_point() 
wt_vs_mut <- FindMarkers(scales, ident.1 = "iv_20", ident.2 = "iv_61", verbose = FALSE)

Idents(scales) <- "orig.ident"

wt_vs_mut <- FindMarkers(scales, ident.1 = c("Hmr_35","Hmr_40","Hmr_45","Hmm_40","iv_20"), ident.2 = c("iv_55","iv_61"), verbose = FALSE)


> head(wt_vs_mut)
                             p_val avg_log2FC pct.1 pct.2     p_val_adj
Hmel.Hmel215003oG537 1.234173e-127 -2.0108807 0.364 0.165 3.702520e-124	# CG15093
Hmel.Hmel220003oG443 1.352425e-124 -4.3023580 0.589 0.196 4.057275e-121	# cg10738
Hmel.Hmel201001oG922  1.295450e-91 -2.5015923 0.248 0.149  3.886349e-88	# aristaless (not sure if 1 or 2)
Hmel.Hmel217001oG338  5.617977e-81 -3.4523949 0.466 0.205  1.685393e-77	# nohit
Hmel.Hmel201001oG825  7.018918e-81 -1.9298380 0.392 0.176  2.105675e-77	# zw glucose-6-phosphate dehydrogenase (NADP
Hmel.Hmel205001oG360  2.163427e-75 -0.7709498 0.188 0.185  6.490280e-72	# rpn11 Metalloprotease component of the 26S proteasome
Hmel.Hmel210001oG295  1.739985e-74 -4.6158845 0.509 0.262  5.219955e-71	# cg9447 main hit
Hmel.Hmel219001oG245  1.847978e-74 -0.3025080 0.313 0.162  5.543934e-71	# 
Hmel.Hmel213001oG870  6.095199e-74  0.7041369 0.114 0.099  1.828560e-70	# 
Hmel.Hmel208001oG377  8.618055e-74  1.7571639 0.306 0.148  2.585416e-70	# 
Hmel.Hmel208001oG292  4.138067e-73  0.7721912 0.169 0.147  1.241420e-69	# 
Hmel.Hmel214004oG598  7.190970e-71 -3.1277568 0.304 0.176  2.157291e-67	#  #cg9447 homolog
Hmel.Hmel206001oG357  5.644472e-70 -1.5023649 0.338 0.198  1.693341e-66	# 
Hmel.Hmel220003oG901  1.033341e-68  0.4293162 0.354 0.229  3.100022e-65	# 
Hmel.Hmel206001oG302  1.170702e-68 -0.5742232 0.391 0.252  3.512105e-65	# 
Hmel.Hmel217001oG276  1.771000e-67 -2.2204097 0.145 0.078  5.313001e-64	# 
Hmel.Hmel212001oG87   2.907347e-67 -8.7556742 0.522 0.258  8.722041e-64	# SIFaR
Hmel.Hmel208001oG231  4.265088e-66  2.1065273 0.338 0.160  1.279527e-62	# 
Hmel.Hmel211001oG794  2.430956e-65 -1.9710843 0.411 0.152  7.292869e-62	# 
Hmel.Hmel204001oG406  9.569063e-65 -3.8955032 0.474 0.226  2.870719e-61	# fructose-1,6-bisphosphatase



## running same de on everything
Idents(wing.sct) <- "orig.ident"
wt_vs_mutWh <- FindMarkers(wing.sct, ident.1 = c("Hmr_35","Hmr_40","Hmr_45","Hmm_40","iv_20"), ident.2 = c("iv_55","iv_61"), verbose = FALSE)



### correlations RUN ON SCT, not the default RNA slot!!!
matrix <- GetAssayData(object = wing.sct[["SCT"]], slot = "counts")
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Ivory",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
IvoryCorr <- sort(correlations)

#### first for all cells, figure out which genes are Differentially correlated
#### between wts and iv-s 
## subset for just the correlations in the NonIvories:
wing.sctWt <- subset(wing.sct, idents=c("Hmr_35","Hmr_40","Hmm_40","iv_20") )

matrix <- GetAssayData(object = wing.sctWt[["SCT"]], slot = "scale.data")
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Ivory",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
IvoryWtCorr <- sort(correlations)


## subset for just the correlations in the Ivory muts:
wing.sctMut <- subset(wing.sct, idents=c("iv_55","iv_61") )

matrix <- GetAssayData(object = wing.sctMut[["SCT"]], slot = "scale.data")
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Ivory",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
IvoryMutCorr <- sort(correlations)

# make a data frame of the correlations in wt vs mut
a <- as.data.frame(IvoryWtCorr)
b <- as.data.frame(IvoryMutCorr)
a$gene <- row.names(a)
b$gene <- row.names(b)
c <- merge(a,b)
c$diff <- c$IvoryWtCorr - c$IvoryMutCorr

# make a correlation plot with an abline showing the x=y line
plot(c$IvoryWtCorr,c$IvoryMutCorr)
abline(0,1,col='red')
abline(0.2,1,col='blue')
abline(-0.2,1,col='blue')
abline(0.3,1,col='green')
abline(-0.3,1,col='green')
abline(0.4,1,col='green')
abline(-0.4,1,col='green')

c[abs(c$IvoryWtCorr - c$IvoryMutCorr) > 0.3,]
                     gene IvoryWtCorr IvoryMutCorr       diff
100  Hmel.Hmel201001oG566 -0.05148062   0.25949805 -0.3109787	# cg2225
1588 Hmel.Hmel212001oG540  0.51564189   0.10699589  0.4086460	# CG16711
1959 Hmel.Hmel215003oG103  0.04686866   0.50978376 -0.4629151	# RNAseZ
1960 Hmel.Hmel215003oG104  0.01500571   0.42522703 -0.4102213	# no good hit (but gene overlaps with G93 in RC). 
2067  Hmel.Hmel215003oG90  0.12921627   0.54154429 -0.4123280	# parn
2068  Hmel.Hmel215003oG91 -0.01966591   0.35484293 -0.3745088	# Manf (evolutionarily conserved neurotrophic factor. It contributes to dopamine neuron survival, regulation of unfolded protein response and retinal repair. [)
2069  Hmel.Hmel215003oG92 -0.01818038   0.51294297 -0.5311234	# Rpn10 (encodes one of the polyubiquitin receptor subunits of the 26S proteasome)
2070  Hmel.Hmel215003oG93  0.03666244   0.74974336 -0.7130809	# either Ror ("related to vertebrate neural-specific Wnt-binding Ror receptors") OR Cadherin 96Ca (RTK, target of the tf Grh). 
2091 Hmel.Hmel216002oG164  0.52971844   0.22740104  0.3023174	# cg4313
2484 Hmel.Hmel218003oG749  0.44053739   0.10034412  0.3401933	# rn rotund (rn) encodes a member of the Kruppel zinc-finger transcription factor family
2717 Hmel.Hmel220003oG241  0.46406504   0.13083238  0.3332327	# Ptr - patched-related


## this was all for every cell, i wil quickly repeat for just scales to see if big diff
scales.Wt <- subset(scales, idents=c("Hmr_35","Hmr_40","Hmm_40","iv_20") )

matrix <- GetAssayData(object = scales.Wt[["SCT"]], slot = "scale.data")
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Ivory",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
IvoryWtCorrSc <- sort(correlations)


## subset for just the correlations in the Ivory muts:
scales.Mut <- subset(scales, idents=c("iv_55","iv_61") )

matrix <- GetAssayData(object = scales.Mut[["SCT"]], slot = "scale.data")
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Ivory",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
IvoryMutCorrSc <- sort(correlations)

a <- as.data.frame(IvoryWtCorrSc)
b <- as.data.frame(IvoryMutCorrSc)
a$gene <- row.names(a)
b$gene <- row.names(b)
d <- merge(a,b)
d$diff <- d$IvoryWtCorr - d$IvoryMutCorr
# make plot
plot(d$IvoryWtCorrSc,d$IvoryMutCorrSc)
abline(0,1,col='red')
abline(0.2,1,col='blue')
abline(-0.2,1,col='blue')
abline(0.3,1,col='green')
abline(-0.3,1,col='green')
abline(0.4,1,col='green')
abline(-0.4,1,col='green')

                     gene IvoryWtCorrSc IvoryMutCorrSc       diff
100  Hmel.Hmel201001oG566   -0.17683610    0.132288299 -0.3091244	# cg2225
188  Hmel.Hmel201001oG987   -0.16605474    0.137518321 -0.3035731	# cg31678
255  Hmel.Hmel202001oG493    0.09997249   -0.261686768  0.3616593	# bi (bifid)
334  Hmel.Hmel203003oG281   -0.28946463    0.145741967 -0.4352066	# neo Involved in apical constriction and regulation of embryonic cell shape
805  Hmel.Hmel207001oG516   -0.09134343    0.306215258 -0.3975587	# Supervillin Predicted to enable actin filament binding activity and phosphatidylinositol-4,5-bisphosphate binding activit
848  Hmel.Hmel207001oG695   -0.09380705    0.294443471 -0.3882505	# Msr-110 Predicted to be located in membrane
1022 Hmel.Hmel209001oG214   -0.08117838    0.230283948 -0.3114623	# Serpin 100A
1054 Hmel.Hmel209001oG383    0.18211024   -0.188629463  0.3707397	# polychaetoid - cytosckeletal dynamics
1091  Hmel.Hmel209001oG81    0.37219455   -0.003130168  0.3753247	# spineless
1164 Hmel.Hmel210001oG307   -0.27426005    0.039433745 -0.3136938	# serrano encodes an apically-enriched cytosolic protein that contributes to trachea development and to the planar cell polarity pathway through binding to the product of dsh. 
1262 Hmel.Hmel210001oG775   -0.04514158    0.262993905 -0.3081355	# pnut the septin family of polymerizing GTPases
1588 Hmel.Hmel212001oG540    0.32758772   -0.122498867  0.4500866	# CG16711
1744 Hmel.Hmel213001oG409    0.19118868   -0.250595590  0.4417843	# dnc dunce (dnc) encodes a cAMP-specific phosphodiesterase, responsible for cAMP degradation
1745 Hmel.Hmel213001oG410    0.10520528   -0.249983467  0.3551887	# NOTE PROXIMITY TO DUNCE - (no hit to this gene)
1779 Hmel.Hmel213001oG630    0.28154647   -0.119756040  0.4013025	# cad96ca
1959 Hmel.Hmel215003oG103    0.04258280    0.441999850 -0.3994171	# RNAseZ
1960 Hmel.Hmel215003oG104   -0.01769747    0.325879910 -0.3435774	# no good hit (but gene overlaps with G93 in RC). 
2067  Hmel.Hmel215003oG90    0.06403477    0.508539626 -0.4445049	# parn
2068  Hmel.Hmel215003oG91   -0.01990207    0.410395494 -0.4302976	# Manf (evolutionarily conserved neurotrophic factor. It contributes to dopamine neuron survival, regulation of unfolded protein response and retinal repair. [)
2069  Hmel.Hmel215003oG92    0.03626049    0.448678649 -0.4124182	# Rpn10 (encodes one of the polyubiquitin receptor subunits of the 26S proteasome)
2070  Hmel.Hmel215003oG93   -0.07041945    0.684282298 -0.7547017	# either Ror ("related to vertebrate neural-specific Wnt-binding Ror receptors") OR Cadherin 96Ca (RTK, target of the tf Grh). 
2168 Hmel.Hmel216002oG578    0.28007680   -0.088121626  0.3681984	# Tubulin tyrosine ligase-like 15
2459  Hmel.Hmel218003oG65    0.03451778   -0.274965018  0.3094828	# no good hit
2717 Hmel.Hmel220003oG241    0.37745187    0.046956186  0.3304957	# Ptr - patched-related
2797 Hmel.Hmel220003oG580    0.06787870   -0.279847803  0.3477265	# RhoGAP19D (involved in imaginal leg disc )
2982 Hmel.Hmel221001oG788    0.17594263   -0.229465019  0.4054076	# fax - failed axon connections


# make a correlation plot with an abline showing the x=y line
plot(d$IvoryWtCorrSc,d$IvoryMutCorrSc)
abline(0,1,col='red')

c[abs(c$IvoryWtCorr - c$IvoryMutCorr) > 0.1,]

# making one big data frame of correlations and expression levels
e <- merge(c,d,by='gene')
f <- merge(e,avg.wing.sct)

#plot the 'differential correlation' for whole wt and whole mut. 
#i also coloured the points by expression level in iv_61 
#and adding an abline. 
ggplot(f,aes(IvoryMutCorr,IvoryWtCorr,col=log(iv_20)))+
   geom_point(alpha=0.4)+
   geom_abline(slope=1, intercept=0, col='red')+
   geom_abline(slope=1, intercept=0.2, col='blue')+
   geom_abline(slope=1, intercept=-0.2, col='blue')+
   theme_classic()

# make a column of the magnitude of difference between wt and mut
f$expDiff <- abs((f$iv_20+f$Hmm_40+f$Hmr_35+f$Hmr_40+f$Hmr_45)/5 - (f$iv_55+f$iv_61)/2)

intercept <- c(-0.2,0.2,-0.4,0.4,-0.6,0.6)
plot1 <- ggplot(f,aes(IvoryMutCorr,IvoryWtCorr))+
   geom_point(alpha=0.7)+
   geom_abline(slope=1, intercept=0, col='red')+
   geom_abline(slope=1, intercept=intercept, col='blue')+
   geom_smooth(method=glm)+
   theme_classic()

plot2 <- ggplot(f,aes(IvoryMutCorrSc,IvoryWtCorrSc,))+
   geom_point(alpha=0.7)+
   geom_abline(slope=1, intercept=0, col='red')+
   geom_abline(slope=1, intercept=intercept, col='blue')+
   geom_smooth(method=glm)+
   theme_classic()
   
plot1+plot2

ggplot(f,aes(iv_20,iv_61,col=abs(IvoryWtCorrSc-IvoryMutCorrSc)))+
   geom_point(alpha=0.7)+
   geom_abline(slope=1, intercept=0, col='red')+
   theme_classic()

### plotting whole vs sc correlations for wts 
ggplot(f,aes(IvoryWtCorr,IvoryWtCorrSc,))+
   geom_point(alpha=0.7)+
   geom_abline(slope=1, intercept=0, col='red')+
   geom_abline(slope=1, intercept=intercept, col='blue')+
   geom_smooth(method=glm)+
   theme_classic()


(iv_20+Hmm_40+Hmr_35+Hmr_40+Hmr_45+iv_55+iv_61)/7)

#### making a plot of the topp diff-corr genes vs exp level. 

genelist <- c[abs(c$IvoryWtCorr - c$IvoryMutCorr) > 0.2,1]

# make data.frame of average gene expression from feature RNA
avg.wing.sct <- as.data.frame(log1p(AverageExpression(wing.sct, verbose = FALSE)$RNA))
avg.wing.sct$gene <- rownames(avg.wing.sct)

avg.wing.sct <- melt(avg.wing.sct)

genelist <- c('Hmel.Hmel215003oG103',
              'Hmel.Hmel215003oG104',
              'Hmel.Hmel215003oG87',
              'Hmel.Hmel215003oG90',
              'Hmel.Hmel215003oG91',
              'Hmel.Hmel215003oG92',
              'Hmel.Hmel215003oG93')
### or alt genelist
genelist <- c[abs(c$diff) > 0.3,1]

wt <- c("Hmr_35","Hmr_40","Hmm_40","iv_20")
mut <- c("iv_55","iv_61","iv_33")


avg.wing.sct$label <- ifelse(avg.wing.sct$variable %in% wt, "wt", 
                           ifelse(avg.wing.sct$variable %in% mut, "mut", "check"))


ggplot(avg.wing.sct[avg.wing.sct$gene %in% genelist & avg.wing.sct$label != 'check',], aes(label,value,group=label, col=label))+
    geom_violin()+
    geom_point()+
    facet_wrap(~gene,scales="free_y")



genelistSc <- d[abs(d$diff) > 0.3,1]

# make data.frame of average gene expression from feature RNA
avg.scales <- as.data.frame(log1p(AverageExpression(scales, verbose = FALSE)$RNA))
avg.scales$gene <- rownames(avg.scales)

avg.scales <- melt(avg.scales)

avg.scales$label <- ifelse(avg.scales$variable %in% wt, "wt", 
                           ifelse(avg.scales$variable %in% mut, "mut", "check"))

ggplot(avg.scales[avg.scales$gene %in% genelistSc & avg.scales$label != 'check',], aes(label,value,col=label))+
    geom_violin()+
    geom_point()+
    facet_wrap(~gene,scales="free_y")






