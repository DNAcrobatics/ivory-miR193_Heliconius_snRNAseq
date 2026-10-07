library(stringi)
library(dplyr)
library(Seurat)
library(patchwork)
library(foreach)
library(doParallel)
library(SeuratDisk)

#object made on pegasus
# for now we're just using 5i - 30

Hm_merged <- LoadH5Seurat('Hm_combined.full.sct.h5seurat')

DefaultAssay(Hm_merged)  <- "SCT"

## look at some plots 
 a <- DimPlot(Hm_merged, label=T,repel=T, shuffle=T) 
 b <- DimPlot(Hm_merged, label=T,repel=T, group.by='orig.ident', shuffle=T)
 c <- DimPlot(Hm_merged, label=T,repel=T, group.by='time', shuffle=T)
 d <- FeaturePlot(Hm_merged, 'Ivory')
 e <- FeaturePlot(Hm_merged, 'Hmel.Hmel212001oG165') #sox15
 a+b+c+d+e

## retaining only cells with over 1 for sv
Scales <- subset(Hm_merged, Hmel.Hmel218003oG779 > 1)

## this first pass retains some sockets. so also filtering out things with sox15+
Scales <- subset(Scales, Hmel.Hmel212001oG165 < 1)

### getting phase
ccgenes <- read.table('../Ivory_snuc_reannot/CellCycleGenes.txt',h=T)
### make objects for each class
s.genes <- ccgenes$HmID[ccgenes$stage == 'S']
g2m.genes <- ccgenes$HmID[ccgenes$stage == 'G2/M']

### assign each cell a cell cycle score
Scales <- CellCycleScoring(Scales, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# switch idents back
Idents(Scales) <- 'seurat_clusters'

## plot pc1 vs time
plt_dat <- FetchData(Scales, c("PC_1","PC_2", "time","orig.ident","seurat_clusters","nCount_RNA","percent.Rb","Phase","percent.MT"))
ggplot(plt_dat, aes(time, PC_1)) +
  geom_jitter(aes(color = seurat_clusters)) +
  labs(y = "PC_1")+
  theme_classic()

## rerunning PCA on just scales: 
DefaultAssay(Scales) <- 'integrated'
#Scales <- FindVariableFeatures(Scales, selection.method = "vst", nfeatures = 5000)
Scales <- SCTransform(Scales, verbose = T, variable.features.n = 4000,vars.to.regress = c('time','orig.ident'))
Scales <- RunPCA(Scales, features = VariableFeatures(object = Scales))
ElbowPlot(Scales)
Scales <- FindNeighbors(Scales, dims = 1:9)
Scales <- FindClusters(Scales, resolution = 1)
Scales <- RunUMAP(Scales, dims = 1:9)
DimPlot(Scales)

a <- DimPlot(Scales, label=T,repel=T, group.by='time')
b <- DimPlot(Scales, label=T,repel=T)
a+b


SaveH5Seurat(Scales, 'Hm_Scales', overwrite = TRUE, verbose = TRUE)

#with this reclustering, I check the cluster IDs
0 = mito and ribo 
1 = ovo+, apical constrictor nyo,
2 = yellow-e+, osiris, 
3 = Hr38, capD2 (a condensin component)
4 = plum, dve, fzr, (plum and myo, which interact).
5 = salm, Sox102F, sns, 
6 = laza, forked-end, 
7 = laza, moody, dy(dusky, "cell size/shape changes during wing morphogenesis")
8 = yellow-e, ci, 
9  = y (i.e. yellow original flavour), Jhbp (juvenlie hormone binding protein)
10 = (nothing interesting). 
11 = dachsous, NPFR
12 = Antp, BarH1
13 = several cuticular porteins
14 = yellow-g
15 = (nothing interesting). 



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
Hmel.Hmel208001oG658 = Ser ; serrate , notch ligand 

#Hr receptors
Hmel.Hmel220003oG328 = e75, interacts with Hr3
Hmel.Hmel211001oG590 = Hr4
Hmel.Hmel207001oG195 = Hr3
Hmel.Hmel201001oG607 = Hr38
Hmel.Hmel211001oG642 = Hr39
Hmel.Hmel221001oG501 = Hr83
Hmel.Hmel219001oG367 = Hr96
Hmel.Hmel214004oG467 = AkhR

Hrlist <- c('Hmel.Hmel220003oG328',
'Hmel.Hmel211001oG590',
'Hmel.Hmel207001oG195',
'Hmel.Hmel201001oG607',
'Hmel.Hmel211001oG642',
'Hmel.Hmel221001oG501',
'Hmel.Hmel219001oG367',
'Hmel.Hmel214004oG467')




## doing some correlation work 
# first optix 
matrix <- GetAssayData(object = Scales)
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Hmel.Hmel218003oG49",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
OptixCorr_scale <- sort(correlations)

          HMEL000025                Ivory Hmel.Hmel217001oG348 
          -0.3428534           -0.2881212           -0.2708222 
Hmel.Hmel220003oG799 Hmel.Hmel221001oG246 Hmel.Hmel202001oG493 
          -0.2581107           -0.2514992           -0.2299996 
 Hmel.Hmel215003oG87 Hmel.Hmel207001oG217 Hmel.Hmel221001oG247 
          -0.1713775           -0.1596497           -0.1406031 
Hmel.Hmel201001oG923 
          -0.1367556 

Hmel.Hmel217001oG348 
Hmel.Hmel220003oG799 CG42741
Hmel.Hmel221001oG246
 Hmel.Hmel202001oG493 bi
 Hmel.Hmel215003oG87 
 Hmel.Hmel207001oG217 Dll
 Hmel.Hmel221001oG247 dar1
Hmel.Hmel201001oG923 

Hmel.Hmel207001oG300 Hmel.Hmel213001oG153 Hmel.Hmel208001oG160 
           0.1864486            0.1978298            0.1996059 
Hmel.Hmel211001oG619  Hmel.Hmel213001oG51 Hmel.Hmel206001oG108 
           0.2106486            0.2149927            0.2224317 
Hmel.Hmel208001oG159 Hmel.Hmel203003oG329 Hmel.Hmel203003oG331 
           0.2448688            0.2504794            0.2985617 
 Hmel.Hmel218003oG49 
 
Hmel.Hmel207001oG300 en
Hmel.Hmel213001oG153 grn
Hmel.Hmel208001oG160 dsx
Hmel.Hmel211001oG619  obst-A
Hmel.Hmel213001oG51 unc-5
Hmel.Hmel206001oG108 mamo
Hmel.Hmel208001oG159 dsx
Hmel.Hmel203003oG329 (next to Antp)
Hmel.Hmel203003oG331 Antp
 Hmel.Hmel218003oG49 Optix


## and with just post 35%
matrix <- GetAssayData(object = subset(Scales,time > 34))
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Hmel.Hmel218003oG49",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
OptixCorr_scale <- sort(correlations)


## with al1
matrix <- GetAssayData(object = Scales)
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Hmel.Hmel201001oG924",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
Al1Corr_scale <- sort(correlations)

 Hmel.Hmel204001oG67 Hmel.Hmel220003oG241 Hmel.Hmel217001oG594 
          -0.2875280           -0.2875208           -0.2700299 
 Hmel.Hmel207001oG31 Hmel.Hmel210001oG474 Hmel.Hmel203003oG133 
          -0.2697183           -0.2632603           -0.2624699 
Hmel.Hmel221001oG555 Hmel.Hmel204001oG515  Hmel.Hmel207001oG36 
          -0.2621733           -0.2618547           -0.2553814 
Hmel.Hmel204001oG524 
          -0.2547854 
> tail(Al1Corr_scale,10)
Hmel.Hmel211001oG115 Hmel.Hmel207001oG199 Hmel.Hmel217001oG591 
           0.2611248            0.2631337            0.2669164 
Hmel.Hmel217001oG262 Hmel.Hmel201001oG814 Hmel.Hmel218003oG666 
           0.2744328            0.2744360            0.2843982 
Hmel.Hmel212001oG121 Hmel.Hmel201001oG922 Hmel.Hmel201001oG923 
           0.3083745            0.4587508            0.7850940 


Hmel.Hmel201001oG814	mbl
Hmel.Hmel204001oG515	futsch
Hmel.Hmel204001oG67	CG42327
Hmel.Hmel207001oG199	Np
Hmel.Hmel207001oG36	CG13676
Hmel.Hmel210001oG474	reb
Hmel.Hmel212001oG121	CG46301
Hmel.Hmel217001oG262	ovo
Hmel.Hmel218003oG666	CG7512
Hmel.Hmel218003oG666	cas
Hmel.Hmel220003oG241	Ptr
Hmel.Hmel221001oG555	Nep2


## with Hmel.Hmel205001oG103 fz3
matrix <- GetAssayData(object = Scales)
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Hmel.Hmel205001oG103",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
Fz3Corr_scale <- sort(correlations)
> tail(Fz3Corr_scale,10)
 Hmel.Hmel204001oG85 Hmel.Hmel207001oG700 Hmel.Hmel207001oG742 
           0.3141679            0.3168997            0.3220543 
Hmel.Hmel201001oG607 Hmel.Hmel210001oG574 Hmel.Hmel207001oG250 
           0.3302637            0.3315583            0.3353316 
Hmel.Hmel207001oG747 Hmel.Hmel219001oG506 Hmel.Hmel212001oG832 
           0.3366050            0.3382878            0.3547568 
Hmel.Hmel205001oG103 
           1.0000000 
           
Hmel.Hmel201001oG607	Hr38
Hmel.Hmel204001oG85	CG11360
Hmel.Hmel207001oG700	unc-13-4A
Hmel.Hmel212001oG832	Cap-D2

> head(Fz3Corr_scale,10)
Hmel.Hmel207001oG375 Hmel.Hmel202001oG271 Hmel.Hmel213001oG260 
          -0.2766926           -0.2380202           -0.2169920 
 Hmel.Hmel204001oG67 Hmel.Hmel203003oG392 Hmel.Hmel216002oG442 
          -0.2165093           -0.2093164           -0.2092295 
Hmel.Hmel204001oG388 Hmel.Hmel209001oG223 Hmel.Hmel212001oG980 
          -0.2075853           -0.2063142           -0.2055207 
 Hmel.Hmel203003oG18 
          -0.1985547 
Hmel.Hmel202001oG271	CG45065
Hmel.Hmel203003oG18	Cad89D
Hmel.Hmel203003oG392	FER
Hmel.Hmel204001oG67	CG42327
Hmel.Hmel207001oG375	mgl
Hmel.Hmel212001oG980	CG13124
Hmel.Hmel213001oG260	Tgi
Hmel.Hmel216002oG442	DAAM


## correlationw ith ss
matrix <- GetAssayData(object = Scales)
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["Hmel.Hmel209001oG81",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
ssCorr_scale <- sort(correlations)

Hmel.Hmel209001oG288 Hmel.Hmel210001oG101 Hmel.Hmel212001oG111 
          -0.7308514           -0.6816751           -0.6618324 
  Hmel.Hmel209001oG9 Hmel.Hmel201001oG635 Hmel.Hmel212001oG261 
          -0.6530595           -0.6503526           -0.6420627 
Hmel.Hmel210001oG648 Hmel.Hmel211001oG504 Hmel.Hmel221001oG797 
          -0.6363945           -0.6303409           -0.6241603 
Hmel.Hmel220003oG820 
          -0.6168807 
Hmel.Hmel201001oG635	CadN
Hmel.Hmel209001oG288	kmr
Hmel.Hmel209001oG9	CG32694
Hmel.Hmel210001oG101	CG7110
Hmel.Hmel210001oG648	a
Hmel.Hmel211001oG504	CG11353
Hmel.Hmel211001oG504	CG30471
Hmel.Hmel212001oG111	B4
Hmel.Hmel212001oG261	peb
Hmel.Hmel220003oG820	CG31997

> tail(ssCorr_scale,10)
Hmel.Hmel213001oG325 Hmel.Hmel218003oG672 Hmel.Hmel207001oG605 
           0.6497101            0.6522116            0.6543820 
 Hmel.Hmel216002oG72 Hmel.Hmel209001oG214 Hmel.Hmel209001oG229 
           0.6633980            0.6871499            0.6983732 
Hmel.Hmel212001oG922  Hmel.Hmel214004oG86 Hmel.Hmel213001oG451 
           0.7050062            0.7130168            0.7830714 
 Hmel.Hmel209001oG81 
Hmel.Hmel207001oG605	pip
Hmel.Hmel209001oG214	Spn100A
Hmel.Hmel209001oG229	Jupiter
Hmel.Hmel209001oG81	ss
Hmel.Hmel212001oG922	Samuel
Hmel.Hmel213001oG325	Blimp-1
Hmel.Hmel213001oG451	heph
Hmel.Hmel214004oG86	Hil
Hmel.Hmel218003oG672	siz







































