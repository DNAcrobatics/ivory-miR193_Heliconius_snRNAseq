################ TURNS OUT millefy won't work so i've got to jerry-rig this 

library(Rsamtools)		## for loading bam files (not sure how this is gonna roll)
library(GenomicRanges)	## for (as I understand it) dealing with gff / gtf annotations in R
library(dplyr)
library(ggplot2)
library(cowplot)
library(khroma)
library(ggrastr)

#### first i need to make an "id file" for each individual, with a list of barcodes from that indiv
# (helped here by ChatGPT)

# Get barcodes (column names of Seurat object)
full_barcodes <- colnames(test)
# Strip the prefix before the last underscore
# This handles cases like "Hmm_10_AAACCCAGTGCGGATA" => "AAACCCAGTGCGGATA"
stripped_barcodes <- sub(".*_", "", full_barcodes)
# Extract cluster assignments
clusters <- data.frame(
  barcode = colnames(test),
  celltype = test$Iv,
  indiv = test$orig.ident,
  time = test$time)

data_frame <- data.frame(
  barcode = stripped_barcodes, clusters)

write.table(data_frame, "bcdf_test_ivhivslo.all.txt",quote=FALSE,row.names=FALSE)

## Identify cells that belong to the specified individual
#cells_of_interest <- wing.sct@meta.data$orig.ident == 'iv_20'
#
## Get the barcodes of cells belonging to the specified individual
#barcodes_of_interest <- rownames(wing.sct@meta.data)[cells_of_interest]

##turn the above into a for loop
#indivs <- c('iv_20',
#            'iv_33',
#            'iv_38',
#            'iv_55',
#            'iv_61',
#            'Hmm_40',
#            'Hmr_35',
#            'Hmr_40',
#            'Hmr_45')
#
#
#for(i in indivs){
#   coi <- wing.sct@meta.data$orig.ident == i
#   boi <- rownames(wing.sct@meta.data)[coi]
#   clus <- wing.sct@meta.data$seurat_clusters[coi]
#   bcdf <- data.frame(CellBarcode = boi, Cluster = clus)
#   write.table(bcdf, paste0(i, "_barcodes.tsv"),row.names=FALSE)
#}

# unlist function, needed below
.unlist <- function (x)
 {
 ## do.call(c, ...) coerces factor to integer, which is undesired
 x1 <- x[[1L]]
 if (is.factor(x1)) {
 structure(unlist(x), class = "factor", levels = levels(x1))
 } else {
 do.call(c, x)
 }
 }

#############forward plots
#### Code for making a figure of bam alignments sorted by cell and faceted by cluster

what <- c("rname", "pos", "qwidth")
param <- ScanBamParam(what=what,tag="CB")
bam1 <- scanBam('bams/iv_20.CortexLocus.fw.bam', param=param)


bam1 <- unname(bam1) # names not useful in unlisted result
#elts <- setNames(bamWhat(param), bamWhat(param))
#lst <- lapply(elts, function(elt) .unlist(lapply(bam1, "[[", elt)))

test <- as.data.frame(do.call(cbind, bam1))
hist_frame <- bind_rows(test, .id = "indiv")
test <- as.data.frame(hist_frame$V1)
test <- test[order(test$CB),] 

### a bit of data management
# some rows have - in the cell barcode, these should be dumped
test <- test[test$CB!='-',]

## keeping all cells that had at least 10 reads (could be more stringent)
test2 <- test %>% group_by(CB) %>% filter(n() > 5)
orig.idents <- read.table('bcdf_test_ivhivslo.all.txt',h=T)
colnames(orig.idents)[colnames(orig.idents) == "barcode"] <- "CB"
merged <- merge(test2, orig.idents)
merged$CB <- as.factor(merged$CB)

## ggplot code for annotation 
## read in your gff
ANNOT <- read.table("Hmel_cortexIvory_region_luca-reannotated.gff",sep="\t")
names(ANNOT) <- c("contig", "HGC", "type", "con_start", "con_end", "dot", "str", "unk", "descr")

## this for-loop makes a list of components that will be used to draw gene lines
# + orientation genes will go above - orientation genes
ANNOTgene <- data.frame(xmin=NA,xmax=NA,ymin=NA,ymax=NA)
# gene loop 
for (g in 1:nrow(ANNOT)){
  if (ANNOT$type[g] == "gene" && ANNOT$str[g] == "-")
    ANNOTgene[g,] <- c(ANNOT$con_start[g], bot+(abs(top+bot))*0.25, ANNOT$con_end[g], bot+(abs(top+bot))*0.25)
  if (ANNOT$type[g] == "gene" && ANNOT$str[g] == "+")
    ANNOTgene[g,] <- c(ANNOT$con_start[g], bot+(abs(top+bot))*0.75, ANNOT$con_end[g], bot+(abs(top+bot))*0.75)
  ANNOTgene <- na.omit(ANNOTgene)
}

## this for-loop makes a list of components that will be used to draw exon rectangles
# again, + orientation genes will go above - orientation genes
ANNOTexon <- data.frame(xmin=NA,xmax=NA,ymin=NA,ymax=NA)
for (g in 1:nrow(ANNOT)){
  if (ANNOT$type[g] == "exon" && ANNOT$str[g] == "-")
    ANNOTexon[g,] <- c(ANNOT$con_start[g], bot+(abs(top+bot))*0.1, ANNOT$con_end[g], bot+(abs(top+bot))*0.4)
  if (ANNOT$type[g] == "exon" && ANNOT$str[g] == "+")
    ANNOTexon[g,] <- c(ANNOT$con_start[g], bot+(abs(top+bot))*0.6, ANNOT$con_end[g], bot+(abs(top+bot))*0.9)
  ANNOTexon <- na.omit(ANNOTexon)
}

## the key points in this ggplot command are the 'annotate' and 'geom_segment' lines, 
## which can either be added into an existing plot or into an empty plot. 
ann <- ggplot(merged,aes(pos,CB))+
    xlim(1300000,1700000)+
    annotate("rect", xmin=ANNOTexon[,1], xmax=ANNOTexon[,3], ymin=ANNOTexon[,2], ymax=ANNOTexon[,4])+
    geom_segment(data=ANNOTgene, aes(x=ANNOTgene[,1],xend=ANNOTgene[,3],y=ANNOTgene[,2],yend=ANNOTgene[,2]))+
    theme_classic()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

merged$CB <- factor(merged$CB, levels = names(sort(table(merged$CB), decreasing = FALSE)))
fw20 <- ggplot(merged,aes(pos,CB,))+
    geom_point(size=0.3,col='blue',shape=95)+
    facet_grid(celltype~.,scales='free',space='free_y')+
    theme_classic()+
    xlim(1300000,1700000)+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(strip.background = element_blank(),strip.text = element_blank())

plot_grid(fw20, ann, ncol=1,rel_heights=c(1,0.1))

### plot also the 'density' summary ribbon
df2 <- as.data.frame(table(merged$pos))
df2$Var1 <- as.numeric(as.character(df2$Var1))
df2$Var1bin <- cut(df2$Var1, 1000)

# Number of rows per bin
bin_size <- 100
# Create grouping index
groups <- ceiling(seq_along(df2$Freq) / bin_size)
# Aggregate Freq and get the 'left' Var1 value (first in each bin)
binned_df <- do.call(rbind, lapply(split(df2, groups), function(chunk) {
  data.frame(
    Var1 = as.numeric(chunk$Var1[1]),  # left value
    Freq = sum(chunk$Freq)
  )
}))

binned_df <- binned_df %>% arrange(Freq)
fw20dens <- ggplot(binned_df, aes(x = Var1, y = 1,fill=log10(Freq))) + 
  geom_rect(xmin=binned_df$Var1,xmax=binned_df$Var1+500, ymin=0,ymax=2,aes(),width=10,colour=NA)+ 
  scale_fill_batlow(reverse=T,limits=c(2,3))+
  xlim(1300000,1700000)+ theme_void()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(strip.background = element_blank(),strip.text = element_blank())+theme(legend.position = "none")


limits=c(0.1,2.6)
### plot all together 
plot_grid(fw20, fw20dens, ann, ncol=1,rel_heights=c(1,0.2,0.1))


###### REPEATING FOR IV_61
param <- ScanBamParam(what=what,tag="CB")
bam1 <- scanBam('bams/iv_61.CortexLocus.fw.bam', param=param)

bam1 <- unname(bam1) # names not useful in unlisted result
#elts <- setNames(bamWhat(param), bamWhat(param))
#lst <- lapply(elts, function(elt) .unlist(lapply(bam1, "[[", elt)))

test <- as.data.frame(do.call(cbind, bam1))
hist_frame <- bind_rows(test, .id = "indiv")
test <- as.data.frame(hist_frame$V1)
test <- test[order(test$CB),] 

### a bit of data management
# some rows have - in the cell barcode, these should be dumped
test <- test[test$CB!='-',]

## keeping all cells that had at least 10 reads (could be more stringent)
test2 <- test %>% group_by(CB) %>% filter(n() > 5)

orig.idents <- read.table('bcdf_test_ivhivslo.all.txt',h=T)
colnames(orig.idents)[colnames(orig.idents) == "barcode"] <- "CB"
merged <- merge(test2, orig.idents)
merged$CB <- as.factor(merged$CB)

merged$CB <- factor(merged$CB, levels = names(sort(table(merged$CB), decreasing = FALSE)))
fw61 <- ggplot(merged,aes(pos,CB,))+
    geom_point(size=0.3,col='blue',shape=95)+
    facet_grid(celltype~.,scales='free',space='free_y')+
    theme_classic()+
    xlim(1300000,1700000)+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(strip.background = element_blank(),strip.text = element_blank())

ann <- ggplot(merged,aes(pos,CB))+
    xlim(1300000,1700000)+
    annotate("rect", xmin=ANNOTexon[,1], xmax=ANNOTexon[,3], ymin=ANNOTexon[,2], ymax=ANNOTexon[,4])+
    geom_segment(data=ANNOTgene, aes(x=ANNOTgene[,1],xend=ANNOTgene[,3],y=ANNOTgene[,2],yend=ANNOTgene[,2]))+
    theme_classic()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

### plot also the 'density' summary ribbon
df2 <- as.data.frame(table(merged$pos))
df2$Var1 <- as.numeric(as.character(df2$Var1))
df2$Var1bin <- cut(df2$Var1, 1000)

# Number of rows per bin
bin_size <- 100
# Create grouping index
groups <- ceiling(seq_along(df2$Freq) / bin_size)
# Aggregate Freq and get the 'left' Var1 value (first in each bin)
binned_df <- do.call(rbind, lapply(split(df2, groups), function(chunk) {
  data.frame(
    Var1 = as.numeric(chunk$Var1[1]),  # left value
    Freq = sum(chunk$Freq)
  )
}))

binned_df <- binned_df %>% arrange(Freq)
fw61dens <- ggplot(binned_df, aes(x = Var1, y = 1,fill=log10(Freq))) + 
  geom_rect(xmin=binned_df$Var1,xmax=binned_df$Var1+500, ymin=0,ymax=2,aes(),width=10,colour=NA)+ 
  scale_fill_batlow(reverse=T,limits=c(2,3))+
  xlim(1300000,1700000)+ theme_void()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(strip.background = element_blank(),strip.text = element_blank())+theme(legend.position = "none")

### plot all together 
plot_grid(fw61, fw61dens, ann, ncol=1,rel_heights=c(1,0.2,0.1))


plot_grid(fw61dens,fw20dens, ann, ncol=1,rel_heights=c(1,1,0.2))

plot_grid(fw61,fw20, ann, ncol=1,rel_heights=c(1,1,0.3))


################ repeating for 33
param <- ScanBamParam(what=what,tag="CB")
bam1 <- scanBam('bams/iv_33.CortexLocus.fw.bam', param=param)

bam1 <- unname(bam1) # names not useful in unlisted result
#elts <- setNames(bamWhat(param), bamWhat(param))
#lst <- lapply(elts, function(elt) .unlist(lapply(bam1, "[[", elt)))

test <- as.data.frame(do.call(cbind, bam1))
hist_frame <- bind_rows(test, .id = "indiv")
test <- as.data.frame(hist_frame$V1)
test <- test[order(test$CB),] 

### a bit of data management
# some rows have - in the cell barcode, these should be dumped
test <- test[test$CB!='-',]

## keeping all cells that had at least 10 reads (could be more stringent)
test2 <- test %>% group_by(CB) %>% filter(n() > 5)

orig.idents <- read.table('../../iv_reanalysis/bcdf.all.txt',h=T)
colnames(orig.idents)[colnames(orig.idents) == "barcodes"] <- "CB"
merged <- merge(test2, orig.idents,by.x = "CB", by.y = "barcode")
merged$CB <- as.factor(merged$CB)

merged$CB <- factor(merged$CB, levels = names(sort(table(merged$CB), decreasing = FALSE)))
fw33 <- ggplot(merged,aes(pos,CB,))+
    geom_point(size=0.3,col='blue',shape=95)+
    facet_grid(celltype~.,scales='free',space='free_y')+
    theme_classic()+
    xlim(1300000,1700000)+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(strip.background = element_blank(),strip.text = element_blank())

ann <- ggplot(merged,aes(pos,CB))+
    xlim(1300000,1700000)+
    annotate("rect", xmin=ANNOTexon[,1], xmax=ANNOTexon[,3], ymin=ANNOTexon[,2], ymax=ANNOTexon[,4])+
    geom_segment(data=ANNOTgene, aes(x=ANNOTgene[,1],xend=ANNOTgene[,3],y=ANNOTgene[,2],yend=ANNOTgene[,2]))+
    theme_classic()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())


plot_grid(fw33, ann, ncol=1,rel_heights=c(1,0.1))

plot_grid(fw61,fw20, fw33, ann, ncol=1,rel_heights=c(1,1,1,0.25))

### plot also the 'density' summary ribbon
df2 <- as.data.frame(table(merged$pos))
df2$Var1 <- as.numeric(as.character(df2$Var1))
df2$Var1bin <- cut(df2$Var1, 1000)

# Number of rows per bin
bin_size <- 100
# Create grouping index
groups <- ceiling(seq_along(df2$Freq) / bin_size)
# Aggregate Freq and get the 'left' Var1 value (first in each bin)
binned_df <- do.call(rbind, lapply(split(df2, groups), function(chunk) {
  data.frame(
    Var1 = as.numeric(chunk$Var1[1]),  # left value
    Freq = sum(chunk$Freq)
  )
}))

binned_df <- binned_df %>% arrange(Freq)
fw33dens <- ggplot(binned_df, aes(x = Var1, y = 1,fill=log10(Freq))) + 
  geom_rect(xmin=binned_df$Var1,xmax=binned_df$Var1+500, ymin=0,ymax=2,aes(),width=10,colour=NA)+ 
  scale_fill_batlow(reverse=T,limits=c(2,3))+
  xlim(1300000,1700000)+ theme_void()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(strip.background = element_blank(),strip.text = element_blank())+theme(legend.position = "none")

### plot all together 
plot_grid(fw33, fw33dens, ann, ncol=1,rel_heights=c(1,0.2,0.1))

plot_grid(fw20dens, fw33dens, fw61dens, ann, ncol=1,rel_heights=c(0.5,0.5,0.5,0.2))

################ getting ribbons for more samples. 
######## Iv_38
param <- ScanBamParam(what=what,tag="CB")
bam1 <- scanBam('bams/iv_38.CortexLocus.fw.bam', param=param)

bam1 <- unname(bam1) # names not useful in unlisted result
#elts <- setNames(bamWhat(param), bamWhat(param))
#lst <- lapply(elts, function(elt) .unlist(lapply(bam1, "[[", elt)))

test <- as.data.frame(do.call(cbind, bam1))
hist_frame <- bind_rows(test, .id = "indiv")
test <- as.data.frame(hist_frame$V1)
test <- test[order(test$CB),] 

### a bit of data management
# some rows have - in the cell barcode, these should be dumped
test <- test[test$CB!='-',]

## keeping all cells that had at least 10 reads (could be more stringent)
test2 <- test %>% group_by(CB) %>% filter(n() > 5)

orig.idents <- read.table('../../iv_reanalysis/bcdf.all.txt',h=T)
colnames(orig.idents)[colnames(orig.idents) == "barcodes"] <- "CB"
merged <- merged <- merge(test2, orig.idents, by.x = "CB", by.y = "barcode")
merged$CB <- as.factor(merged$CB)

merged$CB <- factor(merged$CB, levels = names(sort(table(merged$CB), decreasing = FALSE)))
fw38 <- ggplot(merged,aes(pos,CB,))+
    geom_point(size=0.3,col='blue',shape=95)+
    facet_grid(celltype~.,scales='free',space='free_y')+
    theme_classic()+
    xlim(1300000,1700000)+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(strip.background = element_blank(),strip.text = element_blank())

ann <- ggplot(merged,aes(pos,CB))+
    xlim(1300000,1700000)+
    annotate("rect", xmin=ANNOTexon[,1], xmax=ANNOTexon[,3], ymin=ANNOTexon[,2], ymax=ANNOTexon[,4])+
    geom_segment(data=ANNOTgene, aes(x=ANNOTgene[,1],xend=ANNOTgene[,3],y=ANNOTgene[,2],yend=ANNOTgene[,2]))+
    theme_classic()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

### plot also the 'density' summary ribbon
df2 <- as.data.frame(table(merged$pos))
df2$Var1 <- as.numeric(as.character(df2$Var1))
df2$Var1bin <- cut(df2$Var1, 1000)

# Number of rows per bin
bin_size <- 100
# Create grouping index
groups <- ceiling(seq_along(df2$Freq) / bin_size)
# Aggregate Freq and get the 'left' Var1 value (first in each bin)
binned_df <- do.call(rbind, lapply(split(df2, groups), function(chunk) {
  data.frame(
    Var1 = as.numeric(chunk$Var1[1]),  # left value
    Freq = sum(chunk$Freq)
  )
}))

binned_df <- binned_df %>% arrange(Freq)
fw38dens <- ggplot(binned_df, aes(x = Var1, y = 1,fill=log10(Freq))) + 
  geom_rect(xmin=binned_df$Var1,xmax=binned_df$Var1+500, ymin=0,ymax=2,aes(),width=10,colour=NA)+ 
  scale_fill_batlow(reverse=T,limits=c(2,3))+
  xlim(1300000,1700000)+ theme_void()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(strip.background = element_blank(),strip.text = element_blank())+theme(legend.position = "none")

### plot all together 
plot_grid(fw38, fw38dens, ann, ncol=1,rel_heights=c(1,0.2,0.1))

plot_grid(fw20dens, fw33dens, fw61dens, ann, ncol=1,rel_heights=c(0.5,0.5,0.5,0.2))


######## Iv_55
param <- ScanBamParam(what=what,tag="CB")
bam1 <- scanBam('bams/iv_55.CortexLocus.fw.bam', param=param)

bam1 <- unname(bam1) # names not useful in unlisted result
#elts <- setNames(bamWhat(param), bamWhat(param))
#lst <- lapply(elts, function(elt) .unlist(lapply(bam1, "[[", elt)))

test <- as.data.frame(do.call(cbind, bam1))
hist_frame <- bind_rows(test, .id = "indiv")
test <- as.data.frame(hist_frame$V1)
test <- test[order(test$CB),] 

### a bit of data management
# some rows have - in the cell barcode, these should be dumped
test <- test[test$CB!='-',]

## keeping all cells that had at least 10 reads (could be more stringent)
test2 <- test %>% group_by(CB) %>% filter(n() > 5)

orig.idents <- read.table('../../iv_reanalysis/bcdf.all.txt',h=T)
colnames(orig.idents)[colnames(orig.idents) == "barcodes"] <- "CB"
merged <- merged <- merge(test2, orig.idents, by.x = "CB", by.y = "barcode")
merged$CB <- as.factor(merged$CB)

merged$CB <- factor(merged$CB, levels = names(sort(table(merged$CB), decreasing = FALSE)))
fw55 <- ggplot(merged,aes(pos,CB,))+
    geom_point(size=0.3,col='blue',shape=95)+
    facet_grid(celltype~.,scales='free',space='free_y')+
    theme_classic()+
    xlim(1300000,1700000)+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(strip.background = element_blank(),strip.text = element_blank())

ann <- ggplot(merged,aes(pos,CB))+
    xlim(1300000,1700000)+
    annotate("rect", xmin=ANNOTexon[,1], xmax=ANNOTexon[,3], ymin=ANNOTexon[,2], ymax=ANNOTexon[,4])+
    geom_segment(data=ANNOTgene, aes(x=ANNOTgene[,1],xend=ANNOTgene[,3],y=ANNOTgene[,2],yend=ANNOTgene[,2]))+
    theme_classic()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

### plot also the 'density' summary ribbon
df2 <- as.data.frame(table(merged$pos))
df2$Var1 <- as.numeric(as.character(df2$Var1))
df2$Var1bin <- cut(df2$Var1, 1000)

# Number of rows per bin
bin_size <- 100
# Create grouping index
groups <- ceiling(seq_along(df2$Freq) / bin_size)
# Aggregate Freq and get the 'left' Var1 value (first in each bin)
binned_df <- do.call(rbind, lapply(split(df2, groups), function(chunk) {
  data.frame(
    Var1 = as.numeric(chunk$Var1[1]),  # left value
    Freq = sum(chunk$Freq)
  )
}))

fw55dens <- ggplot(binned_df, aes(x = Var1, y = 1,fill=Freq)) + 
  geom_tile(aes(alpha=0.4))+ 
  scale_fill_gradientn(colours=c('#bebef7','#0000FF'))+
  xlim(1300000,1700000)+ theme_classic()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(strip.background = element_blank(),strip.text = element_blank())+theme(legend.position = "none")

### plot all together 
plot_grid(fw55, fw55dens, ann, ncol=1,rel_heights=c(1,0.2,0.1))

plot_grid(fw20dens, fw33dens, fw55dens, fw61dens, ann, ncol=1,rel_heights=c(0.5,0.5,0.5,0.5,0.2))


  
### plot all together 
plot_grid(fwr30, fwr30dens, fwr30dens2, ncol=1,rel_heights=c(1,0.2,0.2))

plot_grid(fw20dens, fw33dens, fw55dens, fw61dens,fwr10dens, fwr30dens, ann, ncol=1,rel_heights=c(0.5,0.5,0.5,0.5,0.5,0.5,0.2))

######## Hmr10
param <- ScanBamParam(what=what,tag="CB")
bam1 <- scanBam('bams/Hmr_10.CortexLocus.fw.bam', param=param)

bam1 <- unname(bam1) # names not useful in unlisted result
#elts <- setNames(bamWhat(param), bamWhat(param))
#lst <- lapply(elts, function(elt) .unlist(lapply(bam1, "[[", elt)))

test <- as.data.frame(do.call(cbind, bam1))
hist_frame <- bind_rows(test, .id = "indiv")
test <- as.data.frame(hist_frame$V1)
test <- test[order(test$CB),] 

### a bit of data management
# some rows have - in the cell barcode, these should be dumped
test <- test[test$CB!='-',]

## keeping all cells that had at least 10 reads (could be more stringent)
test2 <- test %>% group_by(CB) %>% filter(n() > 5)

orig.idents <- read.table('../../iv_reanalysis/bcdf.all.txt',h=T)
colnames(orig.idents)[colnames(orig.idents) == "barcodes"] <- "CB"
merged <- merged <- merge(test2, orig.idents, by.x = "CB", by.y = "barcode")
merged$CB <- as.factor(merged$CB)

merged$CB <- factor(merged$CB, levels = names(sort(table(merged$CB), decreasing = FALSE)))
fwr10 <- ggplot(merged,aes(pos,CB,))+
    geom_point(size=0.3,col='blue',shape=95)+
    facet_grid(celltype~.,scales='free',space='free_y')+
    theme_classic()+
    xlim(1300000,1700000)+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(strip.background = element_blank(),strip.text = element_blank())

ann <- ggplot(merged,aes(pos,CB))+
    xlim(1300000,1700000)+
    annotate("rect", xmin=ANNOTexon[,1], xmax=ANNOTexon[,3], ymin=ANNOTexon[,2], ymax=ANNOTexon[,4])+
    geom_segment(data=ANNOTgene, aes(x=ANNOTgene[,1],xend=ANNOTgene[,3],y=ANNOTgene[,2],yend=ANNOTgene[,2]))+
    theme_classic()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

### plot also the 'density' summary ribbon
df2 <- as.data.frame(table(merged$pos))
df2$Var1 <- as.numeric(as.character(df2$Var1))

# Number of rows per bin
bin_size <- 100
# Create grouping index
groups <- ceiling(seq_along(df2$Freq) / bin_size)
# Aggregate Freq and get the 'left' Var1 value (first in each bin)
binned_df <- do.call(rbind, lapply(split(df2, groups), function(chunk) {
  data.frame(
    Var1 = as.numeric(chunk$Var1[1]),  # left value
    Freq = sum(chunk$Freq)
  )
}))

binned_df <- binned_df %>% arrange(Freq)
fwr10dens <- ggplot(binned_df, aes(x = Var1, y = 1,fill=log10(Freq))) + 
  geom_rect(xmin=binned_df$Var1,xmax=binned_df$Var1+500, ymin=0,ymax=2,aes(),width=10,colour=NA)+ 
  scale_fill_batlow(reverse=T,limits=c(2,3))+
  xlim(1300000,1700000)+ theme_void()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(strip.background = element_blank(),strip.text = element_blank())+theme(legend.position = "none")

### plot all together 
plot_grid(fwr10, fwr10dens, ann, ncol=1,rel_heights=c(1,0.2,0.1))



######## Hmr15
param <- ScanBamParam(what=what,tag="CB")
bam1 <- scanBam('bams/Hmr_15.CortexLocus.fw.bam', param=param)

bam1 <- unname(bam1) # names not useful in unlisted result
#elts <- setNames(bamWhat(param), bamWhat(param))
#lst <- lapply(elts, function(elt) .unlist(lapply(bam1, "[[", elt)))

test <- as.data.frame(do.call(cbind, bam1))
hist_frame <- bind_rows(test, .id = "indiv")
test <- as.data.frame(hist_frame$V1)
test <- test[order(test$CB),] 

### a bit of data management
# some rows have - in the cell barcode, these should be dumped
test <- test[test$CB!='-',]

## keeping all cells that had at least 10 reads (could be more stringent)
test2 <- test %>% group_by(CB) %>% filter(n() > 5)

orig.idents <- read.table('bcdf_test_ivhivslo.all.txt',h=T)
colnames(orig.idents)[colnames(orig.idents) == "barcodes"] <- "CB"
merged <- merge(test2, orig.idents, by.x = "CB",by.y='barcode')
merged$CB <- as.factor(merged$CB)

merged$CB <- factor(merged$CB, levels = names(sort(table(merged$CB), decreasing = FALSE)))
fwr15 <- ggplot(merged,aes(pos,CB,))+
    geom_point(size=0.3,col='blue',shape=95)+
    facet_grid(celltype~.,scales='free',space='free_y')+
    theme_classic()+
    xlim(1300000,1700000)+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(strip.background = element_blank(),strip.text = element_blank())

ann <- ggplot(merged,aes(pos,CB))+
    xlim(1300000,1700000)+
    annotate("rect", xmin=ANNOTexon[,1], xmax=ANNOTexon[,3], ymin=ANNOTexon[,2], ymax=ANNOTexon[,4])+
    geom_segment(data=ANNOTgene, aes(x=ANNOTgene[,1],xend=ANNOTgene[,3],y=ANNOTgene[,2],yend=ANNOTgene[,2]))+
    theme_classic()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

### plot also the 'density' summary ribbon
df2 <- as.data.frame(table(merged$pos))
df2$Var1 <- as.numeric(as.character(df2$Var1))

# Number of rows per bin
bin_size <- 100
# Create grouping index
groups <- ceiling(seq_along(df2$Freq) / bin_size)
# Aggregate Freq and get the 'left' Var1 value (first in each bin)
binned_df <- do.call(rbind, lapply(split(df2, groups), function(chunk) {
  data.frame(
    Var1 = as.numeric(chunk$Var1[1]),  # left value
    Freq = sum(chunk$Freq)
  )
}))

binned_df <- binned_df %>% arrange(Freq)
fwr15dens <- ggplot(binned_df, aes(x = Var1, y = 1,fill=log10(Freq))) + 
  geom_rect(xmin=binned_df$Var1,xmax=binned_df$Var1+500, ymin=0,ymax=2,aes(),width=10,colour=NA)+ 
  scale_fill_batlow(reverse=T,limits=c(2,3))+
  xlim(1300000,1700000)+ theme_void()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(strip.background = element_blank(),strip.text = element_blank())+theme(legend.position = "none")

### plot all together 
plot_grid(fwr15, fwr15dens, ann, ncol=1,rel_heights=c(1,0.2,0.1))

plot_grid(fw20dens, fw33dens, fw55dens, fw61dens, ann, ncol=1,rel_heights=c(0.5,0.5,0.5,0.5,0.2))

######## Hmr20
param <- ScanBamParam(what=what,tag="CB")
bam1 <- scanBam('bams/Hmr_20.CortexLocus.fw.bam', param=param)

bam1 <- unname(bam1) # names not useful in unlisted result
#elts <- setNames(bamWhat(param), bamWhat(param))
#lst <- lapply(elts, function(elt) .unlist(lapply(bam1, "[[", elt)))

test <- as.data.frame(do.call(cbind, bam1))
hist_frame <- bind_rows(test, .id = "indiv")
test <- as.data.frame(hist_frame$V1)
test <- test[order(test$CB),] 

### a bit of data management
# some rows have - in the cell barcode, these should be dumped
test <- test[test$CB!='-',]

## keeping all cells that had at least 10 reads (could be more stringent)
test2 <- test %>% group_by(CB) %>% filter(n() > 5)

orig.idents <- read.table('bcdf_test_ivhivslo.all.txt',h=T)
colnames(orig.idents)[colnames(orig.idents) == "barcodes"] <- "CB"
merged <- merged <- merge(test2, orig.idents, by.x = "CB", by.y = "barcode")
merged$CB <- as.factor(merged$CB)

merged$CB <- factor(merged$CB, levels = names(sort(table(merged$CB), decreasing = FALSE)))
fwr20 <- ggplot(merged,aes(pos,CB,))+
    geom_point(size=0.3,col='blue',shape=95)+
    facet_grid(celltype~.,scales='free',space='free_y')+
    theme_classic()+
    xlim(1300000,1700000)+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(strip.background = element_blank(),strip.text = element_blank())

ann <- ggplot(merged,aes(pos,CB))+
    xlim(1300000,1700000)+
    annotate("rect", xmin=ANNOTexon[,1], xmax=ANNOTexon[,3], ymin=ANNOTexon[,2], ymax=ANNOTexon[,4])+
    geom_segment(data=ANNOTgene, aes(x=ANNOTgene[,1],xend=ANNOTgene[,3],y=ANNOTgene[,2],yend=ANNOTgene[,2]))+
    theme_classic()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

### plot also the 'density' summary ribbon
df2 <- as.data.frame(table(merged$pos))
df2$Var1 <- as.numeric(as.character(df2$Var1))

# Number of rows per bin
bin_size <- 100
# Create grouping index
groups <- ceiling(seq_along(df2$Freq) / bin_size)
# Aggregate Freq and get the 'left' Var1 value (first in each bin)
binned_df <- do.call(rbind, lapply(split(df2, groups), function(chunk) {
  data.frame(
    Var1 = as.numeric(chunk$Var1[1]),  # left value
    Freq = sum(chunk$Freq)
  )
}))

binned_df <- binned_df %>% arrange(Freq)
fwr20dens <- ggplot(binned_df, aes(x = Var1, y = 1,fill=log10(Freq))) + 
  geom_rect(xmin=binned_df$Var1,xmax=binned_df$Var1+500, ymin=0,ymax=2,aes(),width=10,colour=NA)+ 
  scale_fill_batlow(reverse=T,limits=c(2,3))+
  xlim(1300000,1700000)+ theme_void()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(strip.background = element_blank(),strip.text = element_blank())+theme(legend.position = "none")

### plot all together 
plot_grid(fwr20, fwr20dens, ann, ncol=1,rel_heights=c(1,0.2,0.1))


######## Hmr25
param <- ScanBamParam(what=what,tag="CB")
bam1 <- scanBam('bams/Hmr_25.CortexLocus.fw.bam', param=param)

bam1 <- unname(bam1) # names not useful in unlisted result
#elts <- setNames(bamWhat(param), bamWhat(param))
#lst <- lapply(elts, function(elt) .unlist(lapply(bam1, "[[", elt)))

test <- as.data.frame(do.call(cbind, bam1))
hist_frame <- bind_rows(test, .id = "indiv")
test <- as.data.frame(hist_frame$V1)
test <- test[order(test$CB),] 

### a bit of data management
# some rows have - in the cell barcode, these should be dumped
test <- test[test$CB!='-',]

## keeping all cells that had at least 10 reads (could be more stringent)
test2 <- test %>% group_by(CB) %>% filter(n() > 5)

orig.idents <- read.table('bcdf_test_ivhivslo.all.txt',h=T)
colnames(orig.idents)[colnames(orig.idents) == "barcodes"] <- "CB"
merged <- merged <- merge(test2, orig.idents, by.x = "CB", by.y = "barcode")
merged$CB <- as.factor(merged$CB)

merged$CB <- factor(merged$CB, levels = names(sort(table(merged$CB), decreasing = FALSE)))
fwr25 <- ggplot(merged,aes(pos,CB,))+
    geom_point(size=0.3,col='blue',shape=95)+
    facet_grid(celltype~.,scales='free',space='free_y')+
    theme_classic()+
    xlim(1300000,1700000)+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(strip.background = element_blank(),strip.text = element_blank())

ann <- ggplot(merged,aes(pos,CB))+
    xlim(1300000,1700000)+
    annotate("rect", xmin=ANNOTexon[,1], xmax=ANNOTexon[,3], ymin=ANNOTexon[,2], ymax=ANNOTexon[,4])+
    geom_segment(data=ANNOTgene, aes(x=ANNOTgene[,1],xend=ANNOTgene[,3],y=ANNOTgene[,2],yend=ANNOTgene[,2]))+
    theme_classic()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

### plot also the 'density' summary ribbon
df2 <- as.data.frame(table(merged$pos))
df2$Var1 <- as.numeric(as.character(df2$Var1))

# Number of rows per bin
bin_size <- 100
# Create grouping index
groups <- ceiling(seq_along(df2$Freq) / bin_size)
# Aggregate Freq and get the 'left' Var1 value (first in each bin)
binned_df <- do.call(rbind, lapply(split(df2, groups), function(chunk) {
  data.frame(
    Var1 = as.numeric(chunk$Var1[1]),  # left value
    Freq = sum(chunk$Freq)
  )
}))

binned_df <- binned_df %>% arrange(Freq)
fwr25dens <- ggplot(binned_df, aes(x = Var1, y = 1,fill=log10(Freq))) + 
  geom_rect(xmin=binned_df$Var1,xmax=binned_df$Var1+500, ymin=0,ymax=2,aes(),width=10,colour=NA)+ 
  scale_fill_batlow(reverse=T,limits=c(2,3))+
  xlim(1300000,1700000)+ theme_void()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(strip.background = element_blank(),strip.text = element_blank())+theme(legend.position = "none")

### plot all together 
plot_grid(fwr25, fwr25dens, ann, ncol=1,rel_heights=c(1,0.2,0.1))

######## Hmr30
param <- ScanBamParam(what=what,tag="CB")
bam1 <- scanBam('bams/Hmr_30.CortexLocus.fw.bam', param=param)

bam1 <- unname(bam1) # names not useful in unlisted result
#elts <- setNames(bamWhat(param), bamWhat(param))
#lst <- lapply(elts, function(elt) .unlist(lapply(bam1, "[[", elt)))

test <- as.data.frame(do.call(cbind, bam1))
hist_frame <- bind_rows(test, .id = "indiv")
test <- as.data.frame(hist_frame$V1)
test <- test[order(test$CB),] 

### a bit of data management
# some rows have - in the cell barcode, these should be dumped
test <- test[test$CB!='-',]

## keeping all cells that had at least 10 reads (could be more stringent)
test2 <- test %>% group_by(CB) %>% filter(n() > 5)

orig.idents <- read.table('../../iv_reanalysis/bcdf.all.txt',h=T)
colnames(orig.idents)[colnames(orig.idents) == "barcodes"] <- "CB"
merged <- merged <- merge(test2, orig.idents, by.x = "CB", by.y = "barcode")
merged$CB <- as.factor(merged$CB)

merged$CB <- factor(merged$CB, levels = names(sort(table(merged$CB), decreasing = FALSE)))
fwr30 <- ggplot(merged,aes(pos,CB,))+
    geom_point(size=0.3,col='blue',shape=95)+
    facet_grid(celltype~.,scales='free',space='free_y')+
    theme_classic()+
    xlim(1300000,1700000)+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(strip.background = element_blank(),strip.text = element_blank())

ann <- ggplot(merged,aes(pos,CB))+
    xlim(1300000,1700000)+
    annotate("rect", xmin=ANNOTexon[,1], xmax=ANNOTexon[,3], ymin=ANNOTexon[,2], ymax=ANNOTexon[,4])+
    geom_segment(data=ANNOTgene, aes(x=ANNOTgene[,1],xend=ANNOTgene[,3],y=ANNOTgene[,2],yend=ANNOTgene[,2]))+
    theme_classic()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

### plot also the 'density' summary ribbon
df2 <- as.data.frame(table(merged$pos))
df2$Var1 <- as.numeric(as.character(df2$Var1))

# Number of rows per bin
bin_size <- 100
# Create grouping index
groups <- ceiling(seq_along(df2$Freq) / bin_size)
# Aggregate Freq and get the 'left' Var1 value (first in each bin)
binned_df <- do.call(rbind, lapply(split(df2, groups), function(chunk) {
  data.frame(
    Var1 = as.numeric(chunk$Var1[1]),  # left value
    Freq = sum(chunk$Freq)
  )
}))

binned_df <- binned_df %>% arrange(Freq)
fwr30dens <- ggplot(binned_df, aes(x = Var1, y = 1,fill=log10(Freq))) + 
  geom_rect(xmin=binned_df$Var1,xmax=binned_df$Var1+500, ymin=0,ymax=2,aes(),width=10,colour=NA)+ 
  scale_fill_batlow(reverse=T,limits=c(2,3))+
  xlim(1300000,1700000)+ theme_void()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(strip.background = element_blank(),strip.text = element_blank())+theme(legend.position = "none")

### plot all together 
plot_grid(fwr30, fwr30dens, ann, ncol=1,rel_heights=c(1,0.2,0.1))

######## Hmr35
param <- ScanBamParam(what=what,tag="CB")
bam1 <- scanBam('bams/Hmr_35.CortexLocus.fw.bam', param=param)

bam1 <- unname(bam1) # names not useful in unlisted result
#elts <- setNames(bamWhat(param), bamWhat(param))
#lst <- lapply(elts, function(elt) .unlist(lapply(bam1, "[[", elt)))

test <- as.data.frame(do.call(cbind, bam1))
hist_frame <- bind_rows(test, .id = "indiv")
test <- as.data.frame(hist_frame$V1)
test <- test[order(test$CB),] 

### a bit of data management
# some rows have - in the cell barcode, these should be dumped
test <- test[test$CB!='-',]

## keeping all cells that had at least 10 reads (could be more stringent)
test2 <- test %>% group_by(CB) %>% filter(n() > 5)

orig.idents <- read.table('../../iv_reanalysis/bcdf.all.txt',h=T)
colnames(orig.idents)[colnames(orig.idents) == "barcodes"] <- "CB"
merged <- merged <- merge(test2, orig.idents, by.x = "CB", by.y = "barcode")
merged$CB <- as.factor(merged$CB)

merged$CB <- factor(merged$CB, levels = names(sort(table(merged$CB), decreasing = FALSE)))
fwr35 <- ggplot(merged,aes(pos,CB,))+
    geom_point(size=0.3,col='blue',shape=95)+
    facet_grid(celltype~.,scales='free',space='free_y')+
    theme_classic()+
    xlim(1300000,1700000)+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(strip.background = element_blank(),strip.text = element_blank())

ann <- ggplot(merged,aes(pos,CB))+
    xlim(1300000,1700000)+
    annotate("rect", xmin=ANNOTexon[,1], xmax=ANNOTexon[,3], ymin=ANNOTexon[,2], ymax=ANNOTexon[,4])+
    geom_segment(data=ANNOTgene, aes(x=ANNOTgene[,1],xend=ANNOTgene[,3],y=ANNOTgene[,2],yend=ANNOTgene[,2]))+
    theme_classic()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

### plot also the 'density' summary ribbon
df2 <- as.data.frame(table(merged$pos))
df2$Var1 <- as.numeric(as.character(df2$Var1))

# Number of rows per bin
bin_size <- 100
# Create grouping index
groups <- ceiling(seq_along(df2$Freq) / bin_size)
# Aggregate Freq and get the 'left' Var1 value (first in each bin)
binned_df <- do.call(rbind, lapply(split(df2, groups), function(chunk) {
  data.frame(
    Var1 = as.numeric(chunk$Var1[1]),  # left value
    Freq = sum(chunk$Freq)
  )
}))

binned_df <- binned_df %>% arrange(Freq)
fwr35dens <- ggplot(binned_df, aes(x = Var1, y = 1,fill=log10(Freq))) + 
  geom_rect(xmin=binned_df$Var1,xmax=binned_df$Var1+500, ymin=0,ymax=2,aes(),width=10,colour=NA)+ 
  scale_fill_batlow(reverse=T,limits=c(2,3))+
  xlim(1300000,1700000)+ theme_void()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(strip.background = element_blank(),strip.text = element_blank())+theme(legend.position = "none")

### plot all together 
plot_grid(fwr35, fwr35dens, ann, ncol=1,rel_heights=c(1,0.2,0.1))

######## Hmr40
param <- ScanBamParam(what=what,tag="CB")
bam1 <- scanBam('bams/Hmr_40.CortexLocus.fw.bam', param=param)

bam1 <- unname(bam1) # names not useful in unlisted result
#elts <- setNames(bamWhat(param), bamWhat(param))
#lst <- lapply(elts, function(elt) .unlist(lapply(bam1, "[[", elt)))

test <- as.data.frame(do.call(cbind, bam1))
hist_frame <- bind_rows(test, .id = "indiv")
test <- as.data.frame(hist_frame$V1)
test <- test[order(test$CB),] 

### a bit of data management
# some rows have - in the cell barcode, these should be dumped
test <- test[test$CB!='-',]

## keeping all cells that had at least 10 reads (could be more stringent)
test2 <- test %>% group_by(CB) %>% filter(n() > 5)

orig.idents <- read.table('../../iv_reanalysis/bcdf.all.txt',h=T)
colnames(orig.idents)[colnames(orig.idents) == "barcodes"] <- "CB"
merged <- merged <- merge(test2, orig.idents, by.x = "CB", by.y = "barcode")
merged$CB <- as.factor(merged$CB)

merged$CB <- factor(merged$CB, levels = names(sort(table(merged$CB), decreasing = FALSE)))
fwr40 <- ggplot(merged,aes(pos,CB,))+
    geom_point(size=0.3,col='blue',shape=95)+
    facet_grid(celltype~.,scales='free',space='free_y')+
    theme_classic()+
    xlim(1300000,1700000)+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(strip.background = element_blank(),strip.text = element_blank())

ann <- ggplot(merged,aes(pos,CB))+
    xlim(1300000,1700000)+
    annotate("rect", xmin=ANNOTexon[,1], xmax=ANNOTexon[,3], ymin=ANNOTexon[,2], ymax=ANNOTexon[,4])+
    geom_segment(data=ANNOTgene, aes(x=ANNOTgene[,1],xend=ANNOTgene[,3],y=ANNOTgene[,2],yend=ANNOTgene[,2]))+
    theme_classic()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

### plot also the 'density' summary ribbon
df2 <- as.data.frame(table(merged$pos))
df2$Var1 <- as.numeric(as.character(df2$Var1))

# Number of rows per bin
bin_size <- 100
# Create grouping index
groups <- ceiling(seq_along(df2$Freq) / bin_size)
# Aggregate Freq and get the 'left' Var1 value (first in each bin)
binned_df <- do.call(rbind, lapply(split(df2, groups), function(chunk) {
  data.frame(
    Var1 = as.numeric(chunk$Var1[1]),  # left value
    Freq = sum(chunk$Freq)
  )
}))

binned_df <- binned_df %>% arrange(Freq)
fwr40dens <- ggplot(binned_df, aes(x = Var1, y = 1,fill=log10(Freq))) + 
  geom_rect(xmin=binned_df$Var1,xmax=binned_df$Var1+500, ymin=0,ymax=2,aes(),width=10,colour=NA)+ 
  scale_fill_batlow(reverse=T,limits=c(2,3))+
  xlim(1300000,1700000)+ theme_void()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(strip.background = element_blank(),strip.text = element_blank())+theme(legend.position = "none")

### plot all together 
plot_grid(fwr40, fwr40dens, ann, ncol=1,rel_heights=c(1,0.2,0.1))


######## Hmr45
param <- ScanBamParam(what=what,tag="CB")
bam1 <- scanBam('bams/Hmr_45.CortexLocus.fw.bam', param=param)

bam1 <- unname(bam1) # names not useful in unlisted result
#elts <- setNames(bamWhat(param), bamWhat(param))
#lst <- lapply(elts, function(elt) .unlist(lapply(bam1, "[[", elt)))

test <- as.data.frame(do.call(cbind, bam1))
hist_frame <- bind_rows(test, .id = "indiv")
test <- as.data.frame(hist_frame$V1)
test <- test[order(test$CB),] 

### a bit of data management
# some rows have - in the cell barcode, these should be dumped
test <- test[test$CB!='-',]

## keeping all cells that had at least 10 reads (could be more stringent)
test2 <- test %>% group_by(CB) %>% filter(n() > 5)

orig.idents <- read.table('../../iv_reanalysis/bcdf.all.txt',h=T)
colnames(orig.idents)[colnames(orig.idents) == "barcodes"] <- "CB"
merged <- merged <- merge(test2, orig.idents, by.x = "CB", by.y = "barcode")
merged$CB <- as.factor(merged$CB)

merged$CB <- factor(merged$CB, levels = names(sort(table(merged$CB), decreasing = FALSE)))
fwr45 <- ggplot(merged,aes(pos,CB,))+
    geom_point(size=0.3,col='blue',shape=95)+
    facet_grid(celltype~.,scales='free',space='free_y')+
    theme_classic()+
    xlim(1300000,1700000)+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(strip.background = element_blank(),strip.text = element_blank())

ann <- ggplot(merged,aes(pos,CB))+
    xlim(1300000,1700000)+
    annotate("rect", xmin=ANNOTexon[,1], xmax=ANNOTexon[,3], ymin=ANNOTexon[,2], ymax=ANNOTexon[,4])+
    geom_segment(data=ANNOTgene, aes(x=ANNOTgene[,1],xend=ANNOTgene[,3],y=ANNOTgene[,2],yend=ANNOTgene[,2]))+
    theme_classic()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

### plot also the 'density' summary ribbon
df2 <- as.data.frame(table(merged$pos))
df2$Var1 <- as.numeric(as.character(df2$Var1))

# Number of rows per bin
bin_size <- 100
# Create grouping index
groups <- ceiling(seq_along(df2$Freq) / bin_size)
# Aggregate Freq and get the 'left' Var1 value (first in each bin)
binned_df <- do.call(rbind, lapply(split(df2, groups), function(chunk) {
  data.frame(
    Var1 = as.numeric(chunk$Var1[1]),  # left value
    Freq = sum(chunk$Freq)
  )
}))

binned_df <- binned_df %>% arrange(Freq)
fwr45dens <- ggplot(binned_df, aes(x = Var1, y = 1,fill=log10(Freq))) + 
  geom_rect(xmin=binned_df$Var1,xmax=binned_df$Var1+500, ymin=0,ymax=2,aes(),width=10,colour=NA)+ 
  scale_fill_batlow(reverse=T,limits=c(2,3))+
  xlim(1300000,1700000)+ theme_void()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(strip.background = element_blank(),strip.text = element_blank())+theme(legend.position = "none")

### plot all together 
plot_grid(fwr45, fwr45dens, ann, ncol=1,rel_heights=c(1,0.2,0.1))

######## Hmr50
param <- ScanBamParam(what=what,tag="CB")
bam1 <- scanBam('bams/Hmr_50.CortexLocus.fw.bam', param=param)

bam1 <- unname(bam1) # names not useful in unlisted result
#elts <- setNames(bamWhat(param), bamWhat(param))
#lst <- lapply(elts, function(elt) .unlist(lapply(bam1, "[[", elt)))

test <- as.data.frame(do.call(cbind, bam1))
hist_frame <- bind_rows(test, .id = "indiv")
test <- as.data.frame(hist_frame$V1)
test <- test[order(test$CB),] 

### a bit of data management
# some rows have - in the cell barcode, these should be dumped
test <- test[test$CB!='-',]

## keeping all cells that had at least 10 reads (could be more stringent)
test2 <- test %>% group_by(CB) %>% filter(n() > 5)

orig.idents <- read.table('../../iv_reanalysis/bcdf.all.txt',h=T)
colnames(orig.idents)[colnames(orig.idents) == "barcodes"] <- "CB"
merged <- merged <- merge(test2, orig.idents, by.x = "CB", by.y = "barcode")
merged$CB <- as.factor(merged$CB)

merged$CB <- factor(merged$CB, levels = names(sort(table(merged$CB), decreasing = FALSE)))
fwr50 <- ggplot(merged,aes(pos,CB,))+
    geom_point(size=0.3,col='blue',shape=95)+
    facet_grid(celltype~.,scales='free',space='free_y')+
    theme_classic()+
    xlim(1300000,1700000)+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(strip.background = element_blank(),strip.text = element_blank())

ann <- ggplot(merged,aes(pos,CB))+
    xlim(1300000,1700000)+
    annotate("rect", xmin=ANNOTexon[,1], xmax=ANNOTexon[,3], ymin=ANNOTexon[,2], ymax=ANNOTexon[,4])+
    geom_segment(data=ANNOTgene, aes(x=ANNOTgene[,1],xend=ANNOTgene[,3],y=ANNOTgene[,2],yend=ANNOTgene[,2]))+
    theme_classic()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

### plot also the 'density' summary ribbon
df2 <- as.data.frame(table(merged$pos))
df2$Var1 <- as.numeric(as.character(df2$Var1))

# Number of rows per bin
bin_size <- 100
# Create grouping index
groups <- ceiling(seq_along(df2$Freq) / bin_size)
# Aggregate Freq and get the 'left' Var1 value (first in each bin)
binned_df <- do.call(rbind, lapply(split(df2, groups), function(chunk) {
  data.frame(
    Var1 = as.numeric(chunk$Var1[1]),  # left value
    Freq = sum(chunk$Freq)
  )
}))

binned_df <- binned_df %>% arrange(Freq)
fwr50dens <- ggplot(binned_df, aes(x = Var1, y = 1,fill=log10(Freq))) + 
  geom_rect(xmin=1300000,xmax=1700000,ymin=0.975,ymax=1.025, fill='#F5F5F5')+
  geom_rect(xmin=binned_df$Var1,xmax=binned_df$Var1+500, ymin=0.975,ymax=1.025,aes(),width=10,colour=NA)+ 
  scale_fill_lapaz(reverse=T,limits=c(2,3))+
  xlim(1300000,1700000)+ theme_void()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(strip.background = element_blank(),strip.text = element_blank())+theme(legend.position = "none")

### plot all together 
plot_grid(fwr50, fwr50dens, ann, ncol=1,rel_heights=c(1,0.2,0.1))


## getting appropriate height values for each plot 


Hmm_10 Hmm_25 Hmm_40 Hmr_10 Hmr_15 Hmr_20 Hmr_25 Hmr_30 Hmr_35 Hmr_40 Hmr_45 
   639    580    265    306    571    100    641    579    169    403    110 
Hmr_50  iv_20  iv_33  iv_38  iv_55  iv_61 
   165    708    797    280    786   1232 

plot_grid(fwr10,fwr15,fwr20,fwr25,fwr30,fwr35,fwr40,fwr45,fwr50, ann, ncol=1,rel_heights=c(0.306,0.571,0.1,0.641,0.579,0.169,0.403,0.11,0.16,0.01))


plot_grid(fwr10,fwr15,fwr20,fwr25,fwr30,fwr35,fwr40,fwr45,fwr50, ann, ncol=1,rel_heights=c(0.2,.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.1))

dev.copy(pdf, 'single_reads.pdf', width=10,height=50)
dev.off()

plot_grid(fwr20,fwr25,fwr30,fwr35, ann, ncol=1,rel_heights=c(0.2,0.2,0.2,0.2,0.1))


plot_grid(fwr10dens,fwr15dens,fwr20dens,fwr25dens,fwr30dens,fwr35dens,fwr40dens,fwr45dens,fwr50dens,fw20dens,fw33dens,fw61dens, ann, ncol=1,rel_heights=c(0.2,.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2))


###################################
### LOOP TO DO THIS IN ONE SHOT ###
###################################
## ggplot code for annotation 
## read in your gff
ANNOT <- read.table("Hmel_cortexIvory_region_luca-reannotated.gff",sep="\t")
names(ANNOT) <- c("contig", "HGC", "type", "con_start", "con_end", "dot", "str", "unk", "descr")

## this for-loop makes a list of components that will be used to draw gene lines
# + orientation genes will go above - orientation genes
ANNOTgene <- data.frame(xmin=NA,xmax=NA,ymin=NA,ymax=NA)
# gene loop 
for (g in 1:nrow(ANNOT)){
  if (ANNOT$type[g] == "gene" && ANNOT$str[g] == "-")
    ANNOTgene[g,] <- c(ANNOT$con_start[g], bot+(abs(top+bot))*0.25, ANNOT$con_end[g], bot+(abs(top+bot))*0.25)
  if (ANNOT$type[g] == "gene" && ANNOT$str[g] == "+")
    ANNOTgene[g,] <- c(ANNOT$con_start[g], bot+(abs(top+bot))*0.75, ANNOT$con_end[g], bot+(abs(top+bot))*0.75)
  ANNOTgene <- na.omit(ANNOTgene)
}

## this for-loop makes a list of components that will be used to draw exon rectangles
# again, + orientation genes will go above - orientation genes
ANNOTexon <- data.frame(xmin=NA,xmax=NA,ymin=NA,ymax=NA)
for (g in 1:nrow(ANNOT)){
  if (ANNOT$type[g] == "exon" && ANNOT$str[g] == "-")
    ANNOTexon[g,] <- c(ANNOT$con_start[g], bot+(abs(top+bot))*0.1, ANNOT$con_end[g], bot+(abs(top+bot))*0.4)
  if (ANNOT$type[g] == "exon" && ANNOT$str[g] == "+")
    ANNOTexon[g,] <- c(ANNOT$con_start[g], bot+(abs(top+bot))*0.6, ANNOT$con_end[g], bot+(abs(top+bot))*0.9)
  ANNOTexon <- na.omit(ANNOTexon)
}

ann <- ggplot(merged,aes(pos,CB))+
    xlim(1300000,1700000)+ ylim(-34,4)+
    geom_segment(data=ANNOTgene, aes(x=ANNOTgene[,1],xend=ANNOTgene[,3],y=ANNOTgene[,2],yend=ANNOTgene[,2]))+
    annotate("rect", xmin=ANNOTexon[,1], xmax=ANNOTexon[,3], ymin=ANNOTexon[,2], ymax=ANNOTexon[,4])+
    annotate("text", x=1430470, y=-2, label='Ivory', size=8,hjust = 0)+
    annotate("text", x=1626952, y=-2, label='LMTK', size=8,hjust = 0)+
    annotate("text", x=1569751, y=-2, label='*', size=8,hjust = 0.5)+
    annotate("text", x=1571948, y=-2, label='*', size=8,hjust = 0)+
    annotate("text", x=1532838, y=-28, label='Cortex', size=8,hjust = 1)+
    theme_void()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
annClass <- ggplot(merged,aes(pos,CB))+
    xlim(1300000,1700000)+ ylim(-34,4)+
    geom_segment(data=ANNOTgene, aes(x=ANNOTgene[,1],xend=ANNOTgene[,3],y=ANNOTgene[,2],yend=ANNOTgene[,2]))+
    annotate("rect", xmin=ANNOTexon[,1], xmax=ANNOTexon[,3], ymin=ANNOTexon[,2], ymax=ANNOTexon[,4])+
    annotate("text", x=1430470, y=-2, label='Ivory', size=8,hjust = 0)+
    annotate("text", x=1626952, y=-2, label='LMTK', size=8,hjust = 0)+
    annotate("text", x=1569751, y=-2, label='*', size=8,hjust = 0.5)+
    annotate("text", x=1571948, y=-2, label='*', size=8,hjust = 0)+
    annotate("text", x=1532838, y=-28, label='Cortex', size=8,hjust = 1)+
    theme_classic()+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

samples <- c("Hmr_10",
             "Hmm_10", 
             "Hmr_15", 
             "Hmr_20", 
             "Hmr_25",
             "Hmm_25",
             "Hmr_30",
             "Hmr_35",
             "Hmr_40",
             "Hmm_40",
             "Hmr_45",
             "Hmr_50",
             "iv_20",
             "iv_33",
             "iv_61",
             "iv_55",
             "iv_38")
orig.idents <- read.table('../../iv_reanalysis/bcdf.all.txt', header = TRUE)
colnames(orig.idents)[colnames(orig.idents) == "barcodes"] <- "CB"
#, levels = names(sort(table(merged$CB),
# Loop 
plot_list <- list()
for (sample in samples) {
  
  bam_path <- paste0("bams/", sample, ".CortexLocus.fw.bam")
  param <- ScanBamParam(what = "pos", tag = "CB")
  bam <- scanBam(bam_path, param = param)
  bam <- unname(bam)

  test <- as.data.frame(do.call(cbind, bam))
  hist_frame <- bind_rows(test, .id = "indiv")
  test <- as.data.frame(hist_frame$V1)
  test <- test[order(test$CB),] 
  
  test <- test[order(test$CB),]
  test <- test[test$CB != '-',]
  test2 <- test %>% group_by(CB) %>% filter(n() > 5)
  
  merged <- merge(test2, orig.idents, by.x = "CB", by.y = "barcode")
  merged$CB <- factor(merged$CB, levels = names(sort(table(merged$CB), decreasing = FALSE)))
  # Plot 1: fwr plot
  fwr <- ggplot(merged, aes(pos, CB)) +
    rasterize(geom_point(size = 0.3, col = 'blue', shape = 95),dpi=300)+
    facet_grid(celltype ~ ., scales = 'free', space = 'free_y') +
    theme_classic() +
    xlim(1300000, 1700000) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())  
  # Plot 2: density ribbon
  df2 <- as.data.frame(table(merged$pos))
  df2$Var1 <- as.numeric(as.character(df2$Var1))
  bin_size <- 100
  groups <- ceiling(seq_along(df2$Freq) / bin_size)
  binned_df <- do.call(rbind, lapply(split(df2, groups), function(chunk) {
    data.frame(Var1 = as.numeric(chunk$Var1[1]), Freq = sum(chunk$Freq))
  }))
  binned_df <- binned_df %>% arrange(Freq)

  fwr_dens <- ggplot(binned_df, aes(x = Var1, y = 1, fill = log10(Freq))) +
    geom_rect(xmin = 1300000, xmax = 1700000, ymin = 0.975, ymax = 1.025, fill = '#F5F5F5') +
    geom_rect(xmin = binned_df$Var1, xmax = binned_df$Var1 + 500,
              ymin = 0.975, ymax = 1.025, width = 10, aes(), colour = NA) +
    scale_fill_imola(reverse = TRUE, limits = c(2,3)) +
    xlim(1300000, 1700000) + 
    theme_void() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank())+
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())
  # Store plots
# fwr <- fwr+rasterize(geom_point(size = 0.3, col = 'blue', shape = 95))
  plot_list[[sample]] <- list(fwr = fwr, dens = fwr_dens)
}


plot_grid(plot_list$Hmr_10$dens,
          plot_list$Hmm_10$dens,
          plot_list$Hmr_15$dens,
          plot_list$Hmr_20$dens,
          plot_list$Hmr_25$dens,
          plot_list$Hmm_25$dens,
          plot_list$Hmr_30$dens,
          plot_list$Hmr_35$dens,
          plot_list$Hmr_40$dens,
          plot_list$iv_20$dens,
          plot_list$iv_33$dens,
          plot_list$iv_61$dens,
          plot_list$iv_55$dens,
          plot_list$iv_38$dens,
          annClass,
          ncol=1,rel_heights=c(.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.5))

dev.copy(pdf, 'scRibbons.pdf', width=10,height=10)
dev.off()


plot_list$iv_20$fwr <- plot_list$iv_20$fwr + rasterize(geom_point(size = 0.3, col = 'blue', shape = 95))
plot_list$iv_33$fwr <- plot_list$iv_33$fwr + rasterize(geom_point(size = 0.3, col = 'blue', shape = 95))
plot_list$iv_61$fwr <- plot_list$iv_61$fwr + rasterize(geom_point(size = 0.3, col = 'blue', shape = 95))

## for some reason plotting as one object is failing. 
### switched this from cowplot to patchwork because cowplot was crapping out when trying tp plot rasterized panels 
#plot_list$iv_20$dens+
#          plot_list$iv_20$fwr+
#          plot_list$iv_33$dens+
#          plot_list$iv_33$fwr+
#          plot_list$iv_61$dens+
#          plot_list$iv_61$fwr+
# #         annClass+
#          plot_layout(ncol=1,heights=c(0.2,0.6,0.2,0.6,0.2,0.6))
#
#dev.copy(pdf, 'sc_readPlot_plusRibbons.pdf', width=10,height=10)
#dev.off()

plot_list$iv_20$fwr
dev.copy(pdf, 'iv20fwr.pdf',width=20,height=3)
dev.off()

plot_list$iv_33$fwr
dev.copy(pdf, 'iv33fwr.pdf',width=20,height=3)
dev.off()

plot_list$iv_61fwr
dev.copy(pdf, 'iv61fwr.pdf',width=20,height=3)
dev.off()

plot_list$iv_20$dens
dev.copy(pdf, 'iv20dens.pdf', width=20,height=1)
dev.off()

plot_list$iv_33$dens
dev.copy(pdf, 'iv33dens.pdf', width=20,height=1)
dev.off()

plot_list$iv_38$dens
dev.copy(pdf, 'iv38dens.pdf', width=20,height=1)
dev.off()

annClass
dev.copy(pdf, 'annotation.pdf',width=20, height=2)
dev.off()

###################################
### FOR OPTIX LOCUS  ###
###################################
samples <- c("Hmr_10",
             "Hmm_10", 
             "Hmr_15", 
             "Hmr_20", 
             "Hmr_25",
             "Hmm_25",
             "Hmr_30",
             "Hmr_35",
             "Hmr_40",
             "Hmm_40",
             "Hmr_45",
             "Hmr_50")

plot_listO <- list()

for (sample in samples) {
  
  bam_path <- paste0("bams/", sample, ".OptixLocus.fw.bam")
  param <- ScanBamParam(what = "pos", tag = "CB")
  bam <- scanBam(bam_path, param = param)
  bam <- unname(bam)

  test <- as.data.frame(do.call(cbind, bam))
  hist_frame <- bind_rows(test, .id = "indiv")
  test <- as.data.frame(hist_frame$V1)
  test <- test[order(test$CB),] 
  
  test <- test[order(test$CB),]
  test <- test[test$CB != '-',]
  test2 <- test %>% group_by(CB) %>% filter(n() > 5)
  
  merged <- merge(test2, orig.idents, by.x = "CB", by.y = "barcode")
  merged$CB <- factor(merged$CB, levels = names(sort(table(merged$CB), decreasing = FALSE)))
  # Plot 1: fwr plot
  fwr <- ggplot(merged, aes(pos, CB)) +
    geom_point(size = 0.3, col = 'blue', shape = 95) +
    facet_grid(celltype ~ ., scales = 'free', space = 'free_y') +
    theme_classic() +
    xlim(600000, 900000) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())  
  # Plot 2: density ribbon
  df2 <- as.data.frame(table(merged$pos))
  df2$Var1 <- as.numeric(as.character(df2$Var1))
  bin_size <- 100
  groups <- ceiling(seq_along(df2$Freq) / bin_size)
  binned_df <- do.call(rbind, lapply(split(df2, groups), function(chunk) {
    data.frame(Var1 = as.numeric(chunk$Var1[1]), Freq = sum(chunk$Freq))
  }))
  binned_df <- binned_df %>% arrange(Freq)

  fwr_dens <- ggplot(binned_df, aes(x = Var1, y = 1, fill = log10(Freq))) +
    geom_rect(xmin = 600000, xmax = 900000, ymin = 0.975, ymax = 1.025, fill = '#F5F5F5') +
    geom_rect(xmin = binned_df$Var1, xmax = binned_df$Var1 + 500,
              ymin = 0.975, ymax = 1.025, width = 10, aes(), colour = NA) +
    scale_fill_imola(reverse = TRUE, limits = c(2,3)) +
    xlim(600000, 900000) + 
    theme_void() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank())

  # Store plots
  plot_listO[[sample]] <- list(fwr = fwr, dens = fwr_dens)
}



plot_grid(plot_listO$Hmr_10$dens,
          plot_listO$Hmm_10$dens,
          plot_listO$Hmr_15$dens,
          plot_listO$Hmr_20$dens,
          plot_listO$Hmr_25$dens,
          plot_listO$Hmm_25$dens,
          plot_listO$Hmr_30$dens,
          plot_listO$Hmr_35$dens,
          plot_listO$Hmr_40$dens,
          plot_listO$Hmm_40$dens,
          plot_listO$Hmr_45$dens,
          plot_listO$Hmr_50$dens,
          ncol=1)
















#### experimenting 
library(khroma) 

#fwr30dens <- ggplot(binned_df, aes(x = Var1, y = 1,fill=log(Freq))) + 
#  geom_raster(aes())+ 
#  scale_fill_gradientn(colours=c('white','#0000FF'),limits=c(1,6))+
#  xlim(1300000,1700000)+ theme_classic()+
#  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
#  theme(strip.background = element_blank(),strip.text = element_blank())+theme(legend.position = "none")

binned_df
binned_df <- binned_df %>% arrange(Freq)
fwr30dens2 <- ggplot(binned_df, aes(x = Var1, y = 1,fill=log10(Freq))) + 
  geom_rect(xmin=binned_df$Var1,xmax=binned_df$Var1+500, ymin=0,ymax=2,aes(),width=10,colour=NA)+ 
  scale_fill_batlow(reverse=T)+
  xlim(1300000,1700000)+ theme_classic()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  theme(strip.background = element_blank(),strip.text = element_blank())+theme(legend.position = "none")


plot_grid(fwr30, fwr30dens, fwr30dens2,ann, ncol=1, rel_heights=c(1,0.2,0.2,0.1))






### ggplot code
##required packages
library(gtable,ggplot2,grid)

# set the axes. top and bot can be whatever you want if you're making 
# this as a separate plot, or specific values if you're putting it on top of
# an existing plot. 

start <- 1300000
end <- 1700000
top <- 0
bot <- -30

## read in your gff
ANNOT <- read.table("Hmel_cortexIvory_region_luca-reannotated.gff",sep="\t")
names(ANNOT) <- c("contig", "HGC", "type", "con_start", "con_end", "dot", "str", "unk", "descr")

## this for-loop makes a list of components that will be used to draw gene lines
# + orientation genes will go above - orientation genes
ANNOTgene <- data.frame(xmin=NA,xmax=NA,ymin=NA,ymax=NA)
# gene loop 
for (g in 1:nrow(ANNOT)){
  if (ANNOT$type[g] == "gene" && ANNOT$str[g] == "-")
    ANNOTgene[g,] <- c(ANNOT$con_start[g], bot+(abs(top+bot))*0.25, ANNOT$con_end[g], bot+(abs(top+bot))*0.25)
  if (ANNOT$type[g] == "gene" && ANNOT$str[g] == "+")
    ANNOTgene[g,] <- c(ANNOT$con_start[g], bot+(abs(top+bot))*0.75, ANNOT$con_end[g], bot+(abs(top+bot))*0.75)
  ANNOTgene <- na.omit(ANNOTgene)
}

## this for-loop makes a list of components that will be used to draw exon rectangles
# again, + orientation genes will go above - orientation genes
ANNOTexon <- data.frame(xmin=NA,xmax=NA,ymin=NA,ymax=NA)
for (g in 1:nrow(ANNOT)){
  if (ANNOT$type[g] == "exon" && ANNOT$str[g] == "-")
    ANNOTexon[g,] <- c(ANNOT$con_start[g], bot+(abs(top+bot))*0.1, ANNOT$con_end[g], bot+(abs(top+bot))*0.4)
  if (ANNOT$type[g] == "exon" && ANNOT$str[g] == "+")
    ANNOTexon[g,] <- c(ANNOT$con_start[g], bot+(abs(top+bot))*0.6, ANNOT$con_end[g], bot+(abs(top+bot))*0.9)
  ANNOTexon <- na.omit(ANNOTexon)
}

## the key points in this ggplot command are the 'annotate' and 'geom_segment' lines, 
## which can either be added into an existing plot or into an empty plot. 
ann <- ggplot(merged,aes(pos,CB))+
    xlim(1300000,1700000)+
    annotate("rect", xmin=ANNOTexon[,1], xmax=ANNOTexon[,3], ymin=ANNOTexon[,2], ymax=ANNOTexon[,4])+
    geom_segment(data=ANNOTgene, aes(x=ANNOTgene[,1],xend=ANNOTgene[,3],y=ANNOTgene[,2],yend=ANNOTgene[,2]))+
    theme_classic()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

  fwr_dens <- ggplot(binned_df, aes(x = Var1, y = 1, fill = log10(Freq))) +
    geom_rect(xmin = 1300000, xmax = 1700000, ymin = 0.975, ymax = 1.025, fill = '#F5F5F5') +
    geom_rect(xmin = binned_df$Var1, xmax = binned_df$Var1 + 500,
              ymin = 0.975, ymax = 1.025, width = 10, aes(), colour = NA) +
    scale_fill_imola(reverse = TRUE, limits = c(2,3)) +
    xlim(1300000, 1700000) + 
    theme_void() +
    theme(strip.background = element_blank(),
          strip.text = element_blank())+
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())
dev.copy(pdf,'legend.pdf', width=10,height=3)
dev.off()

