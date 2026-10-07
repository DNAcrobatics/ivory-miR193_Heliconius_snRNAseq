# read the gtf, retaining only the 'gene' entries
gtfAll <- read.table('Ceu_LiftoffCcro.gff',h=F)
gtf <- gtfAll[gtfAll$V3 == 'gene',]

# make a blank list
intervals <- c()

# for each gene line in the file, calculate the distance from the end to the start of the next 
for(i in 1:nrow(gtf)){
  intervals[i] <- gtf[i+1,4] - gtf[i,5]
}

# plot a histogram of intergenic intervals
ggplot(data.frame(intervals),aes(intervals))+
   geom_histogram(bins=200)+
   scale_x_log10()+
   geom_vline(xintercept=c(0,500,200,100), col='red')


## set overlap 
n <- 500

## updating the gff so genes are extended but dont clash. 


for(i in 2:(-1+nrow(gtf))){
 x <- gtf[i+1,4]-gtf[i,5]
  if(x > 1 & x < 1000000){
    if(gtf[i,7] == '+' & x > n) {
    gtf[i,5] <- gtf[i,5]+n
    }
    if(gtf[i,7] == '+' & x < n) {
    gtf[i,5] <- round(gtf[i,5]+x/2-1, digits=0)
    }
  }  
 y <- gtf[i,4]-gtf[i-1,5]
  if(y > 1 & y < 1000000){
    if(gtf[i,7] == '-' & y > n) {
    gtf[i,4] <- gtf[i,4]-n
    }
    if(gtf[i,7] == '-' & y < n) {
    gtf[i,4] <- round(gtf[i,4]-y/2+1, digits=0)
    }
  } 
}


## reconcatenate and sort. 
gtfAll <- gtfAll[gtfAll$V3 != 'gene',]
gtfAll <- rbind(gtf,gtfAll)
gtfAll <- gtfAll[order(gtfAll$V1, gtfAll$V4, gtfAll$V9),]

## write table out. 
write.table(gtfAll, 'gene_lines.extended500.gtf',sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# cl4 is defo sockets.
# cl5  




