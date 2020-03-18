#!/bin/R

### PLOTS ######################################################################
po <- function(cross,nme){
 sm <- scanone(cross, pheno.col=4, model="binary",method="mr")
 Y <- c(0, as.numeric(gsub(".*:","",markernames(cross))))/1000000
 X <- 1:length(Y)
 gt <- geno.table(cross)
 plot_test(nme, width = 5500, height = 750)
 par(mfrow=c(3,1))
  plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
  abline(h=5)
  plot(1:length(gt[,1]), -log10(gt[,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
  abline(h=3)
  plot(c(1,length(X)),c(0,max(Y)),type="n", ylab='physical position')
   points(X,Y)
 dev.off()
}
################################################################################
################################################################################
pop <- 'NBH'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
library(scales)

################################################################################
load(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))
################################################################################
pulls <- markernames(cross)


mpath <- '/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops'
fl <- 'NBH.um.unmapped.f2.csvr'
################################################################################
## read in the QTL cross
cross.all <- read.cross(file = file.path(mpath, fl),
format = "csvr", geno = c(1:3), estimate.map = FALSE)
################################################################################
nw_marks <- grep('NW',markernames(cross.all), value = T)

assign_nw <- pull.markers(cross.all,c(nw_marks,pulls))


#### THIS PVALUE APPEARS TO RETAIN EVEN DISTRIBUTION
gt <- geno.table(assign_nw)
toss <- rownames(gt[which(gt[,'P.value'] < 1.0e-3),])
assign_nw <- drop.markers(assign_nw,toss)
assign_nw <- est.rf(assign_nw)


RF <- 0.05
LOD <- 16
cross_NW <- formLinkageGroups(assign_nw, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
cross_NW2 <- subset(cross_NW,chr= names(which(nmar(cross_NW) > 4)))
nw_chr <- grep('NW_',markernames(cross_NW2), value = T)

nw <- file.path(mpath,'NBH_nw_test.tsv')
write.table(nw_chr,nw)


### READ IN THE CROSS
fl <- paste0(pop,'_filtered_unphased_NW.csv')
cross <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

### READ IN THE MARKER TABLE
mfh <- file.path(mpath,'NBH_high_confid.tsv')
mfl <- file.path(mpath,'NBH_low_confid.tsv')
low_confid <- as.character(read.table(mfl)[,1])
high_confid <- as.character(read.table(mfh)[,1])


cross <- pull.markers(cross.all,c(high_confid,nw_chr))


#nw_marks <- grep('NW_',markernames(cross), value = T)
#nw_chr <- grep('NW_',chrnames(cross), value = T)
#chr_marks <- markernames(cross,1:24)
#
##### HIGH CONFIDENCE MARKERS ONLY
#cross <- pull.markers(cross,high_confid)

#### THIS PVALUE APPEARS TO RETAIN EVEN DISTRIBUTION
gt <- geno.table(cross)
toss <- rownames(gt[which(gt[,'P.value'] < 1.0e-3),])
cross <- drop.markers(cross,toss)


### THIN TO LEAST DISTORTED PER KB
#cross <- thin_by_distortion(cross,10)
#cross <- thin_by_distortion(cross,100)

## cross <- thin_by_distortion(cross,50)

################################################################################
### ASSIGN UNMAPPED SCAFFOLDS TO GROUPS

### Get mean LOD of each unmapped scaffold and the CHRs

nw_chr <- grep('NW_',chrnames(cross), value = T)
chr <- 1:24

cross <- est.rf(cross)

ldm <- function(nw) {
 sapply(chr,function(z){ mean(pull.rf(cross, what='lod')[markernames(cross,nw),markernames(cross,z)]) })
}

#library(doParallel)
#cl <- makeCluster(20)
#registerDoParallel(cl)
#ld <- foreach(nw = nw_chr, .inorder = F, .packages = libs2load) %dopar% ldm(nw)
#ld <- do.call(rbind,ld)
#rownames(ld) <- nw_chr
#reassign <- apply(ld,1,which.max)
#
#nw_marks_assign <- sapply(names(reassign),markernames,cross = cross)
#nw_length <- sapply(nw_marks_assign,length)
#nw_marks_assign <- as.character(unlist(nw_marks_assign))
#nw_ch <- rep(as.numeric(reassign), times = as.numeric(nw_length))
#nw_pos <- 1:length(nw_ch)
#nw_old <- gsub(":.*","",nw_marks_assign)
#
#### REASSING THE UNMAPPED MARKERS
#move <- data.frame(cbind(nw_old,nw_marks_assign,nw_ch,nw_pos), stringsAsFactors=F)
#movefl <- file.path(mpath,'NBH_NW_scaffold_assignments.tsv')
#write.table(move,movefl)

### READ THE UNMAPPED MARKER TABLE
movefl <- file.path(mpath,'NBH_NW_scaffold_assignments.tsv')
move <- read.table(movefl, stringsAsFactors = F, header=T, sep = " ")
move <- move[which(move$nw_marks_assign %in% markernames(cross)),]

cross2 <- cross
for (i in 1:length(move[,1])){
 cross2 <<- movemarker(cross2, marker = move[i,'nw_marks_assign'], newchr = move[i,'nw_ch'], newpos = as.numeric(move[i,'nw_pos']))
}
cross2 <- subset(cross2,chr=1:24)
################################################################################

### REORDER WITH GENOTYPING ERRORS, THEN DROP HIGH DOUBLE CROSSOVERS MARKERS
cross3 <- tspOrder(cross = cross2, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
cross3 <- removeDoubleXO(cross3)

mis <- misg(cross3,0.10)
drop <- names(which(colSums(is.na(pull.geno(cross3))) > mis))
cross4 <- drop.markers(cross3,drop)
################################################################################






cross_13 <- formLinkageGroups(subset(cross,chr=13), max.rf = 0.05, min.lod = 15, reorgMarkers = TRUE)






gt <- geno.table(cross4)
thres <- sapply(1:24, function(X) {
 quantile(-log10(gt[markernames(cross4,X),'P.value']),0.999)
 #rownames(gt[markernames(cross5,X),])[ind]
})


#### REMOVE OUTLIER SEG DIST
chr <- c(1,3:
rm <- sapply(1:24, function(X) {
 ind <- which.max(-log10(gt[markernames(cross6,X),'P.value']))
 rownames(gt[markernames(cross6,X),])[ind]
})


cross5a <- drop.markers(cross5,rm)




### FINAL RE-ORDER
cross5 <- tspOrder(cross = cross4, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
################################################################################


### Drop low lod markers





### WRITE MAP
mapfile <- paste0(pop,'_',sum(nmar(cross5)),'_noimpute_NW_remapped_tsp')
filename <- file.path(mpath,mapfile)
write.cross(cross5,filestem=filename,format="csv")

cross6 <- fill.geno(cross5, method="maxmarginal", error.prob = 0.05, min.prob=0.95)
plotit(cross6)

gt <- geno.table(cross6)
rm <- sapply(1:24, function(X) {
 ind <- which.max(-log10(gt[markernames(cross6,X),'P.value']))
 rownames(gt[markernames(cross6,X),])[ind]
})
cross6a <- drop.markers(cross6,rm)

cross6 <- tspOrder(cross = cross6a, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')


### WRITE MAP
mapfile <- paste0(pop,'_',sum(nmar(cross6)),'_imputed_NW_remapped_tsp')
filename <- file.path(mpath,mapfile)
write.cross(cross6,filestem=filename,format="csv")

### PRUNE ONE MARKER THAT is DISTORTED ON THESE CHRs
gt <- geno.table(cross6)
prune.a.marker <- c(3,4,6,10,12,15,16,21,22,23,24)

rm <- sapply(prune.a.marker, function(X) {
 ind <- which.max(-log10(gt[markernames(cross6,X),'P.value']))
 rownames(gt[markernames(cross6,X),])[ind]
})

cross7 <- drop.markers(cross6,rm)

prune.a.marker <- c(3,4,9,15,16,24)

rm <- sapply(prune.a.marker, function(X) {
 ind <- which.max(-log10(gt[markernames(cross7,X),'P.value']))
 rownames(gt[markernames(cross7,X),])[ind]
})

cross8 <- drop.markers(cross7,rm)
erprob <- 0.01
cross_map8 <- est.map(cross8, error.prob=erprob, map.function="kosambi",maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE, n.cluster=5)
cross8 <- replace.map(cross8,cross_map8)

### WRITE MAP
mapfile <- paste0(pop,'_',sum(nmar(cross8)),'_imputed_NW_est_remapped_pruned_tsp')
filename <- file.path(mpath,mapfile)
write.cross(cross8,filestem=filename,format="csv")



cross8 <- cross6

erprob <- 0.01
cross9 <- fill.geno(cross8, method="no_dbl_XO", error.prob = 0.01, min.prob=0.95)

cross_map9 <- est.map(cross9, error.prob=erprob, map.function="kosambi",maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE, n.cluster=5)
cross9 <- replace.map(cross9,cross_map9)

### WRITE MAP
mapfile <- paste0(pop,'_',sum(nmar(cross9)),'_noxo_imputed_NW_est_remapped_pruned_tsp')
filename <- file.path(mpath,mapfile)
write.cross(cross9,filestem=filename,format="csv")


################################################################################
## THIN TO ONE MARKER PER KB

################################################################################
################################################################################






################################################################################
################################################################################
png(paste0('~/public_html/',pop,'all_phys.png'),height=2500,width=2500)
par(mfrow=c(6,4))
for (B in 1:24){
 fl <- paste0(pop,'_order_impute_',B,'_tsp.csv')
 cross_plot <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
 Y <- c(0, as.numeric(gsub(".*:","",markernames(cross_plot))))/1000000
 X <- 1:length(Y)
 plot(c(1,length(X)),c(0,max(Y)),type="n", xlab=paste('chr',B), ylab='physical position')
 points(X,Y)
}
dev.off()


png(paste0('~/public_html/',pop,'all_rf.png'),height=2500,width=2500)
par(mfrow=c(6,4))
for (B in 1:24){
 fl <- paste0(pop,'_order_impute_',B,'_tsp.csv')
 cross_plot <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
 plotRF(cross_plot)
}
dev.off()
#######################################################################################
#######################################################################################
