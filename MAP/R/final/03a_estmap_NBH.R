#!/bin/R

################################################################################
pop <- 'NBH'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
library(scales)

### READ IN THE CROSS
fl <- paste0(pop,'_filtered_unphased_NW.csv')
cross <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

### READ IN THE MARKER TABLE
mfh <- file.path(mpath,'NBH_high_confid.tsv')
mfl <- file.path(mpath,'NBH_low_confid.tsv')
low_confid <- as.character(read.table(mfl)[,1])
high_confid <- as.character(read.table(mfh)[,1])

nw_marks <- grep('NW_',markernames(cross), value = T)
nw_chr <- grep('NW_',chrnames(cross), value = T)
chr_marks <- markernames(cross,1:24)

#### HIGH CONFIDENCE MARKERS ONLY
##cross1 <- pull.markers(cross,c(unique(low_confid, high_confid)))
cross1 <- cross
#### FOUND LATER TO BE WRONG CROSS TYPE (SEE BELOW)
mfd <- file.path(mpath,'NBH_found_unlinked.tsv')
drop <- as.character(read.table(mfd)[,1])
cross1 <- drop.markers(cross1,drop)

mfd2 <- file.path(mpath,'NBH_found_unlinked_2.tsv')
drop2 <- as.character(read.table(mfd2)[,1])
cross1 <- drop.markers(cross1,drop2)

#### THIS PVALUE APPEARS TO RETAIN EVEN DISTRIBUTION
gt <- geno.table(cross1)
toss <- rownames(gt[which(gt[,'P.value'] < 5.5e-4),])
cross2 <- drop.markers(cross1,toss)

#### UNLINKED ##################################################################
crs <- formLinkageGroups(cross2, max.rf = 0.15, min.lod = 15, reorgMarkers = TRUE)

nmrs <- sapply(1:24, function(D) {
 crs <- formLinkageGroups(subset(cross2,chr = D), max.rf = 0.15, min.lod = 15, reorgMarkers = TRUE)
  markernames(crs,names(nmar(crs)[-1]))
})

cross3 <- drop.markers(cross2,unlist(nmrs))
################################################################################

### FIND HIGH XO MARKERS
cross_chr <- subset(cross3,1:24)

xochr <- sapply(1:24, function(chr) {
 xo <- locateXO(cross_chr, chr=chr, full.info=T)
 xo <- do.call(rbind, xo)
 xo <- xo[which(xo[,'nTypedBetween'] == 0),c('ileft','iright')]
 table(markernames(cross_chr,chr)[as.numeric(c(xo))])
})
### Remove markers that have more than 5 XO events (plotted hist for cuttoff
xochr <- unlist(xochr)[which(as.numeric(unlist(xochr)) > 5)]
cross4 <- drop.markers(cross3,names(xochr))
################################################################################

### THIN TO LEAST DISTORTED PER 5KB
#cross <- thin_by_distortion(cross,10)
#cross <- thin_by_distortion(cross,100)
cross4a <- thin_by_distortion(cross4,50)

### READ THE UNMAPPED MARKER ASSIGNMENT TABLE
movefl <- file.path(mpath,'NBH_NW_scaffold_assignments.tsv')
move <- read.table(movefl, stringsAsFactors = F, header=T, sep = " ")
move <- move[which(move$nw_marks_assign %in% markernames(cross4)),]

### ASSIGN UNMAPPED MARKERS
cross5 <- cross4a
for (i in 1:length(move[,1])){
 cross5 <<- movemarker(cross5, marker = move[i,'nw_marks_assign'], newchr = move[i,'nw_ch'], newpos = as.numeric(move[i,'nw_pos']))
}
cross5 <- subset(cross5,chr=1:24)
################################################################################

### CLEANUP REMOVE HIGHEST DISTORTED FROM EACH LG
gt <- geno.table(cross5)
mark.thres <- sapply(1:24, function(X) {
 qt <- as.numeric(quantile( -log10(gt[markernames(cross5,X),'P.value']),0.99))
 ind <- which(-log10(gt[markernames(cross5,X),'P.value']) > qt)
 rownames(gt[markernames(cross5,X),])[ind]
})

drop <- unlist(mark.thres)
cross6 <- drop.markers(cross5,drop)
################################################################################

### REORDER WITH GENOTYPING ERRORS, THEN DROP HIGH DOUBLE CROSSOVERS MARKERS
cross7 <- est.rf(cross6, maxit=1000, tol=1e-6)
cross7 <- tspOrder(cross = cross7, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
cross7 <- removeDoubleXO(cross7)
################################################################################

### FINAL RE-ORDER
################################################################################
mis <- misg(cross7,0.10)
drop <- names(which(colSums(is.na(pull.geno(cross7))) > mis))
cross8 <- drop.markers(cross7,drop)
cross8 <- tspOrder(cross = cross8, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
################################################################################

### WRITE MAP
mapfile <- paste0(pop,'_',sum(nmar(cross8)),'_noimpute_NW_remapped_tsp')
filename <- file.path(mpath,mapfile)
write.cross(cross8,filestem=filename,format="csv")
################################################################################
################################################################################

#### CLEANUP REMOVE HIGHEST DISTORTED FROM EACH CROSS
#gt <- geno.table(cross8)
#mark.thres <- sapply(1:24, function(X) {
# qt <- as.numeric(quantile( -log10(gt[markernames(cross8,X),'P.value']),0.99))
# ind <- which(-log10(gt[markernames(cross8,X),'P.value']) > qt)
# rownames(gt[markernames(cross8,X),])[ind]
#})
#names(mark.thres) <- 1:24
#drop <- unlist(mark.thres)
#cross9 <- drop.markers(cross8,drop)
#cross9 <- tspOrder(cross = cross9, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

### FIX GENOTYPING ERRORS WITH GENOPROB
cross9 <- fill.geno(cross8, method="maxmarginal", error.prob = 0.05, min.prob=0.95)
mis <- misg(cross9,0.10)
drop <- names(which(colSums(is.na(pull.geno(cross9))) > mis))
cross10 <- drop.markers(cross9,drop)

cross10 <- est.rf(cross10)
cross10 <- tspOrder(cross = cross10, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
cross10 <- fill.geno(cross10, method="no_dbl_XO", error.prob = 0.01, min.prob=0.99)

crs <- formLinkageGroups(cross10, max.rf = 0.15, min.lod = 15, reorgMarkers = TRUE)
drop <- rownames(geno.table(crs,chr=25:26))
cross11 <- drop.markers(cross10,drop)


### WRITE IMPUTED AND CORRECTED MAP
mapfile <- paste0(pop,'_',sum(nmar(cross11)),'_imputed_NW_remapped_tsp')
filename <- file.path(mpath,mapfile)
write.cross(cross11,filestem=filename,format="csv")
################################################################################
################################################################################

#### LINKAGE TEST
 rf <- est.rf(cross9)
 lod.df <- pull.rf(rf, what='lod')
 rf.df <- pull.rf(rf)
 gt <- geno.table(rf)

 ld <- sapply(1:24,function(z){
  mars <- markernames(rf,z)
  sapply(mars,function(X) { mean(lod.df[X,mars], na.rm=T) })
 })

 lodquan <- unlist(lapply(ld,quantile,probs = 0.01, na.rm = T))
 lodquan.50 <- unlist(lapply(ld,quantile,probs = 0.50, na.rm = T))
 ld <- unlist(ld)

 lr <- sapply(1:24,function(z){
  mars <- markernames(rf,z)
  sapply(mars,function(X) { mean(rf.df[X,mars], na.rm=T) })
 })
 rfquan <- unlist(lapply(lr,quantile,probs = 0.99, na.rm = T))
 rfquan.50 <- unlist(lapply(lr,quantile,probs = 0.50, na.rm = T))
 lr <- unlist(lr)

 chrnm <- gt[names(lr),'chr']

 plot_test('sdf', width=2000)
 par(mfrow=c(2,1))
  plot(1:length(ld), ld, col=as.factor(chrnm))
  plot(1:length(lr), lr, col=as.factor(chrnm))
 dev.off()

crs <- formLinkageGroups(rf, max.rf = 0.15, min.lod = 15, reorgMarkers = TRUE)
drop <- rownames(geno.table(crs,chr=25:29))
### FOUND TO BE TRICKY UNLINKED
mfd <- file.path(mpath,'NBH_found_unlinked_2.tsv')
write.table(drop,mfd)
drop_early <- as.character(read.table(mfd)[,1])

mfd <- file.path(mpath,'NBH_found_unlinked.tsv')
drop <- as.character(read.table(mfd)[,1])
cross1 <- drop.markers(cross1,drop)

cross11 <- drop.markers(cross10,drop)
cross11 <- est.rf(cross11)
cross11 <- tspOrder(cross = cross11, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
##plotit(cross11)
################################################################################

## #crs <- formLinkageGroups(rf, max.rf = 0.15, min.lod = 15, reorgMarkers = TRUE)
## THE FINAL SET IS LINKED BY 15 LOD RF < 0.15


##chr8_low_conf <- c('8:2490052','8:1966945','8:2506730','8:1998563','8:2402677','8:1899164','8:4407253')
##chr21_low_conf <- c('21:1740020','21:2111258','21:1447859','21:1373566','21:1729581','21:1412854','21:729095 ','21:770927 ','21:1925351','21:1865364','21:601555 ','21:1811166','21:1049847')
##cross11 <- drop.markers(cross10,c(chr8_low_conf,chr21_low_conf))
##cross11 <- est.rf(cross11)
##cross11 <- tspOrder(cross = cross11, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
##plotit(cross11)


#### WRITE IMPUTED AND CORRECTED MAP
mapfile <- paste0(pop,'_',sum(nmar(cross11)),'_imputed_high_confidence_tsp')
filename <- file.path(mpath,mapfile)
write.cross(cross11,filestem=filename,format="csv")
#################################################################################


source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales")
suppressMessages(sapply(libs2load, require, character.only = TRUE))

################################################################################

#fl <- 'NBH_2897_imputed_NW_remapped_tsp.csv'
fl <- 'NBH_2897_imputed_high_confidence_tsp.csv'
#fl <- file.path(mpath,fl)

cross <- read.cross(file=fl , format = "csv", dir=mpath,
 genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

################################################################################

#################################################################################

png(paste0('~/public_html/',pop,'high_conf_imputed_rf.png'),height=2500,width=2500)
par(mfrow=c(6,4))
for (B in 1:24){
  plotRF(cross11,chr=B)
}
dev.off()
################################################################################
################################################################################
plot_test('nbh_high_confid_rf')
 plotRF(cross10,zmax=10,col.scheme="redblue")
dev.off()
################################################################################
################################################################################

### ESTIMATE MAP ###############################################################

loglik <- err <- c(0.0001, 0.001, 0.01, 0.05)

update.lik <- function(z){
  cat(z, "of", length(err), "\n")
  tempmap <- est.map(cross,maxit=100, error.prob=err[z])
  loglik[z] <- sum(sapply(tempmap, attr, "loglik"))
}

library(qtl)
library(doParallel)
cl <- makeCluster(5)
registerDoParallel(cl)

loglik <- foreach(z = seq(along=err), .inorder = T, .export = c("loglik"), .packages = c("qtl")) %dopar% update.lik(z)

loglik <- unlist(loglik)

lod <- (loglik - max(loglik))/log(10)

erprob <- err[which.max(lod)]

print(paste('error lod =',erprob))

cross_map <-  est.map(cross, error.prob=erprob, map.function="kosambi",maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)

cross <- qtl:::replace.map(cross,cross_map)

direc <- sapply(1:24,function(i) {
 pos <- as.numeric(gsub(".*:","",markernames(cross,i)))
 map <- as.numeric(pull.map(cross)[[i]])
 cor(pos,map, use="complete.obs")
})

if(any(direc < 0)) cross <- flip.order(cross,which(direc < 0))

mapfile <- paste0(pop,'_',sum(nmar(cross)),'_imputed_high_confidence_tsp_mapped')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")

print(paste(filename, 'cross written'))
