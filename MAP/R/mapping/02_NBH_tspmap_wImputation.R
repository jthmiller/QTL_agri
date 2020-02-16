#!/bin/R

i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'


libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap")
suppressMessages(sapply(libs2load, require, character.only = TRUE))

library(scales)



indHET <- cbind(1:length(cross$pheno$ID),rowSums(pull.geno(cross) == 2, na.rm = T))
indHET <- indHET[order(as.numeric(indHET[,2]), decreasing=T),]

drop1 <- names(which(! pull.geno(cross)[indHET[1,1],] == 2))
drop2 <- names(which(! pull.geno(cross)[indHET[2,1],] == 2))
cross <- drop.markers(cross,c(drop1,drop2))
print(nmar(cross))


################################################################################

sm <- scanone(cross, pheno.col=4, model="binary",method="mr")
ord <- order(as.numeric(sm$lod))

cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
 maxit = 1, tol = 0.1, sex.sp = F)

crs <- fill.geno(crs, method="maxmarginal", error.prob = 0.05, min.prob=0.9975)
crs <- tspOrder(cross = crs,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

drop <- names(which(colSums(is.na(pull.geno(crs))) > (nind(crs)*0.15)))
crs <- drop.markers(crs,drop)


crs <- cross_flg
crs <- cross30
Y <- c(0, as.numeric(gsub(".*:","",markernames(crs))))/1000000
X <- 1:length(Y)
gt <- geno.table(crs)
sm <- scanone(crs, pheno.col=4, model="binary",method="mr")

png(paste0('~/public_html/',pop,'_gts_test',i,'.png'),height=1500,width=2500)
par(mfrow=c(4,1))
 plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,20), cex =1)
 plot(1:length(gt[,1]), -log10(gt[,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,18), cex =1)
 crs$pheno$gtps <- (as.numeric(rowSums(pull.geno(crs) == 1 | pull.geno(crs) == 1, na.rm = T))*10) + (as.numeric(rowSums(pull.geno(crs) == 3, na.rm = T))*5)
 #crs$pheno$gtps <- rowSums(pull.geno(cross))
 geno.image(crs, reorder=6, cex=2)
 plot(c(1,length(X)),c(0,max(Y)),type="n", xlab=paste('chr',i), ylab='physical position')
 abline(h=27.504907, col='red')
 points(X,Y)
dev.off()
################################################################################



RF <- 0.05
LOD <- 17
cross_flg <- formLinkageGroups(crs, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

chr <- as.character(i)
map <- pull.map(cross_flg)
newpos <- lapply(map,function(X) { setNames(as.numeric(gsub(".*:","",markernames(cross_flg)))/100000,markernames(cross_flg))  } )
attr(newpos,'class') <- 'map'
class(newpos[[chr]]) <- 'A'
attr(newpos[[chr]], "loglik") <- attr(map[[chr]], "loglik")
names(newpos) <- chr
cross_flg <- replace.map(cross_flg,newpos)
print(summary(pull.map(cross_flg)))

### GET ONLY 1 MARKER PER RAD TAG
mrks <- as.numeric(gsub(".*:","",markernames(cross_flg)))/100
names(mrks) <- markernames(cross_flg)
n.missing <- nmissing(cross_flg, what="mar")
wts <- -log( (n.missing+1) / (nind(cross_flg)+1) )
a <- pickMarkerSubset(mrks, 2, wts)
cross <- pull.markers(cross,a)
print(nmar(cross))






cross_flg1 <- removeDoubleXO(cross_flg, chr=1:5)
cross_flg2 <- tspOrder(cross = subset(cross_flg1, chr = 1:50),hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
cross_flg3 <- fill.geno(cross_flg2, method="maxmarginal", error.prob = 0.05, min.prob=0.99)

RF <- 0.05
LOD <- 17
cross_reflg <- formLinkageGroups(cross_flg3, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
cross_reflg2 <- tspOrder(cross = cross_reflg ,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')


plm(cross_flg)
################################################################################
################################################################################
################################################################################
################################################################################
##  8  12 23

fl <- file.path(paste0(pop,'_unmapped_filtered.csv'))
cross_all <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

for (i in 1:24){

mapfile <- paste0(pop,'_all_mark_',i,'_tsp')
filename <- file.path(mpath,mapfile)

cross  <- subset(cross_all,chr=i)
sex <- read.table(file.path(mpath,'sex.txt'),stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross$pheno$ID), 'sex']
cross$pheno$sex <- sex.vec

gt <- geno.table(cross)
sm <- scanone(cross, pheno.col=4, model="binary",method="mr")
Y <- c(0, as.numeric(gsub(".*:","",markernames(cross))))/1000000
X <- 1:length(Y)

png(paste0('~/public_html/',pop,'_gts_preclean',i,'.png'),height=1500,width=2500)
par(mfrow=c(4,1))
 plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,18), cex =1)
 plot(1:length(gt[,1]), -log10(gt[,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,18), cex =1)
 cross$pheno$gtps <- as.numeric(rowSums(pull.geno(cross) == 3 | pull.geno(cross) == 1, na.rm = T))
 geno.image(cross, chr=i, reorder=6, cex=2)
 plot(c(1,length(X)),c(0,max(Y)),type="n", xlab=paste('chr',i), ylab='physical position')
 points(X,Y)
dev.off()

################################################################################
################################################################################

## SET MAP TO RESONABLE DIST TO CLEAN
chr <- as.character(i)
map <- pull.map(cross)
newpos <- lapply(map,function(X) { setNames(as.numeric(gsub(".*:","",markernames(cross)))/100000,markernames(cross))  } )
attr(newpos,'class') <- 'map'
class(newpos[[chr]]) <- 'A'
attr(newpos[[chr]], "loglik") <- attr(map[[chr]], "loglik")
names(newpos) <- chr
cross <- replace.map(cross,newpos)
print(summary(pull.map(cross)))

### Remove single crossovers (not informative for ordering chuncks
cross <- removeDoubleXO(cross, chr=chr)
drop <- names(which(colSums(is.na(pull.geno(cross))) > (nind(cross)*0.15)))
cross <- drop.markers(cross,drop)
print(nmar(cross))
plm(cross)

### GET ONLY 1 MARKER PER RAD TAG
mrks <- as.numeric(gsub(".*:","",markernames(cross)))/100
names(mrks) <- markernames(cross)
n.missing <- nmissing(subset(cross, chr=chr), what="mar")
wts <- -log( (n.missing+1) / (nind(cross)+1) )
a <- pickMarkerSubset(mrks, 2, wts)
cross <- pull.markers(cross,a)
print(nmar(cross))

cross <- removeDoubleXO(cross, chr=chr)
drop <- names(which(colSums(is.na(pull.geno(cross))) > (nind(cross)*0.15)))
cross <- drop.markers(cross,drop)
print(nmar(cross))
plm(cross)

### Smooth over errors
cross <- fill.geno(cross, method="maxmarginal", error.prob = 0.05, min.prob=0.9975)
cross <- fill.geno(cross, method="no_dbl_XO", error.prob = 0.05, min.prob=0.9975)
plm(cross)


## PRELIM ORDER w all errors
cross <- tspOrder(cross = cross,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

pos <- as.numeric(gsub(".*:","",markernames(cross)))
map <- as.numeric(pull.map(cross)[[1]])
if(cor(pos,map, use="complete.obs") < 0) cross <- flip.order(cross, i)

### Remove single crossovers
cross <- removeDoubleXO(cross, chr=chr)
drop <- names(which(colSums(is.na(pull.geno(cross))) > (nind(cross)*0.20)))
length(drop)
cross <- drop.markers(cross,drop)
print(nmar(cross))

gt <- geno.table(cross)
sm <- scanone(cross, pheno.col=4, model="binary",method="mr")
Y <- c(0, as.numeric(gsub(".*:","",markernames(cross))))/1000000
X <- 1:length(Y)

png(paste0('~/public_html/',pop,'_gts_postclean_mapped',i,'.png'),height=1500,width=2500)
par(mfrow=c(4,1))
 plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 1)
 plot(1:length(gt[,1]), -log10(gt[,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 1)
 cross$pheno$gtps <- as.numeric(rowSums(pull.geno(cross) == 3 | pull.geno(cross) == 1, na.rm = T))
 geno.image(cross, chr=i, reorder=6, cex=2)
 plot(c(1,length(X)),c(0,max(Y)),type="n", xlab=paste('chr',i), ylab='physical position')
 points(X,Y)
dev.off()


#### MAP #######################################################################

cross <- tspOrder(cross = cross,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

pos <- as.numeric(gsub(".*:","",markernames(cross)))
map <- as.numeric(pull.map(cross)[[1]])

if(cor(pos,map, use="complete.obs") < 0){
 cross <<- flip.order(cross, i)
}

cross <- shiftmap(cross, offset=0)

write.cross(cross,chr=i,filestem=filename,format="csv")


gt <- geno.table(cross)
sm <- scanone(cross, pheno.col=4, model="binary",method="mr")
Y <- c(0, as.numeric(gsub(".*:","",markernames(cross))))/1000000
X <- 1:length(Y)

png(paste0('~/public_html/',pop,'_gts_phenosort_mapped',i,'.png'),height=1500,width=2500)
par(mfrow=c(4,1))
 plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 1)
 plot(1:length(gt[,1]), -log10(gt[,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 1)
 geno.image(cross, chr=i, reorder=1, cex=2)
 plot(c(1,length(X)),c(0,max(Y)),type="n", xlab=paste('chr',i), ylab='physical position')
 points(X,Y)
dev.off()

}
