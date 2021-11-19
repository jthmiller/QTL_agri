#!/bin/R
pop <- 'NEW'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales","doParallel")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
mpath <- '/home/jmiller1/QTL_agri/data'
library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)

## Map a smaller subset to id qtls
## go back and remap chrs with a qtl with a larger sample set
#################################################################################
### read in the QTL cross
cross <- read.cross(file = file.path(mpath, paste0(pop, ".unphased.f2.csvr")),
format = "csvr", geno = c(1:3), estimate.map = FALSE)
#################################################################################

################################################################################
### Pull sample names from plinkfile
path <- file.path(mpath, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross$pheno$ID <- paste(popname, indname, sep = "_")
################################################################################

#### PHENO #####################################################################
cross$pheno$bin <- ifelse(cross$pheno$Pheno > 2, 1 , 0)
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno))
################################################################################

#### SEX #######################################################################
sex <- read.table(file.path(mpath,'sex.txt'),stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross$pheno$ID), 'sex']
cross$pheno$sex <- sex.vec
################################################################################

## DROP PARENTS ################################################################
pars <- c("NEW_NEW911M","NEW_NEW911F")
parc <- subset(cross, ind = cross$pheno$ID %in% pars)
#plot_test('par_nbh_keep_many', width=4000, height = 5000)
#par(mfrow = c(24,1)) ; for(i in 1:24){ geno.image(parc, chr=i)} ; dev.off()
cross <- subset(cross, ind=!cross$pheno$ID %in% c("NEW_NEW911M","NEW_NEW911F"))
################################################################################

################################################################################
### FILTER
################################################################################
## Drop samples that have high missing data or discordant genotype
toss.missing <- names(which(nmissing(cross)/(sum(nmar(cross))) > 0.25))
cross <- subset(cross, ind=!cross$pheno$ID %in% toss.missing)
################################################################################

################################################################################
### TOSS MARKERS WITH HIGH PERCENTAGE OF MISSING DATA ##########################
misg <- function(X, perc) { nind(cross) * perc }
## Drop markers with greater than 12.5% missing data
mis <- misg(cross,0.125)
drop <- names(which(colSums(is.na(pull.geno(cross))) > mis))
print(paste('dropping',length(drop),'markers'))
cross <- drop.markers(cross,drop)
################################################################################

sum(nmar(cross))

################################################################################
### DROP DISTORTED UNMAPPED (pvalues later shown to retain good markers on all LGs)
gt <- geno.table(cross)
toss <- rownames(gt[which(gt[,'P.value'] < 1.0e-2),])
cross <- drop.markers(cross,toss)
i <- 5 ; plotit(cross,'cons')
################################################################################

sum(nmar(cross))

################################################################################
mapfile <- paste0(pop,'_unmapped_filtered')
fl <- file.path(mpath,mapfile)
write.cross(cross,filestem=fl,format="csv")
cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

################################################################################
linked_marks <- function(cross, X, LOD = 10, RF = 1){
 crossX <- est.rf(subset(cross, chr=X))
 crossX <- formLinkageGroups(crossX, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 markernames(crossX, chr=1)
}

linked <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% linked_marks(cross = cross, X)
cross <- pull.markers(cross, unlist(linked))
i <- 95 ; plotit(cross,'cross_linked')
################################################################################

################################################################################
mapfile <- paste0(pop,'_unmapped_filtered_linked')
fl <- file.path(mpath,mapfile)
write.cross(cross,filestem=fl,format="csv")
cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

################################################################################
linked_marks <- function(cross, X, LOD = 15, RF = 0.075){
 crossX <- est.rf(subset(cross, chr=X))
 crossX <- formLinkageGroups(crossX, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 markernames(crossX, chr=1)
}

linked <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% linked_marks(cross = cross, X)
cross <- switchAlleles(cross, unlist(linked))
i <- 95 ; plotit(cross,'cross_linked_new')
linked <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% linked_marks(cross = cross, X)
cross <- switchAlleles(cross, unlist(linked))
i <- 96 ; plotit(cross,'cross_linked_new')
################################################################################

cross.sub <- subset(cross,chr=1)
cross.sub <- est.rf(cross.sub, maxit=1000, tol=1e-6)
rf <- pull.rf(cross.sub, what='lod')
freq <- rowSums(rf, na.rm=T)/dim(rf)[1]
plot_test('chr1_new')
hist(as.numeric(freq), breaks = 50)
dev.off()

################################################################################
switch.phase <- function(i, cross){
 cross.sub <- subset(cross,chr=i)
 cross.sub <- est.rf(cross.sub, maxit=1000, tol=1e-6)
 rf <- pull.rf(cross.sub)
 freq <- rowSums(rf, na.rm=T)/dim(rf)[1]
 return(names(freq)[which(freq > 0.40)])
}
################################################################################

################################################################################
switch.many <- foreach(i = 1:24, .inorder = F, .packages = libs2load) %dopar% switch.phase(i, cross = cross)
cross <- switchAlleles(cross, as.character(unlist(switch.many)))
i <- 95 ; plotit(cross,'cross_linked_sw')
################################################################################

################################################################################
linked_marks <- function(cross, X, LOD = 8, RF = 0.1){
 crossX <- est.rf(subset(cross, chr=X))
 crossX <- formLinkageGroups(crossX, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 markernames(crossX, chr=1)
}

linked <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% linked_marks(cross = cross, X)
cross <- pull.markers(cross, unlist(linked))
i <- 95 ; plotit(cross,'cross_linked')
################################################################################

cross.imp <- fill.geno(cross, error.prob = 0.01, map.function= 'kosambi')
i <- 9 ; plotit(cross.imp,'cross.cons')

cross.rm <- findDupMarkers(cross.imp, exact.only=FALSE, adjacent.only=TRUE) # finds 6 pairs
cross.imp <- drop.markers(cross.imp, unlist(cross.rm))
i <- 10 ; plotit(cross.imp,'cross.many')

################################################################################
mapfile <- paste0(pop,'_cons_filt_refined_switch')
filename <- file.path(mpath,mapfile)
write.cross(cross.imp, filestem=filename, format="csv")
################################################################################

cross.imp  <- est.rf(cross.imp)
cross.imp  <- tspOrder(cross = cross.imp , hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
i <- 100 ; plotit(cross.imp,'cross.imp')

################################################################################
mapfile <- paste0(pop,'_cons_filt_refined_switch_mapped')
filename <- file.path(mpath,mapfile)
write.cross(cross.imp, filestem=filename, format="csv")
################################################################################
