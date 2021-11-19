#!/bin/R

### Split dataset: A: mapping cross B: Genotype validation for exploring causative variation
pop <- 'NBH'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales","doParallel")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
mpath <- '/home/jmiller1/QTL_agri/data'
library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)

## Map a smaller subset to id qtls
## go back and remap chrs with a qtl with a larger sample set

################################################################################
## read in the QTL cross
umpath <- '/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops'
fl <- 'NBH.um.unmapped.f2.csvr'
cross <- read.cross(file = file.path(umpath, fl),
format = "csvr", geno = c(1:3), estimate.map = FALSE)
################################################################################

#################################################################################
### read in the QTL cross
#cross <- read.cross(file = file.path(mpath, paste0(pop, ".unphased.f2.csvr")),
#format = "csvr", geno = c(1:3), estimate.map = FALSE)
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

### RAW DATA SET ###############################################################
mapfile <- paste0(pop,'_raw')
filename <- file.path(mpath,mapfile)
write.cross(cross, filestem=filename, format="csv")
#cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

################################################################################
i <- 1 ; plotit(cross,'_raw_nopar')
################################################################################

################################################################################
#cg <- comparegeno(cross)
#plot_test('raw_rela_nbh')
#hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
#rug(cg[lower.tri(cg)])
#dev.off()
################################################################################

################################################################################
### FILTER
################################################################################
## Drop samples that have high missing data or discordant genotype
toss.missing <- names(which(nmissing(cross)/(sum(nmar(cross))) > 0.25))
### is "NBH_5646" another grandparent sample??
### after filtering, NBH_6137 appears to have high allelic dropout
toss.missing <- c(toss.missing,"NBH_5646")
cross <- subset(cross, ind=!cross$pheno$ID %in% toss.missing)
################################################################################

################################################################################
## Parent markers
## AA, AB, BB are displayed in the colors red, blue, and green,

## DROP PARENTS ################################################################
pars <- c('NBH_NBH1M','NBH_NBH1F')
parc <- subset(cross, ind = cross$pheno$ID %in% pars)
#plot_test('par_nbh_keep_many', width=4000, height = 5000)
#par(mfrow = c(24,1)) ; for(i in 1:24){ geno.image(parc, chr=i)} ; dev.off()
cross <- subset(cross, ind=!cross$pheno$ID %in% pars)
################################################################################

################################################################################
mapfile <- paste0(pop,'_nopar_nofilt')
filename <- file.path(mpath,mapfile)
write.cross(cross, filestem=filename, format="csv")
# cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
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

################################################################################
mapfile <- paste0(pop,'_nopar_nofilt_refined')
filename <- file.path(mpath,mapfile)

# write.cross(cross, filestem=filename, format="csv")
# cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath,
#  genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

################################################################################
### DROP DISTORTED UNMAPPED (pvalues later shown to retain density of markers on all LGs)
gt <- geno.table(cross)
toss <- rownames(gt[which(gt[,'P.value'] < 1.0e-3),])
cross <- drop.markers(cross,toss)
i <- 5 ; plotit(cross,'filt_0.001')
################################################################################

################################################################################
mapfile <- paste0(pop,'_nopar_filt')
filename <- file.path(mpath,mapfile)
#write.cross(cross, filestem=filename, format="csv")
cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

################################################################################
sum(nmar(cross))
# [1] 42930
################################################################################

## Keep highest lod marker in final crosses to get the markers close to genomic positions
prescan_bin <- scanone(cross, pheno.col=4, method="mr", model="bin", n.cluster = cores)
keeps <- summary(prescan_bin)
keeps <- rownames(keeps[keeps$chr %in% 1:24,])
################################################################################

################################################################################
## Drop unlinked markers
linked_marks <- function(cross, X, LOD = 10, RF = 0.1){
 crossX <- est.rf(subset(cross, chr=X))
 crossX <- formLinkageGroups(crossX, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 a <- markernames(crossX, chr=1:2)
 b <- markernames(crossX, chr=2)
 return(list(keep=a,switch=b))
}

linked <- foreach(X = 1:24, .inorder = F, .packages = libs2load) %dopar% linked_marks(cross = cross, X)
pull <- lapply(linked,"[[",1)
switch <- lapply(linked,"[[",2)
cross <- switchAlleles(cross, unlist(switch))
cross <- pull.markers(cross, c(unlist(pull), keeps))
i <- 95 ; plotit(cross,'linked_10lod')
################################################################################

################################################################################
sum(nmar(cross))
# [1] 35665
################################################################################

### WRITE THE ABOVE CROSS OBJECT ###############################################
mapfile <- paste0(pop,'_filtered_linked_nopar')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")
#cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

################################################################################
## Imputed
cross_imp <- fill.geno(cross, error.prob = 0.05, map.function= 'kosambi')
i <- 9 ; plotit(cross_imp,'nopar_imp')
################################################################################

################################################################################
cross_rm_exact <- findDupMarkers(cross_imp, exact.only=TRUE, adjacent.only=TRUE) # finds 6 pairs
cross_rm_exact <- unlist(cross_rm_exact)
cross_rm_exact <- cross_rm_exact[!cross_rm_exact %in% keeps]
cross_rm_exact <- drop.markers(cross_imp, cross_rm_exact)
i <- 9 ; plotit(cross_rm_exact,'cross_exact_nopar')
################################################################################

################################################################################
sum(nmar(cross_rm_exact))
# [1]  1774
################################################################################

####################################################################################
#cross_reorg <- formLinkageGroups(cross_rm_exact, max.rf = 0.05, min.lod = 10, reorgMarkers = TRUE)
#i <- 100 ; plotit(cross_reorg,'cross_reorg_nopar')
#drop <- markernames(cross_reorg, chr = 24:nchr(cross_reorg))
# cross_reorder <- drop.markers(cross_rm_exact, drop)
####################################################################################

####################################################################################
sum(nmar(cross_rm_exact))
# [1] 1757
cross_reorder <- cross_rm_exact
cross_reorder  <- est.rf(cross_reorder)
cross_reorder  <- tspOrder(cross = cross_reorder , hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
i <- 100 ; plotit(cross_reorder,'cross_imp_exact_nopar')
####################################################################################

################################################################################
## Imputed
cross_reorder_imp <- fill.geno(cross_reorder, error.prob = 0.05, map.function= 'kosambi')
i <- 9 ; plotit(cross_reorder_imp ,'cross_reorder_imp_nopar')
################################################################################

################################################################################
## Flip if the chr is backwards
direc <- sapply(c(1:24),function(i) {
 pos <- as.numeric(gsub(".*:","",markernames(cross_reorder_imp,i)))
 map <- as.numeric(pull.map(cross_reorder_imp)[[i]])
 cor(pos,map, use="complete.obs")
})

if(any(direc < 0)) cross_reorder_imp <- flip.order(cross_reorder_imp,which(direc < 0))

i <- 100 ; plotit(cross_reorder_imp,'cross_reori')
################################################################################

################################################################################
mapfile <- paste0(pop,'_reorder_imp_nopar')
filename <- file.path(mpath,mapfile)
write.cross(cross_reorder_imp,filestem=filename,format="csv")
#cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

################################################################################
png(paste0('~/public_html/',pop,'all_rf.png'),height=2500,width=2500)
par(mfrow=c(6,4))
for (B in 1:24){
 plotRF(cross_reorder_imp, chr = B)
}
dev.off()
################################################################################

################################################################################
loglik <- err <- c(0.0001, 0.001, 0.01, 0.05)

update.lik <- function(z, cross){
  cat(z, "of", length(err), "\n")
  tempmap <- est.map(cross,maxit=100, error.prob=err[z])
  loglik[z] <- sum(sapply(tempmap, attr, "loglik"))
}

cl <- makeCluster(5)
registerDoParallel(cl)
loglik <- foreach(z = seq(along=err), .inorder = T, .export = c("loglik"), .packages = c("qtl")) %dopar% update.lik(z, cross = cross_reorder_imp)
loglik <- unlist(loglik)
lod <- (loglik - max(loglik))/log(10)
erprob <- err[which.max(lod)]
print(paste('error lod =',erprob))
################################################################################

################################################################################
cross <- cross_reorder_imp
newmap <- est.map(cross, error.prob=erprob, map.function="kosambi", maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)
cross <- replace.map(cross, newmap)
summary(pull.map(cross))
################################################################################

################################################################################
## FINALMAP
mapfile <- paste0(pop,'_reorder_imp_nopar')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")
#cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

################################################################################
################################################################################

### READ IN PRE PERMUTED CROSS OB ###############################################
mapfile <- paste0(pop,'_filtered_linked_nopar')
filename <- file.path(mpath,mapfile)
noperms <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

################################################################################
noperms <- pull.markers(noperms,markernames(cross))
noperms  <- est.rf(noperms)
noperms  <- tspOrder(cross = noperms , hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
i <- 100 ; plotit(noperms, 'cross_noimp_exact_nopar')
################################################################################

################################################################################
loglik <- err <- c(0.0001, 0.001, 0.01, 0.05)

update.lik <- function(z, cross){
  cat(z, "of", length(err), "\n")
  tempmap <- est.map(cross,maxit=100, error.prob=err[z])
  loglik[z] <- sum(sapply(tempmap, attr, "loglik"))
}

cl <- makeCluster(5)
registerDoParallel(cl)
loglik <- foreach(z = seq(along=err), .inorder = T, .export = c("loglik"), .packages = c("qtl")) %dopar% update.lik(z, cross = noperms)
loglik <- unlist(loglik)
lod <- (loglik - max(loglik))/log(10)
erprob <- err[which.max(lod)]
print(paste('error lod =',erprob))
################################################################################

################################################################################
newmap <- est.map(noperms, error.prob=erprob, map.function="kosambi", maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)
noperms <- replace.map(noperms, newmap)
summary(pull.map(noperms))
################################################################################

################################################################################
mapfile <- paste0(pop,'_reorder_noimp_nopar')
filename <- file.path(mpath,mapfile)
write.cross(cross_reorder_imp,filestem=filename,format="csv")
#cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

#save.image(file.path(mpath,paste0(pop,'_imputed.rsave')))


################################################################################
## done
################################################################################
