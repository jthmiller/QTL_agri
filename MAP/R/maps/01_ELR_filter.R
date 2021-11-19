#!/bin/R

pop <- 'ELR'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales","doParallel")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
mpath <- '/home/jmiller1/QTL_agri/data'
library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)
library(scales)

################################################################################
## read in the QTL cross
cross <- read.cross(file = file.path(mpath, paste0(pop, ".unphased.f2.csvr")),
format = "csvr", geno = c(1:3), estimate.map = FALSE)
################################################################################

################################################################################
### Pull names from plinkfile
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

################################################################################
### The parents are wrong in this population
pars <- c('BLI_BI1124M','ELR_ER1124F')
cross.par <- subset(cross,ind=cross$pheno$ID %in% pars)
################################################################################

################################################################################
### Filter samples that are duplicates of one another
#cg <- comparegeno(cross)
#plot_test('raw_rela_elr')
#hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
#rug(cg[lower.tri(cg)])
#dev.off()
#wh <- which(cg > 0.75, arr=TRUE)
#wh <- wh[wh[,1] < wh[,2],]
#sim <- cg[cbind(wh[,1],wh[,2])]
#matches <- cbind(rownames(wh), rownames(cg)[wh[,2]])
#set_a <- cross$pheno[cross$pheno$ID %in% matches[,1],]
#set_b <- cross$pheno[cross$pheno$ID %in% matches[,2],]
#mis <- nmissing(cross)/sum(nmar(cross))
#set_a <- cbind(set_a, percent_missing = mis[set_a$ID])
#set_b <- cbind(set_b, percent_missing = mis[set_b$ID])
#set_a <- cbind(sim = sim, set_a)
#diff_pheno <- cbind(set_a,set_b)[which(!set_a$bin == set_b$bin),]
#same_pheno <- cbind(set_a,set_b)[which(set_a$bin == set_b$bin),]
#same_pheno[,c(1,4,7,10,13)]
## ELR_10972 = ELR_10982 (different phenotypes)
## ELR_10983 = ELR_10981 (different phenotypes)
## ELR_10869 might be a parent
## MATCHING, SAME PHENO (only throw one of two)
#toss.badata <- c(,"ELR_10987","ELR_11580")

toss.badata <- c('ELR_10869','ELR_10972','ELR_10982','ELR_10983','ELR_10981')
toss.repdata <- c("ELR_10977","ELR_10987",'ELR_10980','ELR_10984','ELR_10974','ELR_10971','ELR_10988')

################################################################################

################################################################################
### Toss individuals that have high missing data
## Many missing ELR_10998 ELR_10924 ELR_10953
toss.missing <- names(which(nmissing(cross)/(sum(nmar(cross))) > 0.25))
################################################################################

################################################################################
## Drop ########################################################################
#cross <- subset(cross, ind=!cross$pheno$ID %in% c(toss.badata,toss.missing,pars))
cross <- subset(cross, ind=!cross$pheno$ID %in% c(toss.repdata,toss.badata,toss.missing,pars))
################################################################################

################################################################################
mapfile <- paste0(pop,'_samples_checked')
filename <- file.path(mpath,mapfile)
#write.cross(cross, filestem=filename, format="csv")
cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

################################################################################
### TOSS MARKERS WITH HIGH PERCENTAGE OF MISSING DATA ##########################
misg <- function(X,perc) { nind(cross) * perc }
mis <- misg(cross,0.125)
drop <- names(which(colSums(is.na(pull.geno(cross))) > mis))
cross <- drop.markers(cross,drop)
################################################################################

sum(nmar(cross))
#[1] 151344

################################################################################
### DROP DISTORTED UNMAPPED (pvalues later shown to retain density of markers on all LGs)
gt <- geno.table(cross)
toss <- rownames(gt[which(gt[,'P.value'] < 1.0e-2),])
cross <- drop.markers(cross,toss)
i <- 5 ; plotit(cross,'filt_0.01')
################################################################################

sum(nmar(cross))
##[1] 58687

################################################################################
mapfile <- paste0(pop,'_checked_filt_0.01')
filename <- file.path(mpath,mapfile)

write.cross(cross, filestem=filename, format="csv")
#cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

## Keep highest lod marker in final crosses to get the markers close to genomic positions
prescan_bin <- scanone(cross, pheno.col=4, method="mr", model="bin", n.cluster = cores)
keeps <- summary(prescan_bin)
keeps <- rownames(keeps[keeps$chr %in% 1:24,])
################################################################################

#################################################################################
#### ALL BUT 17 can be filtered down to 1e-2. Truncates these LGS (see plots)
#gt.sub <- geno.table(cross,chr=c(1:16,18:24))
#toss.sub <- rownames(gt.sub[which(gt.sub[,'P.value'] < 5.0e-2),])
#cross <- drop.markers(cross,toss.sub)
#################################################################################

#cross_chr1 <- est.rf(subset(cross, chr=8))
#rf <- pull.rf(cross_chr1, what = 'lod')
#plot_test('chr_8_elr_lod')
#hist(as.numeric(rf), breaks = 100)
#dev.off()
#
#rf <- pull.rf(cross_chr1)
#lod <- pull.rf(cross_chr1, what = 'lod')
#
#plot_test('chr_8_elr_lodvrf')
#plot(as.numeric(rf),as.numeric(lod))
#dev.off()

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
## [1] 56069
################################################################################

### WRITE THE ABOVE CROSS OBJECT ###############################################
mapfile <- paste0(pop,'_filtered_linked')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")
#cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

################################################################################
## Imputed
cross_imp <- fill.geno(cross, error.prob = 0.05, map.function= 'kosambi')
i <- 9 ; plotit(cross_imp,'elr_cross')
################################################################################

################################################################################
sum(nmar(cross_imp))
# [1] 56069
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
# [1] 1461
################################################################################

####################################################################################
cross_reorg <- formLinkageGroups(cross_rm_exact, max.rf = 0.05, min.lod = 10, reorgMarkers = TRUE)
i <- 100 ; plotit(cross_reorg,'cross_reorg')
drop <- markernames(cross_reorg, chr = 25:nchr(cross_reorg))
cross_reorder <- drop.markers(cross_rm_exact, drop)
####################################################################################

####################################################################################
cross_reorder <- cross_rm_exact
sum(nmar(cross_reorder))
#[1] 1402
cross_reorder  <- est.rf(cross_reorder)
cross_reorder  <- tspOrder(cross = cross_reorder , hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
i <- 100 ; plotit(cross_reorder,'cross_imp_exact')
####################################################################################

################################################################################
## Imputed
cross_reorder_imp <- fill.geno(cross_reorder, error.prob = 0.05, map.function= 'kosambi')
i <- 9 ; plotit(cross_reorder_imp ,'cross_reorder_imp')
################################################################################

################################################################################
mapfile <- paste0(pop,'_reorder_imp')
filename <- file.path(mpath,mapfile)
write.cross(cross_reorder_imp,filestem=filename,format="csv")
#cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
cross <- cross_reorder_imp
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
loglik <- foreach(z = seq(along=err), .inorder = T, .export = c("loglik"), .packages = c("qtl")) %dopar% update.lik(z, cross = cross)
loglik <- unlist(loglik)
lod <- (loglik - max(loglik))/log(10)
erprob <- err[which.max(lod)]
print(paste('error lod =',erprob))
################################################################################

################################################################################
newmap <- est.map(cross, error.prob=erprob, map.function="kosambi", maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)
cross <- replace.map(cross, newmap)
################################################################################

################################################################################
## Flip if the chr is backwards
direc <- sapply(c(1:24),function(i) {
 pos <- as.numeric(gsub(".*:","",markernames(cross,i)))
 map <- as.numeric(pull.map(cross)[[i]])
 cor(pos,map, use="complete.obs")
})

if(any(direc < 0)) cross <- flip.order(cross,which(direc < 0))

i <- 10 ; plotit(cross ,'cross_new_ori')
################################################################################

################################################################################
mapfile <- paste0(pop,'_reorder_imp_mapped')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")
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
## NON IMPUTED CROSS
################################################################################

################################################################################
### READ IN PRE PERMUTED CROSS OB ##############################################
mapfile <- paste0(pop,'_filtered_linked')
filename <- file.path(mpath,mapfile)
noperms <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

################################################################################
mapfile <- paste0(pop,'_reorder_imp_mapped')
filename <- file.path(mpath,mapfile)
cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
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
newmap <- est.map(noperms, error.prob = erprob, map.function="kosambi", maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)
noperms <- replace.map(noperms, newmap)
summary(pull.map(noperms))
################################################################################

################################################################################
mapfile <- paste0(pop,'_reorder_noimp')
filename <- file.path(mpath,mapfile)
write.cross(noperms, filestem = filename, format="csv")
i <- 9 ; plotit(noperms,'cross_reorder_noimp')
################################################################################

#save.image(file.path(mpath,paste0(pop,'_imputed.rsave')))

################################################################################
## done
################################################################################
