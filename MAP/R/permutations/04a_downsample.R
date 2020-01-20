#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
library('qtl')
##library('parallel')
library('snow')

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)
##cores <- detectCores() - 2

################################################################################
################################################################################
erp <- 0.0025
cores <- 22
################################################################################
################################################################################
## Read cross
cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno),5)
cross <- jittermap(cross)
dups <- findDupMarkers(cross, exact.only=F, adjacent.only=F)
dups <- names(dups)
if(pop == 'ELR.missing') dups <- c(dups,"AHR2a_del")
cross <- pull.markers(cross, dups)
cross
################################################################################
fl <- file.path(mpath,paste0(pop,'_downsampled'))
write.cross(cross,filestem=fl,format="csv")
################################################################################

cross$pheno <- as.data.frame(cross$pheno)
cross <- sim.geno(cross, stepwidth="fixed", step=1,off.end=5, error.prob=erp ,map.function="kosambi", n.draws=100)
cross <- calc.genoprob(cross, stepwidth="fixed", step=1, error.prob=erp, off.end=5)
cross <- reduce2grid(cross)

################################################################################
save.image(file.path(mpath,paste0(pop,'_downsampled.rsave')))
################################################################################
