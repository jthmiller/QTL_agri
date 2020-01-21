#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
perm_count <- as.numeric(commandArgs(TRUE)[3])
cores <- as.numeric(commandArgs(TRUE)[4])

print(commandArgs(TRUE))
print(paste(pop,perm_count))

library('qtl')
##library('parallel')
library('snow')

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################

print(paste(cores,'cores'))
erp <- 0.0025

################################################################################

################################################################################
load(file.path(mpath,paste0(pop,'_downsampled.rsave')))
################################################################################

sex.phen <- pull.pheno(cross, "sex")
names(cross$geno) <- ifelse(names(cross$geno) == "5","X",names(cross$geno))

bin.em.2 <- scantwo(cross, pheno.col=4, model="binary", method="em",
 incl.markers=F,clean.output=T, clean.nmar=15, clean.distance=15,
 assumeCondIndep=T, n.cluster=cores, intcovar=sex.phen, maxit=2000)

################################################################################
save.image(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
################################################################################
