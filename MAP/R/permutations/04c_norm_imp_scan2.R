#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
cores <- as.numeric(commandArgs(TRUE)[3])

print(commandArgs(TRUE))
print(paste(pop))

library('qtl')
##library('parallel')
library('snow')

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################
load(file.path(mpath,paste0(pop,'_downsampled.rsave')))
################################################################################

################################################################################

print(paste(cores,'cores'))
erp <- 0.0025
sex.phen <- pull.pheno(cross, "sex")
names(cross$geno) <- ifelse(names(cross$geno) == "5","X",names(cross$geno))

################################################################################
cov <- ifelse(pop == 'ELR',18,2)
so <- summary(scanone(cross,pheno.col=5, model="normal", method="imp", intcovar=sex.phen))[cov,]
mar <- find.marker(cross, so$chr, so$lod)
g <- pull.geno(fill.geno(cross))[,mar]
g <- cbind(as.numeric(g==1), as.numeric(g==2))
summary(scanone(cross,pheno.col=5, model="normal", method="imp",addcovar=g))
################################################################################

norm.imp.2 <- scantwo(cross, pheno.col=5, model="normal", method="imp",
 incl.markers=F, chr = c(1:4,6:24),clean.output=T, clean.nmar=200, clean.distance=200,
 assumeCondIndep=T, n.cluster=cores, intcovar=sex.phen, addcovar=g)

################################################################################
save.image(file.path(mpath,paste0(pop,'_scan2_norm_imp.rsave')))
################################################################################
