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

norm.hk.perms.2 <- scantwo(cross, pheno.col=4, model="binary", method="hk",
 incl.markers=F, clean.output=T, clean.nmar=10, clean.distance=10,
 n.perm=perm_count, assumeCondIndep=T, n.cluster=cores)

norm.hk.perms.pens <- calc.penalties(norm.hk.perms.2, alpha=0.1)

norm.hk.perms.1 <- scanone(cross, pheno.col=5, model='normal', method = "hk",
 n.perm = 10000, n.cluster=cores)

lod <- summary(norm.hk.perms.1)[1]

summary(norm.hk.perms.2)
summary(norm.hk.perms.pens)
print(norm.hk.perms.pens)
print(paste('done with', perm_count, 'scan 2 permutations'))
print(summary(norm.hk.perms.1))

################################################################################
save.image(file.path(mpath,paste0(pop,'_scan_perms_bin_hk.rsave')))
################################################################################
