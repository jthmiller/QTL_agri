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
################################################################################
print(paste(cores,'cores'))
erp <- 0.0025
sex.phen <- pull.pheno(cross, "sex")
names(cross$geno) <- ifelse(names(cross$geno) == "5","X",names(cross$geno))

################################################################################
cov <- ifelse(pop == 'ELR',18,2)
so <- summary(scanone(cross,pheno.col=4, model="binary", method="hk", intcovar=sex.phen))[cov,]
mar <- find.marker(cross, so$chr, so$lod)
g <- pull.geno(fill.geno(cross))[,mar]
g <- cbind(as.numeric(g==1), as.numeric(g==2))
summary(scanone(cross,pheno.col=4, model="binary", method="hk",addcovar=g))
################################################################################

bin.hk.perms.2 <- scantwo(cross, pheno.col=4, model="binary", method="hk",
 incl.markers=F, clean.output=T, clean.nmar=200, clean.distance=200, maxit=1000,
 assumeCondIndep=T, n.cluster=cores, intcovar=sex.phen, addcovar=g, n.perm=perm_count)

bin.hk.perms.pens <- calc.penalties(bin.hk.perms.2, alpha=0.1)

bin.hk.perms.1 <- scanone(cross, pheno.col=4, model='binary', method = "hk",
 n.perm = 10000, n.cluster=cores, intcovar=sex.phen, addcovar=g)

lod <- summary(bin.hk.perms.1)[1]

summary(bin.hk.perms.2)
summary(bin.hk.perms.pens)
print(bin.hk.perms.pens)
print(paste('done with', perm_count, 'scan 2 permutations'))
print(summary(bin.hk.perms.1))
################################################################################
save.image(file.path(mpath,paste0(pop,'_scan_perms_bin_hk.rsave')))
################################################################################
