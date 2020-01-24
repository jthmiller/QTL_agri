#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]

library('qtl')
##library('parallel')
library('snow')

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################

################################################################################
load(file.path(mpath,paste0(pop,'_downsampled.rsave')))
################################################################################

perm_count <- as.numeric(commandArgs(TRUE)[3])
cores <- as.numeric(commandArgs(TRUE)[4])
arraynum <- as.numeric(commandArgs(TRUE)[5])

print(commandArgs(TRUE))
print(paste('pop =',pop,', perm = ',perm_count,', cores =', cores))

erp <- 0.0025
sex.phen <- pull.pheno(cross, "sex")
names(cross$geno) <- ifelse(names(cross$geno) == "5","X",names(cross$geno))

################################################################################
cov <- ifelse(pop == 'ELR',18,2)
so <- summary(scanone(cross,pheno.col=4, model="binary", method="em", intcovar=sex.phen))[cov,]
mar <- find.marker(cross, so$chr, so$lod)
g <- pull.geno(fill.geno(cross))[,mar]
g <- cbind(as.numeric(g==1), as.numeric(g==2))
summary(scanone(cross,pheno.col=4, model="binary", method="em",addcovar=g))
################################################################################

bin.em.perms.2 <- scantwo(cross, pheno.col=4, model="binary", method="em",
 incl.markers=F, clean.output=T, clean.nmar=200, clean.distance=200, maxit=1000,
 assumeCondIndep=T, n.cluster=cores, intcovar=sex.phen, addcovar=g, n.perm=perm_count
 verbose=2)

bin.em.perms.pens <- calc.penalties(bin.em.perms.2, alpha=0.1)

bin.em.perms.1 <- scanone(cross, pheno.col=4, model='binary', method = "em",
 n.perm = 10000, n.cluster=cores, intcovar=sex.phen, addcovar=g)

lod <- summary(bin.em.perms.1)[1]

assign(paste0("bin.em.perms.2.",arraynum), bin.em.perms.2)

summary(bin.em.perms.2)
summary(bin.em.perms.pens)
print(bin.em.perms.pens)
print(paste('done with', perm_count, 'scan 2 permutations'))
print(summary(bin.em.perms.1))
################################################################################
save.image(file.path(mpath,paste0(pop,arraynum,'_scan_perms_bin_em.rsave')))
################################################################################
