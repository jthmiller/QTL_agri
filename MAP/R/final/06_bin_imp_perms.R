#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
library('qtl')
##library('parallel')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")

mpath <- '/home/jmiller1/QTL_agri/data'

#fl <- paste0(pop,'_imp.mapped.tsp.csv')
#fl <- file.path(mpath,fl)

################################################################################
load(file.path(mpath,paste0(pop,'_imputed.rsave')))
################################################################################

perm_count <- as.numeric(commandArgs(TRUE)[3])
arraynum <- as.numeric(commandArgs(TRUE)[5])
cores <- as.numeric(commandArgs(TRUE)[4])
batch <- round(nind(cross)/2)

print(commandArgs(TRUE))
print(paste('pop =',pop,', perm = ',perm_count,', cores =', cores,', array =',arraynum))
set.seed(arraynum)
print(paste(cores,'cores'))
################################################################################

sex.phen <- pull.pheno(cross, "sex")
names(cross$geno) <- ifelse(names(cross$geno) == "5","X",names(cross$geno))
attr(cross$geno[["X"]], 'class') <- 'X'
cross$pheno$pheno_norm <- nqrank(cross$pheno$Pheno)
################################################################################

bin.mr.perms.2 <- scantwo(cross, pheno.col=4, model="binary", method="mr",
 clean.output=T, clean.nmar=50, clean.distance=50, maxit=1000,
 assumeCondIndep=T, n.cluster=cores, n.perm=perm_count, perm.Xsp=F,
 verbose=2, batchsize=batch)

bin.mr.perms.pens <- calc.penalties(bin.mr.perms.2, alpha=0.1)

bin.mr.perms.1 <- scanone(cross, pheno.col=4, model='binary', method = "mr",
 n.perm = 10, n.cluster=cores)

lod <- summary(bin.mr.perms.1)[1]

assign(paste0("bin.mr.perms.2.",arraynum), bin.mr.perms.2)
assign(paste0("bin.mr.perms.1.",arraynum), bin.mr.perms.1)

summary(bin.mr.perms.2)
summary(bin.mr.perms.pens)
print(bin.mr.perms.pens)
print(paste('done with', perm_count, 'scan 2 permutations'))
print(summary(bin.mr.perms.1))

################################################################################
save.image(file.path(mpath,paste0(pop,arraynum,'_scan_perms_bin_mr.rsave')))
################################################################################
