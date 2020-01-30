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
load(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))

perm_count <- as.numeric(commandArgs(TRUE)[3])
arraynum <- as.numeric(commandArgs(TRUE)[5])
cores <- as.numeric(commandArgs(TRUE)[4])
batch <- round(nind(cross)/2)

print(commandArgs(TRUE))
print(paste('pop =',pop,', perm = ',perm_count,', cores =', cores,', array =',arraynum))

################################################################################

print(paste(cores,'cores'))
erp <- 0.0025
sex.phen <- pull.pheno(cross, "sex")
names(cross$geno) <- ifelse(names(cross$geno) == "5","X",names(cross$geno))
attr(cross$geno[["X"]], 'class') <- 'X'

(summary(pull.map(cross))['overall','length']) / (length(colnames(pull.genoprob(cross)))/3)
print('markers per CM')

length(colnames(pull.genoprob(cross)))/3
print('markers')

################################################################################
if(pop == 'ELR'){
 cross <- subset(cross, chr=c(1:4,6:17,19:24))
} else {
 cross <- subset(cross, chr=c(1,3:4,6:24))
}
################################################################################

bin.em.perms.2 <- scantwo(cross, pheno.col=4, model="binary", method="em",
 incl.markers=F, clean.output=T, clean.nmar=25, clean.distance=25, maxit=2000,
 assumeCondIndep=T, n.cluster=cores, addcovar=g, n.perm=perm_count, perm.Xsp=F,
 verbose=2, batchsize=batch)

bin.em.perms.pens <- calc.penalties(bin.em.perms.2, alpha=0.1)

bin.em.perms.1 <- scanone(cross, pheno.col=4, model='binary', method = "em",
 n.perm = 10, n.cluster=cores, addcovar=g)

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
