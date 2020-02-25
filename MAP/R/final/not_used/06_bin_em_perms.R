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
load(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))
################################################################################

perm_count <- as.numeric(commandArgs(TRUE)[3])
cores <- as.numeric(commandArgs(TRUE)[4])
arraynum <- as.numeric(commandArgs(TRUE)[5])
method <- as.numeric(commandArgs(TRUE)[6])
model <- as.numeric(commandArgs(TRUE)[7])
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

################################################################################
if(pop == 'NBH'){
 mar <- '2:27373969'
 g <- pull.geno(fill.geno(cross))[,mar]
 g <- cbind(as.numeric(g==1), as.numeric(g==2))
} else {
 mar <- '18:20422142'
 g <- pull.geno(fill.geno(cross))[,mar]
 g <- cbind(as.numeric(g==1), as.numeric(g==2))
}
################################################################################

perms.2 <- scantwo(cross, pheno.col=4, model=model, method = method,
 clean.output=T, clean.nmar=50, clean.distance=50, maxit=100,
 assumeCondIndep=T, n.cluster=cores, n.perm=perm_count, perm.Xsp=F,
 verbose=2, batchsize=batch, addcovar=g)

pens <- calc.penalties(perms.2, alpha=0.1)

perms.1 <- scanone(cross, pheno.col=4, model=model, method = method,
 n.perm = 10, n.cluster=cores)

lod <- summary(perms.1)[1]

assign(paste0("perms.2.",arraynum), perms.2)
assign(paste0("perms.1.",arraynum), perms.1)

summary(perms.2)
summary(pens)
print(pens)
print(paste('done with', perm_count, 'scan 2 permutations'))
print(summary(perms.1))

################################################################################
save.image(file.path(mpath,paste0(pop,arraynum,'_scan_perms_',model,'_',method,'.rsave')))
################################################################################
