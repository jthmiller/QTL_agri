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
##load(file.path(mpath,paste0(pop,'_downsampled.rsave')))
################################################################################

################################################################################

print(paste(cores,'cores'))
erp <- 0.0025
sex.phen <- pull.pheno(cross, "sex")
names(cross$geno) <- ifelse(names(cross$geno) == "5","X",names(cross$geno))
attr(cross$geno[["X"]], 'class') <- 'X'

################################################################################

sone <- scanone(cross,pheno.col=4, model="binary", method="mr")
sone.perms <- scanone(cross,pheno.col=4, model="binary", method="mr", perms=10000)
summary(sone.perms)


################################################################################
save.image(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
################################################################################

################################################################################

bin.em.2 <- scantwo(cross, pheno.col=4, model="binary", method="mr",
 clean.output=T, clean.nmar=25, clean.distance=25, maxit=1000,
 assumeCondIndep=T, n.cluster=1)

################################################################################
save.image(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
################################################################################



(summary(pull.map(cross))['overall','length']) / (length(colnames(pull.genoprob(cross)))/3)
print('markers per CM')

length(colnames(pull.genoprob(cross)))/3
print('markers')



cov <- rownames(summary(sone, perms=sone.perms, alpha=0.1))
so <- summary(sone)[cov,]
top_2 <- order(so$lod,decreasing =T)[1]
mar <- find.pseudomarker(cross, so$chr[top_2], so$pos[top_2])

g <- lapply(mar,function(X){ pull.argmaxgeno(cross)[,X] } )
names(g) <- mar
g <- lapply(g, function(X,Y){ cbind(as.numeric(X==1), as.numeric(X==2))} )
g <- data.frame(do.call(cbind,g))

sone <- scanone(cross,pheno.col=4, model="binary", method="em", addcovar=g)
summary(sone)
