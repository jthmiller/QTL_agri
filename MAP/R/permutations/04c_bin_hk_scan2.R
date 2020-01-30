#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]

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
cores <- as.numeric(commandArgs(TRUE)[3])
print(paste(cores,'cores'))
################################################################################

erp <- 0.0025
sex.phen <- pull.pheno(cross, "sex")
names(cross$geno) <- ifelse(names(cross$geno) == "5","X",names(cross$geno))
attr(cross$geno[["X"]], 'class') <- 'X'

################################################################################

(summary(pull.map(cross))['overall','length']) / (length(colnames(pull.genoprob(cross)))/3)
print('markers per CM')

length(colnames(pull.genoprob(cross)))/3
print('markers')

sone <- scanone(cross,pheno.col=4, model="binary", method="hk")

if(pop == 'ELR'){
 sone.perms <- scanone(subset(cross, chr=c(1:4,'X',6:17,19:24)), pheno.col=4, model="binary", method="hk", n.perm=1000, n.cluster = cores, perm.Xsp=T)
} else {
 sone.perms <- scanone(subset(cross, chr=c(1,3,4,'X',6:24)), pheno.col=4, model="binary", method="hk", n.perm=1000, n.cluster = cores, perm.Xsp=T)
}

summary(sone.perms)

cov <- rownames(summary(sone, perms=sone.perms, alpha=0.1))
so <- summary(sone)[cov,]
top_2 <- order(so$lod,decreasing =T)[1]
mar <- find.pseudomarker(cross, so$chr[top_2], so$pos[top_2])

g <- lapply(mar,function(X){ pull.argmaxgeno(cross)[,X] } )
names(g) <- mar
g <- lapply(g, function(X,Y){ cbind(as.numeric(X==1), as.numeric(X==2))} )
g <- data.frame(do.call(cbind,g))

sone <- scanone(cross,pheno.col=4, model="binary", method="hk", addcovar=g)
summary(sone)

################################################################################
save.image(file.path(mpath,paste0(pop,'_scan2_bin_hk.rsave')))
################################################################################

################################################################################
batch <- round(nind(cross)/2)

bin.hk.2 <- scantwo(cross, pheno.col=4, model="binary", method="hk",
 incl.markers=T, clean.output=T, clean.nmar=25, clean.distance=25, maxit=1000,
 assumeCondIndep=T, n.cluster=cores, addcovar=g, batchsize=batch)

################################################################################
save.image(file.path(mpath,paste0(pop,'_scan2_bin_hk.rsave')))
##
