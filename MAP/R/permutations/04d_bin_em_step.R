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
load(file.path(mpath,paste0(pop,'_all_perms_bin_em.rsave')))
## perms = scan2 permutations
## perms_1 = scanone permutations
## bin.em.perms.pens = 0.1 penalties
load(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
## sone = scanone
## bin.em.2 = scantow
## addcovar = g
################################################################################
##bin.add.em <- scanone(cross, pheno.col=5, model='binary', method = "em")
#qtl <- summary(bin.add.em,lod)
#qtl <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="draws")
#qtl

qtl <- makeqtl(cross, chr=so[top_2,'chr'], pos=so[top_2,'pos'], what="prob")
qtl

cross <- subset(cross,chr=c(1:4,6:24))

full.bin.em.step <- stepwiseqtl(cross, model='binary', method = "hk", pheno.col = 4,
 penalties=pens, incl.markers=F, qtl=qtl, additive.only = F,  scan.pairs = T, max.qtl=8)

summary(full.bin.em.step)
################################################################################
save.image(file.path(mpath,paste0(pop,'_step_bin_em.rsave')))
################################################################################

bin.em.perms.2 <- scantwo(cross, pheno.col=4, model="binary", method="em",
 incl.markers=F, clean.output=T, clean.nmar=25, clean.distance=25, maxit=2000,
 assumeCondIndep=T, n.cluster=cores, addcovar=g, n.perm=perm_count, perm.Xsp=F,
 verbose=2, batchsize=batch)

bin.em.perms.pens <- calc.penalties(bin.em.perms.2, alpha=0.1)

bin.em.perms.1 <- scanone(cross, pheno.col=4, model='binary', method = "em",
 n.perm = 10, n.cluster=cores, addcovar=g)
