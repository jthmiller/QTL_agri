#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
perm_count <- as.numeric(commandArgs(TRUE)[3])
cores <- as.numeric(commandArgs(TRUE)[4])

print(commandArgs(TRUE))
print(paste(pop,perm_count))

library('qtl')
##library('parallel')
library('snow')

source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################

print(paste(cores,'cores'))
erp <- 0.001

load(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))
################################################################################

norm.add.imp.perms <- scanone(cross, pheno.col=5, model='normal', method = "imp", n.perm = 100, n.cluster = 10)
lod <- summary(norm.add.imp.perms)[2]

cross <- sim.geno(cross, n.draws=160, map.function="kosambi", stepwidth="fixed",error.prob=0.001)

norm.add.imp <- scanone(cross, pheno.col=5, model='normal', method = "imp")

qtl <- summary(norm.add.imp,lod)
qtl <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="draws")
qtl

#full.bin.em.step <- stepwiseqtl(cross, model='normal', method = "imp", pheno.col = 5,
 #penalties=pens, incl.markers=T, qtl=qtl, additive.only = F,  scan.pairs = T, max.qtl=8)

full.bin.em.step <- stepwiseqtl(cross, model='normal', method = "imp", pheno.col = 5,
 incl.markers=T, qtl=qtl, additive.only = F,  scan.pairs = T, max.qtl=8)



summary(full.norm.imp.step)
################################################################################
save.image(file.path(mpath,paste0(pop,'_step_norm_imp.rsave')))
################################################################################
