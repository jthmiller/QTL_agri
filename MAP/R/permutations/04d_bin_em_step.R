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
load(file.path(mpath,paste0(pop,'_scan_perms_bin_em.rsave')))
################################################################################

bin.add.em <- scanone(cross, pheno.col=5, model='binary', method = "em")
qtl <- summary(bin.add.em,lod)
qtl <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="draws")
qtl

full.bin.em.step <- stepwiseqtl(cross, model='binary', method = "em", pheno.col = 4,
 penalties=pens, incl.markers=F, qtl=qtl, additive.only = F,  scan.pairs = T, max.qtl=8)

summary(full.bin.em.step)
################################################################################
save.image(file.path(mpath,paste0(pop,'_step_bin_em.rsave')))
################################################################################
