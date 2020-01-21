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
load(file.path(mpath,paste0(pop,'_scan_perms_bin_hk.rsave')))
################################################################################

bin.add.hk <- scanone(cross, pheno.col=5, model='binary', method = "hk")
qtl <- summary(bin.add.hk,lod)
qtl <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="draws")
qtl

full.bin.hk.step <- stepwiseqtl(cross, model='binary', method = "hk", pheno.col = 4,
 penalties=pens, incl.markers=F, qtl=qtl, additive.only = F,  scan.pairs = T, max.qtl=8)

summary(full.bin.hk.step)
################################################################################
save.image(file.path(mpath,paste0(pop,'_step_bin_hk.rsave')))
################################################################################
