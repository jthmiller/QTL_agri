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
load(file.path(mpath,paste0(pop,'_scan_perms_norm_imp.rsave')))
################################################################################

norm.add.imp <- scanone(cross, pheno.col=5, model='normal', method = "imp")
qtl <- summary(norm.add.imp,lod)
qtl <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="draws")
qtl


summary(full.norm.imp.step)
################################################################################
save.image(file.path(mpath,paste0(pop,'_step_norm_imp.rsave')))
################################################################################
