#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
library('qtl')
library('snow')

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

load(file.path(mpath,paste0(pop,1,'_scan_perms_bin_em.rsave')))
perms.2 <- get(paste0('bin.em.perms.2.',1))
#perms.1 <- get(paste0('bin.em.perms.1.',1))

for (i in 2:94){
 arraynum <- i
 load(file.path(mpath,paste0(pop,arraynum,'_scan_perms_bin_em.rsave')))
 perms.2 <- c(perms.2,get(paste0('bin.em.perms.2.',i)))
 #perms_1 <- c(perms,get(paste0('bin.em.perms.2.',i)))
}

pens <- calc.penalties(perms.2, alpha=0.1)

save.image(file.path(mpath,paste0(pop,'_all_perms_bin_em.rsave')))

################################################################################
