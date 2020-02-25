#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
library('qtl')
library('snow')

permname <- '_all_perms_bin_em.rsave'


source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

load(file.path(mpath,paste0(pop,1,permname)))
perms.2 <- get(paste0('bin.em.perms.2.',1))
perms.1 <- get(paste0('bin.em.perms.1.',1))

for (i in 2:200){
 arraynum <- i
 load(file.path(mpath,paste0(pop,arraynum,permname)))

 print(paste('done with array',i))

 nm <- paste0('bin.em.perms.2.',i)
 perms.2 <- c(perms.2,get(nm))
 rm(nm)

 nm <- paste0('bin.em.perms.1.',i)
 perms.1 <- c(perms.1,get(nm))
 rm(nm)
}

summary(perms.2)

save.image(file.path(mpath,paste0(pop,permname)))

################################################################################
