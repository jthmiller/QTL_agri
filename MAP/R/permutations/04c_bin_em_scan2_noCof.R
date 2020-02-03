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
load(file.path(mpath,paste0(pop,'_downsampled.rsave')))
################################################################################

################################################################################

print(paste(cores,'cores'))
erp <- 0.0025

cross <- subset(cross, chr = c(1:4,6:24))

################################################################################

################################################################################

bin.em.2 <- scantwo(cross, pheno.col=4, model="binary", method="em",
 incl.markers=T, clean.output=T, clean.nmar=100, clean.distance=100, maxit=2000,
 assumeCondIndep=T, n.cluster=cores)

bin.em.nofac <- bin.em.2

################################################################################
save.image(file.path(mpath,paste0(pop,'_scan2_bin_em_noCof.rsave')))
################################################################################
