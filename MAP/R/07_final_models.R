#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]

pop <- 'ELR'

library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)


# single scans
load(file.path(mpath,'single_scans.elr.rsave'))
#short
load(file.path(mpath,'scantwo.scans.elr.short.rsave'))


png(paste0('~/public_html/ELR_model.png'))
 plotModel(full.norm.add_only)
dev.off()


full.norm.add_only
