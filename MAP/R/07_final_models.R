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
sw_elr <- full.norm.add_only
pni_elr <- perms.norm.imp
pbe_elr <- perms.bin.em
sbe_elr <- scan.bin.em
sbm_elr <- scan.bin.mr

load(file.path(mpath,'single_scans.nbh.rsave'))
sw_nbh <- full.norm.add_only
pni_nbh <- perms.norm.imp
pbe_nbh <- perms.bin.em
sbe_nbh <- scan.bin.em
sbm_nbh <- scan.bin.mr

load(file.path(mpath,'single_scans.new.rsave'))
sw_new <- full.norm.add_only
pni_new <- perms.norm.imp
pbe_new <- perms.bin.em
sbe_new <- scan.bin.em
sbm_new <- scan.bin.mr

load(file.path(mpath,'single_scans.brp.rsave'))
sw_brp <- full.norm.add_only
pni_brp <- perms.norm.imp
pbe_brp <- perms.bin.em
sbe_brp <- scan.bin.em
sbm_brp <- scan.bin.mr

#short
load(file.path(mpath,'scantwo.scans.elr.short.rsave'))


png(paste0('~/public_html/ELR_model.png'))
 plotModel(full.norm.add_only)
dev.off()


full.norm.add_only
