#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]
library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.missing_mapped.tsp.csv')
fl <- file.path(mpath,fl)

load(file.path(mpath,'single_scans.elr_missing.rsave'))

bins <- data.frame(
 em=summary(scan.bin.em),
 mr=summary(scan.bin.mr)[,'lod'],
 np=summary(scan.np.em.b)[,'lod'])

binpo <- data.frame(
 empo=rownames(summary(scan.bin.em)),
 mrpo=rownames(summary(scan.bin.mr)),
 nppo=rownames(summary(scan.np.em.b)))

norms <- data.frame(
 em=summary(scan.norm.em),
 imp=summary(scan.norm.imp)[,'lod'],
 mr=summary(scan.norm.mr)[,'lod'],
 np=summary(scan.np.em.n)[,'lod'],
 ehk=summary(scan.norm.ehk)[,'lod'])

normpo <- data.frame(
 empo=rownames(summary(scan.norm.em)),
 impo=rownames(summary(scan.norm.imp)),
 mrpo=rownames(summary(scan.norm.mr)),
 nppo=rownames(summary(scan.np.em.n)))

png(paste0('~/public_html/ELR_np.png'), width=1000)
 plot(scan.np.em.b)
dev.off()

png(paste0('~/public_html/ELR_scan.norm.mr.png'))
 plot(scan.norm.mr)
dev.off()

png(paste0('~/public_html/ELR_full.norm.add_only.png'))
 plot(full.norm.add_only)
dev.off()
