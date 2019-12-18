#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]

pop <- 'ELR'

library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)


# single scans
load(file.path(mpath,'single_scans.elr_missing.rsave'))
#short
load(file.path(mpath,'scantwo.scans.elr.short.rsave'))

## Scan with the AHR locus
qc <- c(1, 13, 18)
qp <- c(0.635, 58, 46)
qtl <- makeqtl(gg_step2, chr=qc, pos=qp, what="draws")
ahr_interaction <- addint(gg_step2, qtl=qtl, pheno.col=4, method="imp", formula=y~Q1+Q2+Q3, model='binary')
################################################################################

add_int <- refineqtl(gg_step2, qtl= full.norm.add_only, formula=y ~ Q1+Q2, method="imp", model='binary')
qtl_interaction <- addint(gg_step2, qtl=add_int, pheno.col=4, method="imp", formula=y~Q1+Q2, model='binary')


################################################################################
qc <- c(1, 13, 18)
qp <- c(0.635, 58, 46)
qtl <- makeqtl(gg_step2, chr=qc, pos=qp, what="draws")
ahr_interaction <- addint(gg_step2, qtl=qtl, pheno.col=4, method="imp", formula=y~Q1+Q2+Q3, model='binary')
################################################################################




png(paste0('~/public_html/ELR_model.png'))
 plotModel(full.norm.add_only)
dev.off()


full.norm.add_only
