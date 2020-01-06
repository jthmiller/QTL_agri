#!/bin/R

#pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]

pop <- 'ELR'

library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

# single scans
load(file.path(mpath,'single_scans.elr_missing.rsave'))
#short
load(file.path(mpath,'scantwo.scans.elr.short.missing.rsave'))
################################################################################
## 2 additive qtls
qtls <- summary(full.norm.add_only)
qtls <- makeqtl(gg_step2, chr=as.character(qtls$chr), pos=as.numeric(qtls$pos), what="draws")
qtls <- refineqtl(gg_step2, qtl=qtls, pheno.col=4, formula=y ~ Q1+Q2, method="imp", model='binary')

## DOES 13 or 18 interact with a third qtl or one another?
qtls_int <- addint(gg_step2, qtl=qtls, pheno.col=4, method="imp", formula=y~Q1+Q2, model='binary')
qtls_int_13 <- addqtl(gg_step2, qtl=qtls, formula=y~Q1+Q2+Q1*Q3)
qtls_int_18 <- addqtl(gg_step2, qtl=qtls, formula=y~Q1+Q2+Q2*Q3)



AHR2a_del
qtls_interaction <- addint(gg_step2, qtl=qtls, pheno.col=4, method="imp", formula=y~Q1+Q2+Q1*Q3, model='binary')
qtls_interaction <- addint(gg_step2, qtl=qtls, pheno.col=4, method="imp", formula=y~Q1+Q2, model='binary')


## Scan with the AHR locus as a QTL
qc <- c(1, 13, 18)
qp <- c(0, 24, 48)
qtls_int <- makeqtl(gg_step2, chr=qc, pos=qp, what="draws")
qtls_int <- refineqtl(gg_step2, qtl=qtls_int, pheno.col=4, formula=y ~ Q1+Q2+Q1*Q3, method="imp", model='binary')
ahr_int <- addint(gg_step2, qtl=qtls_int, pheno.col=4, method="imp", formula=y~Q1+Q2+Q3, model='binary')
################################################################################

add_int <- refineqtl(gg_step2, qtl= full.norm.add_only, formula=y ~ Q1+Q2, method="imp", model='binary')
qtl_interaction <- addint(gg_step2, qtl=add_int, pheno.col=4, method="imp", formula=y~Q1+Q2, model='binary')
################################################################################

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
