#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]
library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################
## SCANTWO ON SUBSET

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

cross <- sim.geno(cross)
cross <- calc.genoprob(cross,step=1,error.prob=0.01,off.end=5)

gg <- sim.geno(cross,step=4)
gg_step2 <- reduce2grid(gg)

################################################################################

full.norm.add_only <- stepwiseqtl(gg_step2, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=5)

#add_int <- refineqtl(gg_step2, qtl= full.norm.add_only, formula=y ~ Q1+Q2+Q3+Q4, method="imp", model='binary')

#qtl_interaction <- addint(gg_step2, qtl=add_int, pheno.col=4, method="imp", formula=y~Q1+Q2+Q3+Q4, model='binary')

################################################################################

## PERMS 5% LOD 4.07
sc2_normal_imp <- scantwo(gg_step2, pheno.col=5, model="normal",
             method="imp",
             addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs",
             incl.markers=FALSE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=12)

sc2_normal_imp_perms <- scantwo(gg_step2, pheno.col=5, model="normal",
             method="imp", addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs", n.perm=1000,
             incl.markers=FALSE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=12)

save.image(file.path(mpath,'scantwo.scans.new.short.rsave'))
