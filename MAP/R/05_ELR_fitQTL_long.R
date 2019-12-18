#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]
library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

load(file.path(mpath,'single_scans.elr.rsave'))
################################################################################

################################################################################
## SCANTWO ON SUBSET (this map made in tspmap.R

fl <- file.path(mpath,'elr.mapped.tsp.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE)

cross <- sim.geno(cross)
cross <- calc.genoprob(cross,step=1,error.prob=0.01,off.end=5)

gg <- sim.geno(cross,step=4)
gg_step2 <- reduce2grid(gg)

## Scan with the AHR locus
qc <- c(1, 13, 18)
qp <- c(0.635, 58, 46)
qtl <- makeqtl(gg_step2, chr=qc, pos=qp, what="draws")
ahr_interaction <- addint(gg_step2, qtl=qtl, pheno.col=4, method="imp", formula=y~Q1+Q2+Q3, model='binary')
################################################################################

add_int <- refineqtl(gg_step2, qtl= full.norm.add_only, formula=y ~ Q1+Q2, method="imp", model='binary')
qtl_interaction <- addint(gg_step2, qtl=add_int, pheno.col=4, method="imp", formula=y~Q1+Q2, model='binary')

################################################################################

sc2_bin_em <- scantwo(gg_step2, pheno.col=4, model="binary",
             method="em",addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs",
             incl.markers=FALSE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=6)

sc2_bin_em_perms <- scantwo(gg_step2, pheno.col=4, model="binary",
             method="em",addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs", n.perm=1000,
             incl.markers=FALSE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=4)

## sm <- summary(object, thresholds,what=c("best", "full", "add", "int"),
##             perms=sc2_normal_imp_perms, alphas, lodcolumn=1,
##             pvalues=FALSE, allpairs=TRUE)
##
save.image(file.path(mpath,'scantwo.scans.elr.long.rsave'))

################################################################################
full.bin.add_only <- stepwiseqtl(gg_step2, additive.only = T, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=4)
add.norm <- stepwiseqtl(gg_step2, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=6)
full.norm <- stepwiseqtl(gg_step2, incl.markers=F, scan.pairs=T, keeptrace=T, additive.only = F, model='normal', method = "imp", pheno.col = 5, max.qtl=6)

save.image(file.path(mpath,'pairwise_scans.elr.rsave'))
################################################################################

fit.full.norm <- fitqtl(gg_step2,qtl=full.norm, pheno.col=4)

save.image(file.path(mpath,'pairwise_scans.elr.rsave'))
