#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]
library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

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

################################################################################

sc2_normal_imp <- scantwo(gg_step2, pheno.col=5, model="normal",
             method="imp",
             addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs",
             incl.markers=FALSE, clean.output=TRUE,
             clean.nmar=20, clean.distance=20,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=TRUE, batchsize=250, n.cluster=12)

sc2_normal_imp_perms <- scantwo(gg_step2, pheno.col=5, model="normal",
             method="imp", addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs", n.perm=1000,
             incl.markers=FALSE, clean.output=TRUE,
             clean.nmar=20, clean.distance=20,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=12)

sc2_normal_imp_penalties <- calc.penalties(sc2_normal_imp_perms, alpha=0.1, lodcolumn=1)

## sm <- summary(object, thresholds,what=c("best", "full", "add", "int"),
##             perms=sc2_normal_imp_perms, alphas, lodcolumn=1,
##             pvalues=FALSE, allpairs=TRUE)
##
save.image(file.path(mpath,'scantwo.scans.elr.long.rsave'))

sc2_bin_em <- scantwo(gg_step2, pheno.col=4, model="binary",
             method="em",addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs",
             incl.markers=FALSE, clean.output=TRUE,
             clean.nmar=20, clean.distance=20,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=12)

sc2_bin_em_perms <- scantwo(gg_step2, pheno.col=4, model="binary",
             method="em",addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs", n.perm=1000,
             incl.markers=FALSE, clean.output=TRUE,
             clean.nmar=20, clean.distance=20,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=12)

sc2_bin_em_penalties <- calc.penalties(sc2_bin_em_perms, alpha=0.1, lodcolumn=1)

save.image(file.path(mpath,'scantwo.scans.elr.long.rsave'))

################################################################################
full.bin.add_only <- stepwiseqtl(gg_step2, additive.only = T, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=4)
add.norm <- stepwiseqtl(gg_step2, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=6)
full.norm <- stepwiseqtl(gg_step2, incl.markers=F, scan.pairs=T, keeptrace=T, additive.only = F, model='normal', method = "imp", pheno.col = 5, max.qtl=6)

save.image(file.path(mpath,'pairwise_scans.elr.rsave'))
################################################################################

fit.full.norm <- fitqtl(gg_step2,qtl=full.norm, pheno.col=4)

save.image(file.path(mpath,'pairwise_scans.elr.rsave'))
