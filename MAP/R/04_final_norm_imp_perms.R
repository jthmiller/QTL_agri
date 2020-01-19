pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
perm_count <- commandArgs(TRUE)[is.numeric(commandArgs(TRUE))]

library('qtl')
library('parallel')
library('snow')

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)
cores <- detectCores() - 2

################################################################################
load(file.path(mpath,paste0(pop,'_norm_imp.rsave')))
################################################################################

perms <- scantwo(gg_step2, incl.markers=F, chr = c(1:4,6:24), pheno.col=5, model="normal", method="imp", clean.output=T, clean.nmar=10,clean.distance=10,n.perm=1,assumeCondIndep=T,n.cluster=22)
pens <- calc.penalties(perms, alpha=0.10)
print(paste('done with', perm_count, 'permutations'))

################################################################################
save.image(file.path(mpath,paste0(pop,'_perms_norm_imp.rsave')))
################################################################################

full.norm.imp <- stepwiseqtl(gg_step2, penalties=pens, incl.markers=F, qtl=norm.add.qtls, additive.only = F, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=8, chr=qtls_chr)
print(paste('done with', perm_count, 'full.norm.imp'))

################################################################################
save.image(file.path(mpath,paste0(pop,'_perms_norm_imp.rsave')))
################################################################################

full.norm.hk <- stepwiseqtl(gg_step2, incl.markers=T, additive.only = F, model='normal', method = "hk", pheno.col = 5, scan.pairs = T, max.qtl=8)
print(paste('done with', perm_count, 'full.norm.hk'))

################################################################################
save.image(file.path(mpath,paste0(pop,'_perms_norm_imp.rsave')))
################################################################################
