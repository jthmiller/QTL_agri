#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
dens <- as.numeric(commandArgs(TRUE)[3])

print(dens)

library('qtl')
##library('parallel')
library('snow')

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)
##cores <- detectCores() - 2

################################################################################
################################################################################
erp <- 0.0025
cores <- 22
cores <- 5
################################################################################
################################################################################
## Read cross
cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno),5)
cross <- jittermap(cross)

gt <- geno.table(cross)
mis <- 5
bfixA <- rownames(gt[which(gt$missing < mis),])

cross <- pull.markers(cross, bfixA)
cross_imp <- fill.geno(cross, method="maxmarginal", error.prob = 0.08, min.prob=0.995)



## Estimate gt prob and impute before downsample
dups <- findDupMarkers(cross, exact.only=T, adjacent.only=T)
#dups <- findDupMarkers(cross, exact.only=F, adjacent.only=F)
dups <- names(dups)
if(pop == 'ELR.missing') dups <- c(dups,"AHR2a_del")
cross <- drop.markers(cross, dups)
cross

crossbk <- cross

cross_imp <- fill.geno(cross, method="argmax", error.prob = 0.08)



cross_imp <- fill.geno(cross, method="maxmarginal", error.prob = 0.08, min.prob=0.998)
cross_imp <- fill.geno(cross_imp, method="no_dbl_XO", error.prob = 0.08)
cross_imp <- fill.geno(cross_imp, method="maxmarginal", error.prob = 0.08, min.prob=0.995)

png(paste0('~/public_html/',pop,'_gts_all.png'),height=2500,width=4500)
 geno.image(cross_imp, reorder=1, cex=2)
dev.off()



################################################################################
fl <- file.path(mpath,paste0(pop,'_downsampled'))
write.cross(cross,filestem=fl,format="csv")
################################################################################

cross$pheno <- as.data.frame(cross$pheno)

cross <- sim.geno(cross, stepwidth="fixed", step=dens,  error.prob=erp, off.end=1, map.function="kosambi", n.draws=100)
cross <- calc.genoprob(cross, stepwidth="fixed", step=dens, error.prob=erp, off.end=1, map.function="kosambi")
cross <- argmax.geno(cross, stepwidth="fixed", step=dens, error.prob=erp, off.end=1, map.function="kosambi")

cross <- reduce2grid(cross)

(summary(pull.map(cross))['overall','length']) / (length(colnames(pull.genoprob(cross)))/3)
print('markers per CM')

################################################################################
save.image(file.path(mpath,paste0(pop,'_downsampled.rsave')))
################################################################################
