#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
dens <- as.numeric(commandArgs(TRUE)[3])

print(dens)

library('qtl')
##library('parallel')
library('snow')

source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)
##cores <- detectCores() - 2

################################################################################
################################################################################
erp <- 0.0025
cores <- 22
################################################################################
################################################################################
fl <- paste0(pop,'_imp.mapped.tsp.csv')
fl <- file.path(mpath,fl)

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno),5)
cross <- jittermap(cross)

cross <- removeDoubleXO(cross)

gt <- geno.table(cross)
mis <- 5
bfixA <- rownames(gt[which(gt$missing > mis),])
cross <- drop.markers(cross, bfixA)
print(paste('markers dropped due to missing =',length(bfixA)))

dups <- findDupMarkers(cross, exact.only=F, adjacent.only=T)
dups <- unlist(dups)
print(paste('markers dropped due to duplicate =',length(dups)))
cross <- drop.markers(cross, dups)
cross <- removeDoubleXO(cross)

png(paste0('~/public_html/',pop,'_gts_dwsmpl.png'),height=2500,width=4500)
 geno.image(cross, reorder=1, cex=2)
dev.off()


cross_imp <- fill.geno(cross, method="maxmarginal", error.prob = 0.01, min.prob=0.99)

png(paste0('~/public_html/',pop,'_gts_all.png'),height=2500,width=4500)
 geno.image(cross_imp, reorder=1, cex=2)
dev.off()

################################################################################
fl <- file.path(mpath,paste0(pop,'_downsampled.csv'))
write.cross(cross,filestem=fl,format="csv")
################################################################################

################################################################################
save.image(file.path(mpath,paste0(pop,'_downsampled.rsave')))
################################################################################
