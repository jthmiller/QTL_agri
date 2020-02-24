#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
library('qtl')
##library('parallel')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################
##load(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
mapfile <- "NBH_2172_imputed_high_confidence_tsp_mapped.csv"
filename <- file.path(mpath,mapfile)
cross <- read.cross(file=mapfile , format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
cross$pheno <- as.data.frame(cross$pheno)
################################################################################

perm_count <- as.numeric(commandArgs(TRUE)[3])
arraynum <- as.numeric(commandArgs(TRUE)[5])
cores <- as.numeric(commandArgs(TRUE)[4])
batch <- round(nind(cross)/2)

print(commandArgs(TRUE))
print(paste('pop =',pop,', perm = ',perm_count,', cores =', cores,', array =',arraynum))
set.seed(arraynum)

################################################################################

print(paste(cores,'cores'))
erp <- 0.001
sex.phen <- pull.pheno(cross, "sex")
names(cross$geno) <- ifelse(names(cross$geno) == "5","X",names(cross$geno))
attr(cross$geno[["X"]], 'class') <- 'X'
##cross <- drop.markers(cross, markernames(cross,'X'))

################################################################################
if(pop == 'ELR'){
 cross <- subset(cross, chr=c(1:4,6:17,19:24))
} else {
 cross <- subset(cross, chr=c(1,3:4,6:24))
}
cov <- ifelse(pop == 'ELR',18,2)
################################################################################

################################################################################
so <- summary(scanone(cross,pheno.col=4, model="normal", method="imp"))[cov,]
mar <- find.marker(cross, so$chr, so$lod)
g <- pull.geno(fill.geno(cross))[,mar]
g <- cbind(as.numeric(g==1), as.numeric(g==2))
summary(scanone(cross,pheno.col=4, model="normal", method="imp",addcovar=g))
################################################################################

norm.imp.perms.2 <- scantwo(cross, pheno.col=1, model="normal", method="imp",
 incl.markers=T, clean.output=T, clean.nmar=100, clean.distance=100,
 assumeCondIndep=T, n.cluster=cores,  addcovar=g, n.perm=perm_count)

norm.imp.perms.pens <- calc.penalties(norm.imp.perms.2, alpha=0.1)

norm.imp.perms.1 <- scanone(cross, pheno.col=1, model='normal', method = "imp",
 n.perm = 10000, n.cluster=cores, addcovar=g)

lod <- summary(norm.imp.perms.1)[1]

summary(norm.imp.perms.2)
summary(norm.imp.perms.pens)
print(norm.imp.perms.pens)
print(paste('done with', perm_count, 'scan 2 permutations'))
print(summary(norm.imp.perms.1))
################################################################################
save.image(file.path(mpath,paste0(pop,'_scan_perms_norm_imp.rsave')))
################################################################################
