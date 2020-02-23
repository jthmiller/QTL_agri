32
#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
cores <- as.numeric(commandArgs(TRUE)[3])

print(commandArgs(TRUE))
print(paste(pop))

library('qtl')
##library('parallel')
library('snow')

source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")

mpath <- '/home/jmiller1/QTL_agri/data'
#fl <- paste0(pop,'_imp.mapped.tsp.csv')
#fl <- file.path(mpath,fl)


################################################################################
##load(file.path(mpath,paste0(pop,'_downsampled.rsave')))
################################################################################
################################################################################
################################################################################

##TOO MANY MARKERS
##fl <- paste0(pop,'_imp.mapped.tsp.csv')

#fl <- file.path(paste0(pop,'_downsampled.csv'))

fl <- 'NBH_2172_imputed_high_confidence_tsp_mapped.csv'
cross <- read.cross(file=fl , format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
#cores <- 20
################################################################################

print(paste(cores,'cores'))
erp <- 0.001
sex.phen <- pull.pheno(cross, "sex")
names(cross$geno) <- ifelse(names(cross$geno) == "5","X",names(cross$geno))
attr(cross$geno[["X"]], 'class') <- 'X'

################################################################################

################################################################################

bin.em.2 <- scantwo(cross, pheno.col=4, model="binary", method="mr",
 clean.output=T, clean.nmar=50, clean.distance=50, maxit=1000,
 assumeCondIndep=T, n.cluster=cores, use="complete.obs")

################################################################################
save.image(file.path(mpath,paste0(pop,'_scan2_bin_mr.rsave')))
################################################################################
