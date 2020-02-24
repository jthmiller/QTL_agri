#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
perm_count <- as.numeric(commandArgs(TRUE)[3])
cores <- as.numeric(commandArgs(TRUE)[4])

print(commandArgs(TRUE))
print(paste(pop,perm_count))

library('qtl')
##library('parallel')
library('snow')

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################
################################################################################
mapfile <- "NBH_2172_imputed_high_confidence_tsp_mapped.csv"
filename <- file.path(mpath,mapfile)
cross <- read.cross(file=mapfile , format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
cross$pheno <- as.data.frame(cross$pheno)
cross <- calc.genoprob(cross, step=0, off.end=0, error.prob=0.001, map.function="kosambi", stepwidth="fixed")

print(paste(cores,'cores'))
erp <- 0.0025

################################################################################

################################################################################
load(file.path(mpath,paste0(pop,'_all_perms_bin_em.rsave')))
## perms = scan2 permutations
## perms_1 = scanone permutations
## bin.em.perms.pens = 0.1 penalties
load(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
## sone = scanone
## bin.em.2 = scantow
## addcovar = g
################################################################################
##bin.add.em <- scanone(cross, pheno.col=5, model='binary', method = "em")
#qtl <- summary(bin.add.em,lod)
#qtl <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="draws")
#qtl
top_2 <- order(so$lod,decreasing =T)[1:2]

qtl <- makeqtl(cross, chr=so[top_2,'chr'], pos=so[top_2,'pos'], what="prob")
qtl

pens <- calc.penalties(perms.2, alpha=0.05)

full.bin.em.step <- stepwiseqtl(cross, model='binary', method = "hk", pheno.col = 4,
 penalties=pens, incl.markers=T, qtl=qtl, additive.only = T,  scan.pairs = F, max.qtl=8)

summary(full.bin.em.step)
################################################################################
save.image(file.path(mpath,paste0(pop,'_step_bin_em.rsave')))
################################################################################
