#!/bin/R
### first run combine pops for multi-pop cross objects
pop <- 'NBH'
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
library("ggridges")
library("plyr")
library("scales")
library("ggrepel")
library('qtl')
library('RColorBrewer')
################################################################################
mpath <- '/home/jmiller1/QTL_agri/data'
setwd(mpath)
##load(file.path(mpath,'08_phys_plots_pos.rsave'))
####################################################################################
fl <- paste0('ELR.mapped.tsp.csv')
fl <- file.path(mpath,fl)
cross_ELR <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)
cross_ELR$pheno$pheno_norm <- round(nqrank(cross_ELR$pheno$Pheno),5)
#gt.elr <- geno.table(cross_ELR)

fl <- paste0('NBH.mapped.tsp.csv')
fl <- file.path(mpath,fl)
cross_NBH <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)
cross_NBH$pheno$pheno_norm <- round(nqrank(cross_NBH$pheno$Pheno))

#gt.nbh <- geno.table(cross_NBH)
################################################################################

################################################################################
#### AHRs #####
AHR.bed <- read.table(file.path(mpath,"lift_AHR_genes.bed"), stringsAsFactors = F, header = F)
colnames(AHR.bed) <- c("chrom", "str", "stp", "gene")
AHR.bed$chrom <- as.numeric(gsub("chr", "", AHR.bed$chrom))
AHR.bed$str <- as.numeric(AHR.bed$str)
AHR.bed$stp <- as.numeric(AHR.bed$stp)
AHR.notmap <- AHR.bed[is.na(AHR.bed$chrom), ]
AHR.bed <- AHR.bed[!is.na(AHR.bed$chrom), ]
AHR.bed$gene <- gsub(":158640", "", AHR.bed$gene)
# add arnts (forgot to scan for them)
################################################
nbh.gens <- cnv.ahrs(cross_NBH, AHRdf = AHR.bed, EXP = F)
#new.gens <- cnv.ahrs(cross.new, AHRdf = AHR.bed, EXP = F)
elr.gens <- cnv.ahrs(cross_ELR, AHRdf = AHR.bed, EXP = F)
#brp.gens <- cnv.ahrs(cross.brp, AHRdf = AHR.bed, EXP = F)
################################################
cands <- c("AHR1","aip","ARNT","ARNT2","ahrr","ahr1b","AHR2b")
ahr_nbh <- nbh.gens[which(nbh.gens$gene %in% cands),]
ahr_elr <- nbh.gens[which(elr.gens$gene %in% cands),]
################################################################################
### ggplot popgen locations
nbh.popgen <- read.table(file.path(mpath,"outliersNBH.txt.ncbi.lifted"), sep = "\t", header = T)
new.popgen <- read.table(file.path(mpath,"outliersNYC.txt.ncbi.lifted"), sep = "\t", header = T)
elr.popgen <- read.table(file.path(mpath,"outliersER.txt.ncbi.lifted"), sep = "\t", header = T)
brp.popgen <- read.table(file.path(mpath,"outliersBP.txt.ncbi.lifted"), sep = "\t", header = T)
################################################################################
### Use nbh coords but elr and new popgen
#new.rank <- cnv.popgen(cross.nbh, new.popgen, top = 50)
nbh.rank <- cnv.popgen(cross_NBH, nbh.popgen, top = 75)
elr.rank <- cnv.popgen(cross_ELR, elr.popgen, top = 75)
#brp.rank <- cnv.popgen(cross.nbh, brp.popgen, top = 50)
################################################################################

## ALL GENES
genes.bed <- read.table(file.path(mpath,"lifted_genes.bed"), stringsAsFactors = F, header = T)
genes.bed$chr <- gsub('chr','',genes.bed$chr)
genes.bed <- genes.bed[genes.bed$chr %in% c(1:24),]
genes.bed$mid <- round(apply(genes.bed[c('start','end')],1,mean))

nbh_gene_models <- conv_popstat(cross_NBH, popgen=genes.bed, whichcol='start',newname='cm_start')
nbh_gene_models$cm_end <- conv_popstat(cross_NBH, popgen=genes.bed, whichcol='end',newname='cm_end')[,'cm_end']
nbh_gene_models$cm_mid <- conv_popstat(cross_NBH, popgen=genes.bed, whichcol='end',newname='cm_mid')[,'cm_mid']

elr_gene_models <- conv_popstat(cross_ELR, popgen=genes.bed, whichcol='start',newname='cm_start')
elr_gene_models$cm_end <- conv_popstat(cross_ELR, popgen=genes.bed, whichcol='end',newname='cm_end')[,'cm_end']
elr_gene_models$cm_mid <- conv_popstat(cross_ELR, popgen=genes.bed, whichcol='end',newname='cm_mid')[,'cm_mid']

## get_genes_cm(chr=2, start = 32809365,stop = 32962365, models = nbh_gene_models, colm = 'start')
## get_genes_cm(chr=1, start = 20,stop = 30,models = nbh_gene_models, colm = 'cm_mid')
 get_genes_cm(chr=1, start = 3000000,stop = 5000000,models = nbh_gene_models, colm = 'start')
################################################################################
################################################################################

pbs <- file.path(mpath, 'pbstat.txt.ncbi.lifted')
pbs <- read.table(pbs, sep = "\t", header = T)
pbs$mid <- pbs$V2 + (abs(pbs$V3 - pbs$V2) * .5)
pbs$V1 <- gsub('chr',"",pbs$V1)
colnames(pbs)[1:3] <- c('chr','start','end')
pbs <- pbs[!is.na(as.numeric(pbs$chr)),]
pbs$nbh_cm <- conv_popstat(cross_NBH, popgen=pbs, whichcol='mid',newname='nbh_cm')$nbh_cm
pbs$elr_cm <- conv_popstat(cross_ELR, popgen=pbs, whichcol='mid',newname='elr_cm')$elr_cm

pfst <- file.path(mpath, 'pfst.txt.ncbi.lifted')
pfst <- read.table(pfst, sep = "\t", header = T)
pfst$mid <- pfst$start + (abs(pfst$end - pfst$start) * .5)
pfst$Scaffold <- gsub('chr',"",pfst$Scaffold)
colnames(pfst)[1] <- 'chr'
pfst <- pfst[!is.na(as.numeric(pfst$chr)),]
pfst$nbh_cm <- conv_popstat(cross_NBH, popgen=pfst, whichcol='mid',newname='nbh_cm')$nbh_cm
pfst$elr_cm <- conv_popstat(cross_ELR, popgen=pfst, whichcol='mid',newname='elr_cm')$elr_cm

###taj <- file.path(mpath, 'tajstat.txt.ncbi.lifted')
###taj <- read.table(taj, sep = "\t", header = T)
###taj$mid <- taj$start + (abs(taj$end - taj$start) * .5)
###taj$Scaffold <- gsub('chr',"",taj$Scaffold)
###colnames(taj)[1] <- 'chr'
###taj <- conv_popstat(cross.nbh, popgen=taj, whichcol='mid',newname='nbh_mp')
#####taj <- conv_popstat(cross.elr, popgen=taj, whichcol='mid',newname='elr_mp')
###
###pi <- file.path(mpath, 'piper.txt.ncbi.lifted')
###pi <- read.table(pi, sep = "\t", header = T)
###pi$mid <- pi$start + (abs(pi$end - pi$start) * .5)
###pi$Scaffold <- gsub('chr',"",pi$Scaffold)
###colnames(pi)[1] <- 'chr'
###pi <- conv_popstat(cross.nbh, popgen=pi, whichcol='mid',newname='nbh_mp')
#####pi <- conv_popstat(cross.elr, popgen=pi, whichcol='mid',newname='elr_mp')
###
################################################################################
################################################################################

cross_NBH <- sim.geno(cross_NBH, step=1, error.prob=0.025, off.end=5, map.function="kosambi", n.draws=160)
cross_ELR <- sim.geno(cross_ELR, step=1, error.prob=0.025, off.end=5, map.function="kosambi", n.draws=160)

cross_NBH <- calc.genoprob(cross_NBH, step=1, error.prob=0.025, off.end=5, map.function="kosambi")
cross_ELR <- calc.genoprob(cross_ELR, step=1, error.prob=0.025, off.end=5, map.function="kosambi")

#scan_nbh <- scanone(cross_NBH, method = "mr", model = "binary", pheno.col = 4)
#scan_elr <- scanone(cross_ELR, method = "mr", model = "binary", pheno.col = 4)
scan_nbh <- scanone(cross_NBH, method = "imp", model = "binary", pheno.col = 4)
scan_elr <- scanone(cross_ELR, method = "imp", model = "binary", pheno.col = 4)
################################################################################

################################################################################
## plot_ef in rQTL2 ############################################################
col <- c("slateblue", "violetred", "green3")

nbh <- convert2cross2(cross_NBH)
nbh_map <- insert_pseudomarkers(nbh$gmap, step=1)
nbh_pr <- calc_genoprob(nbh, nbh_map, error_prob=0.025, cores=4)

elr <- convert2cross2(cross_ELR)
elr_map <- insert_pseudomarkers(elr$gmap, step=1)
elr_pr <- calc_genoprob(elr, elr_map, error_prob=0.025, cores=4)

################################################################################
################################################################################

################################################################################
## Correlate lod and segregation distortion

elr_c2eff <- lapply(1:24,function(X) {
 scan1coef(elr_pr[,as.character(X)], elr$pheno[,"bin"])
})
elr_c2eff <- do.call(rbind,elr_c2eff)
elr_seg <- geno.table(cross_ELR)[rownames(elr_c2eff),'P.value']

nbh_c2eff <- lapply(1:24,function(X) {
 scan1coef(nbh_pr[,as.character(X)], nbh$pheno[,"bin"])
})
nbh_c2eff <- do.call(rbind,nbh_c2eff)
nbh_seg <- geno.table(cross_NBH)[rownames(nbh_c2eff),'P.value']

################################################################################
################################################################################

## rf interaction plot #########################################################
################################################################################
a <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross_ELR)[[X]], 0.25)} ))
a <- pull.markers(cross_ELR,a)
a <- subset(a,ind=a$pheno$bin==0)
rf <- est.rf(a, maxit=10000, tol=1e-6)

png("/home/jmiller1/public_html/elr_rf.png", width=1000, height=1000)
 plotRF(rf, zmax=8, col.scheme="redblue",what='lod')
dev.off()

################################################################################
a <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross_NBH)[[X]], 0.25)} ))
a <- pull.markers(cross_NBH,a)
#a <- subset(a,ind=a$pheno$bin==0)
rf <- est.rf(a, maxit=10000, tol=1e-6)

png("/home/jmiller1/public_html/nbh_rf.png", width=1000, height=1000)
 plotRF(rf, zmax=10, col.scheme="redblue",what='lod')
dev.off()
################################################################################
################################################################################
