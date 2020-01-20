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

mpath <- '/home/jmiller1/QTL_agri/data'
setwd(mpath)

## Phenotypes
################################################
cross_ELR <- read.cross(format = "csv", dir = mpath, file = 'ELR.mapped.tsp.csv', genotypes=c("1","2","3"), estimate.map = FALSE)
cross_NBH <- read.cross(format = "csv", dir = mpath, file = 'NBH.mapped.tsp.csv', genotypes=c("1","2","3"), estimate.map = FALSE)

fl <- file.path('BRP_unmapped_filtered.csv')
cross_BRP <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
#cross_BRP <- read.cross(format = "csv", dir = mpath, file = 'BRP.mapped.tsp.csv', genotypes=c("1","2","3"), estimate.map = FALSE)

fl <- file.path('NEW_unmapped_filtered.csv')
cross_NEW <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
#cross_NEW <- read.cross(format = "csv", dir = mpath, file = 'NEW.mapped.tsp.csv', genotypes=c("1","2","3"), estimate.map = FALSE)

################################################################################
remap_orders <- function(cross_tomap, cross_usemap){

 pmap_names <- lapply(pull.map(cross_tomap), function(X) {
    names(X)
  })

 usemap_names <- lapply(pull.map(cross_usemap), function(X) {
    names(X)
  })

 mapply(match,pmap_names,usemap_names)
}
################################################################################

NEW_orders <- remap_orders(cross_NEW,cross_NBH)
NEW_orders <- lapply(NEW_orders,order)

BRP_orders <- remap_orders(cross_BRP,cross_NBH)
BRP_orders <- lapply(BRP_orders,order)

reorder_BRP <- cross_BRP
reorder_NEW <- cross_NEW

 for (i in 1:24){
  i <- as.character(i)
  reorder_NEW <<- switch.order(reorder_NEW, i, NEW_orders[[i]])
  reorder_BRP <<- switch.order(reorder_BRP, i, BRP_orders[[i]])
 }

NBH_NEW <- pull.markers(cross_NBH, markernames(reorder_NEW))
NBH_BRP <- pull.markers(cross_NBH, markernames(reorder_BRP))

reorder_NEW <- replace.map(reorder_NEW, pull.map(NBH_NEW))
reorder_BRP <- replace.map(reorder_BRP, pull.map(NBH_BRP))

################################################################################
################################################################################
cross.new <- sim.geno(reorder_NEW, n.draws = 500, step = 5, off.end = 10, error.prob = 0.1,
  map.function = "kosambi", stepwidth = "fixed")

cross.brp <- sim.geno(reorder_BRP, n.draws = 500, step = 5, off.end = 10, error.prob = 0.1,
  map.function = "kosambi", stepwidth = "fixed")

cross.nbh <- sim.geno(cross_NBH, n.draws = 500, step = 5, off.end = 10, error.prob = 0.025,
  map.function = "kosambi", stepwidth = "fixed")

cross.elr <- sim.geno(cross_ELR, n.draws = 500, step = 5, off.end = 10, error.prob = 0.025,
  map.function = "kosambi", stepwidth = "fixed")
################################################################################
################################################################################

################################################################################
##keeping colors consistent####################
all.pops <- c("NBH", "BRP", "ELR", "NEW")
popcol <- brewer.pal(8, "Paired")[c(2, 4, 6, 8)]
names(popcol) <- all.pops
popcol <- popcol[c('NBH','BRP','NEW','ELR')]

popgen <- popcol
names(popgen) <- c('NBH','BP','NYC','ER')

popout <- c(popgen,'grey')
names(popout) <- c('NBH','BP','NYC','ER','BI')

### Color for stat comparisons
statcol <- popcol
names(statcol) <- c('BI.NBH','BP.F','NYC.SH','ER.KC')
################################################

pbs <- file.path(mpath, 'pbstat.txt.ncbi.lifted')
pbs <- read.table(pbs, sep = "\t", header = T)
pbs$mid <- pbs$V2 + (abs(pbs$V3 - pbs$V2) * .5)
pbs$V1 <- gsub('chr',"",pbs$V1)
colnames(pbs)[1:3] <- c('chr','start','end')
pbs <- conv_popstat(cross.nbh, popgen=pbs, whichcol='mid',newname='nbh_mp')
#pbs <- conv_popstat(cross.elr, popgen=pbs, whichcol='mid',newname='elr_mp')

pfst <- file.path(mpath, 'pfst.txt.ncbi.lifted')
pfst <- read.table(pfst, sep = "\t", header = T)
pfst$mid <- pfst$start + (abs(pfst$end - pfst$start) * .5)
pfst$Scaffold <- gsub('chr',"",pfst$Scaffold)
colnames(pfst)[1] <- 'chr'
pfst <- conv_popstat(cross.nbh, popgen=pfst, whichcol='mid',newname='nbh_mp')
#pfst <- conv_popstat(cross.elr, popgen=pfst, whichcol='mid',newname='elr_mp')

taj <- file.path(mpath, 'tajstat.txt.ncbi.lifted')
taj <- read.table(taj, sep = "\t", header = T)
taj$mid <- taj$start + (abs(taj$end - taj$start) * .5)
taj$Scaffold <- gsub('chr',"",taj$Scaffold)
colnames(taj)[1] <- 'chr'
taj <- conv_popstat(cross.nbh, popgen=taj, whichcol='mid',newname='nbh_mp')
##taj <- conv_popstat(cross.elr, popgen=taj, whichcol='mid',newname='elr_mp')

pi <- file.path(mpath, 'piper.txt.ncbi.lifted')
pi <- read.table(pi, sep = "\t", header = T)
pi$mid <- pi$start + (abs(pi$end - pi$start) * .5)
pi$Scaffold <- gsub('chr',"",pi$Scaffold)
colnames(pi)[1] <- 'chr'
pi <- conv_popstat(cross.nbh, popgen=pi, whichcol='mid',newname='nbh_mp')
##pi <- conv_popstat(cross.elr, popgen=pi, whichcol='mid',newname='elr_mp')


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

#################################################
#scan.norm.imp.NBH <- scanone(cross.nbh, method = "imp", model = "normal", pheno.col = 5)
#scan.norm.imp.ELR <- scanone(cross.elr, method = "imp", model = "normal", pheno.col = 5)
#scan.norm.imp.NEW <- scanone(cross.new, method = "imp", model = "normal", pheno.col = 5)
#scan.norm.imp.BRP <- scanone(cross.brp, method = "imp", model = "normal", pheno.col = 5)
#################################################
#### marker regression plots
scan.bin.mr.NBH <-  scanone(cross.nbh, method = "mr", model = "binary", pheno.col = 4)
scan.bin.mr.ELR <-  scanone(cross.elr, method = "mr", model = "binary", pheno.col = 4)
scan.bin.mr.BRP <-  scanone(cross.brp, method = "mr", model = "binary", pheno.col = 4)
scan.bin.mr.NEW <-  scanone(cross.new, method = "mr", model = "binary", pheno.col = 4)
#################################################

### use scanone for plots
scan.bin.imp.NBH <-  scanone(cross.nbh, method = "em", model = "binary", pheno.col = 4)
scan.bin.imp.ELR <-  scanone(cross.elr, method = "em", model = "binary", pheno.col = 4)
scan.bin.imp.NEW <-  scanone(cross.new, method = "em", model = "binary", pheno.col = 4)
scan.bin.imp.BRP <-  scanone(cross.brp, method = "em", model = "binary", pheno.col = 4)

themelt.nbh <- scan.bin.imp.NBH
themelt.new <- scan.bin.imp.NEW
themelt.elr <- scan.bin.imp.ELR
themelt.brp <- scan.bin.imp.BRP

themelt.nbh$pop <- "NBH"
themelt.new$pop <- "NEW"
themelt.elr$pop <- "ELR"
themelt.brp$pop <- "BRP"

themelt.brp.mr <- scan.bin.mr.BRP
themelt.elr.mr <- scan.bin.mr.ELR
themelt.nbh.mr <- scan.bin.mr.NBH
themelt.new.mr <- scan.bin.mr.NEW

themelt.brp.mr$pop <- 'BRP'
themelt.elr.mr$pop <- 'ELR'
themelt.nbh.mr$pop <- 'NBH'
themelt.new.mr$pop <- 'NEW'

save.image(file.path(mpath,'08_phys_plots_pos.rsave'))
################################################

## cross.BRP <- read.cross(format = "csv", dir = mpath, file = 'BRP.mapped.tsp.csv', genotypes=c("1","2","3"), estimate.map = FALSE)
#cross.ELR <- read.cross(format = "csv", dir = mpath, file = 'ELR.mapped.tsp.csv', genotypes=c("1","2","3"), estimate.map = FALSE)
#cross.NBH <- read.cross(format = "csv", dir = mpath, file = 'NBH.mapped.tsp.csv', genotypes=c("1","2","3"), estimate.map = FALSE)
##cross.NEW <- read.cross(format = "csv", dir = mpath, file = 'NEW.mapped.tsp.csv', genotypes=c("1","2","3"), estimate.map = FALSE)

################################################
#cor_nbh <- get_cor(cross_NBH)
#cor_elr <- get_cor(cross_ELR)
#cor_brp <- get_cor(cross_BRP)
#cor_new <- get_cor(cross_NEW)
################################################

#cross_BRP <- flip.order(cross_BRP, names(cor_brp)[which(cor_brp < 0)])
#cross_NBH <- flip.order(cross_NBH, names(cor_nbh)[which(cor_nbh < 0)])
#cross_NEW <- flip.order(cross_NEW, names(cor_new)[which(cor_new < 0)])
#cross_ELR <- flip.order(cross_ELR, names(cor_elr)[which(cor_elr < 0)])


#################################################
### NBH map info
#nbh_gmap <- pull.map(cross_NBH)
#nbh_pmap <- pull.map(cross_NBH)
#
#nbh_names <- lapply(nbh_pmap, function(X) {
#   names(X)
# })
#nbh_pmap <- lapply(nbh_pmap, function(X) {
#   return(as.numeric(gsub("[0-9]+:", "", names(X))))
# })
#
# for (i in 1:24) {
#   names(nbh_pmap[[i]]) <- names(nbh_gmap[[i]])
# }
#################################################


#remap_df <- function(cross_tomap, cross_usemap){
#
# gmap <- pull.map(cross_tomap)
#
# nms <- lapply(gmap, function(X) {
#    names(X)
#  })
#
# chr <- lapply(gmap, function(X) {
#    return(as.numeric(gsub(":[0-9]+", "", names(X))))
#  })
#
# pmap <- pull.map(cross_tomap)
# pmap <- lapply(pmap, function(X) {
#    return(as.numeric(gsub("[0-9]+:", "", names(X))))
#  })
#
# for (i in 1:24) {
#   names(pmap[[i]]) <- names(gmap[[i]])
# }
#
# nms <- unlist(nms,use.names = F)
# pos <- unlist(pmap,use.names = F)
# chr <- unlist(chr,use.names = F)
# data.frame(chr=chr, pos=pos,row.names=nms)
#}
#NEW_info <- remap_df(reorder_NEW, NBH_NEW)
#BRP_info <- remap_df(reorder_BRP, NBH_BRP)
