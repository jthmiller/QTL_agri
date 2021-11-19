#!/bin/R
### first run combine pops for multi-pop cross objects
pop <- 'ELR'

pop <- 'NBH'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
library("ggridges")
library("plyr")
library("scales")
library("ggrepel")
library('qtl')
library('RColorBrewer')
################################################################################
mpath <- '/home/jmiller1/QTL_agri/data'
setwd(mpath)
#load(file.path(mpath,'08_phys_plots_pos.rsave'))
####################################################################################


################################################################################
## FUNCTIONS ###
################################################################################
################################################################################

plot_ef <- function(crs,map,pr,ahr,popgen,chs,main,qtl,model=c("bin","pheno_norm"),...){

 for (chr in chs){

  c2eff <- scan1coef(pr[,as.character(chr)], crs$pheno[,model])

  plot(c2eff, map[as.character(chr)], columns=1:3, col=col, ylim=c(0,1), cex.axis = 2,main=main,...)

    if(any( chr %in% ahr$chr )) {
      indx <- which(ahr$chr %in% chr)
      abline(v=as.numeric(ahr[indx,'pos']), col='red',lwd=0.5)
      #xleft, ybottom, xright, ytop,
    }
    #if(any( chr %in% qtl$chr )) {
    #  indx <- which(qtl$chr %in% chr)
    #  abline(v=as.numeric(qtl[indx,'pos']), col='green')
    #}


  last_coef <- unclass(c2eff)[nrow(c2eff),] # pull out last coefficients

  for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
  }

}


################################################################################
################################################################################

plot_pgen <- function(crs, chrs, stat, map, ahr, ahr_clm, colnm, popgen, ylimo, rank_clm, stat_name,...){

 for (chr in chrs){

  xl <- summary(pull.map(crs))[chr,'length']
  ind <- which(stat$chr == chr)

  Y <- stat[ind,colnm]
  X <- stat[ind,map]/1000000
##  plot(X, Y, col='blue', cex.axis = 2, ylim = ylimo, xlim = c(0,xl), main=paste('CHR',chr), cex.main=2)

  plot(X, Y, col='black',type="n",xlim=c(0,max(X)), ylim = ylimo, main=NULL,
   xlab='physical position', ylab=stat_name, xaxs="i",yaxs="i", mgp = c(1, 0.5, 0),...)

    if(any( chr %in% ahr$chr )) {
      indx <- which(ahr$chr %in% chr)
      #rect(ahr[indx,ahr_clm]/1000000,ylimo[1],ahr[indx,'stp']/1000000,ylimo[2],lwd=0.5,col=alpha('lightgrey',.5))
      abline(v=as.numeric(ahr[indx,ahr_clm])/1000000,
       col='red',lwd=0.5)
    }

    if(any( chr %in% popgen$chr )) {
      indx <- which(popgen$chr %in% chr)
      rect(popgen[indx,'start']/1000000,ylimo[1],popgen[indx,'end']/1000000,ylimo[2],
       border = NA,lwd=0,col=alpha('lightgrey',.5))

      #abline(v=as.numeric(popgen[indx,rank_clm])/1000000, col='grey',lwd=2)
    }
  points(X, Y, col='black',...)
 }
}
################################################################################
################################################################################
################################################################################
################################################################################
## END FUNCTIONS
################################################################################
################################################################################
################################################################################

################################################################################
pop <- 'NBH'
erprob <- 0.001
mapfile <- paste0(pop,'_reorder_imp_nopar')
filename <- file.path(mpath,mapfile)
cross_NBH <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
newmap <- est.map(cross_NBH, error.prob=erprob, map.function="kosambi", maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)
cross_NBH <- replace.map(cross_NBH, newmap)
summary(pull.map(cross_NBH))
################################################################################

################################################################################
pop <- 'ELR'
erprob <- 0.001
mapfile <- paste0(pop,'_reorder_imp_mapped')
filename <- file.path(mpath,mapfile)
cross_ELR <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
newmap <- est.map(cross_ELR, error.prob=erprob, map.function="kosambi", maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE)
cross_ELR <- replace.map(cross_ELR, newmap)
summary(pull.map(cross_NBH))
################################################################################

################################################################################
pop <- 'NBH'
################################################################################

names(cross_NBH$geno) <- ifelse(names(cross_NBH$geno) == "X","5",names(cross_NBH$geno))
attr(cross_NBH$geno[["5"]], 'class') <- 'A'
names(cross_ELR$geno) <- ifelse(names(cross_ELR$geno) == "X","5",names(cross_ELR$geno))
attr(cross_ELR$geno[["5"]], 'class') <- 'A'

################################################################################

scan_nbh <- scanone(cross_NBH, method = "mr", model = "binary", pheno.col = 4)
scan_elr <- scanone(cross_ELR, method = "mr", model = "binary", pheno.col = 4)
#scan_nbh <- scanone(cross_NBH, method = "em", model = "binary", pheno.col = 4)
#scan_elr <- scanone(cross_ELR, method = "em", model = "binary", pheno.col = 4)

################################################################################
################################################################################

#nbh_marks <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross_NBH)[[X]], 1)} ))
#dwnsmpl_NBH <- pull.markers(cross_NBH, nbh_marks)
#elr_marks <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross_ELR)[[X]], 1)} ))
#dwnsmpl_ELR <- pull.markers(cross_ELR, elr_marks)

################################################################################
################################################################################
###### plot_ef rQTL2

col <- c("slateblue", "violetred", "green3")

nbh <- convert2cross2(cross_NBH)
nbh_map <- insert_pseudomarkers(nbh$gmap, step=1)
nbh_pr <- calc_genoprob(nbh, nbh_map, error_prob=0.001, cores=4)

elr <- convert2cross2(cross_ELR)
elr_map <- insert_pseudomarkers(elr$gmap, step=1)
elr_pr <- calc_genoprob(elr, elr_map, error_prob=0.001, cores=4)

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
cands <- c("AHR1","aip","ARNT","ARNT2","ahrr","ahr1b","AHR2b")
# add arnts (forgot to scan for them)

################################################

nbh_gens <- cnv.ahrs(cross_NBH, AHRdf = AHR.bed, EXP = F)
elr_gens <- cnv.ahrs(cross_ELR, AHRdf = AHR.bed, EXP = F)
ahr_nbh <- nbh_gens[which(nbh_gens$gene %in% cands),]
ahr_elr <- elr_gens[which(elr_gens$gene %in% cands),]

################################################

#nbh_gens_dwnsmpl <- cnv.ahrs(dwnsmpl_NBH, AHRdf = AHR.bed, EXP = F)
#elr_gens_dwnsmpl <- cnv.ahrs(dwnsmpl_ELR, AHRdf = AHR.bed, EXP = F)
#ahr_nbh_dwnsmpl <- nbh_gens[which(nbh_gens$gene %in% cands),]
#ahr_elr_dwnsmpl <- elr_gens[which(elr_gens$gene %in% cands),]

################################################################################
### ggplot popgen locations

nbh.popgen <- read.table(file.path(mpath,"outliersNBH.txt.ncbi.lifted"), sep = "\t", header = T)
new.popgen <- read.table(file.path(mpath,"outliersNYC.txt.ncbi.lifted"), sep = "\t", header = T)
elr.popgen <- read.table(file.path(mpath,"outliersER.txt.ncbi.lifted"), sep = "\t", header = T)
brp.popgen <- read.table(file.path(mpath,"outliersBP.txt.ncbi.lifted"), sep = "\t", header = T)

################################################################################
### Use nbh coords but elr and new popgen
#new.rank <- cnv.popgen(cross.nbh, new.popgen, top = 50)
#nbh.rank <- cnv.popgen(cross_NBH, nbh.popgen, top = 62)
nbh.rank <- cnv.popgen(cross_NBH, nbh.popgen, top = 121)
# dim(elr.popgen[which(elr.popgen$chrom %in% paste0('chr',1:24) & elr.popgen$rank < 100),])
# top 100 have 50 mapped ranked outliers
elr.rank <- cnv.popgen(cross_ELR, elr.popgen, top = 200)
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
################################################################################

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

#pbs$nbh_cm_dwnsmpl <- conv_popstat(dwnsmpl_NBH, popgen=pbs, whichcol='mid', newname='nbh_cm_dwnsmpl')$nbh_cm_dwnsmpl
#pbs$elr_cm_dwnsmpl <- conv_popstat(dwnsmpl_ELR, popgen=pbs, whichcol='mid', newname='elr_cm_dwnsmpl')$elr_cm_dwnsmpl
################################################################################

pfst <- file.path(mpath, 'pfst.txt.ncbi.lifted')
pfst <- read.table(pfst, sep = "\t", header = T)
pfst$mid <- pfst$start + (abs(pfst$end - pfst$start) * .5)
pfst$Scaffold <- gsub('chr',"",pfst$Scaffold)
colnames(pfst)[1] <- 'chr'
pfst <- pfst[!is.na(as.numeric(pfst$chr)),]
pfst$nbh_cm <- conv_popstat(cross_NBH, popgen=pfst, whichcol='mid',newname='nbh_cm')$nbh_cm
pfst$elr_cm <- conv_popstat(cross_ELR, popgen=pfst, whichcol='mid',newname='elr_cm')$elr_cm

#pfst$nbh_cm_dwnsmpl <- conv_popstat(dwnsmpl_NBH, popgen=pfst, whichcol='mid',newname='nbh_cm_dwnsmpl')$nbh_cm_dwnsmpl
#pfst$elr_cm_dwnsmpl <- conv_popstat(dwnsmpl_ELR, popgen=pfst, whichcol='mid',newname='elr_cm_dwnsmpl')$elr_cm_dwnsmpl
################################################################################

taj <- file.path(mpath, 'tajstat.txt.ncbi.lifted')
taj <- read.table(taj, sep = "\t", header = T)
taj$mid <- taj$start + (abs(taj$end - taj$start) * .5)
taj$Scaffold <- gsub('chr',"",taj$Scaffold)
colnames(taj)[1] <- 'chr'
taj <- taj[!is.na(as.numeric(taj$chr)),]
taj$nbh_cm <- conv_popstat(cross_NBH, popgen=taj, whichcol='mid',newname='nbh_cm')$nbh_cm
taj$elr_cm <- conv_popstat(cross_ELR, popgen=taj, whichcol='mid',newname='elr_cm')$elr_cm

#taj$nbh_cm_dwnsmpl <- conv_popstat(dwnsmpl_NBH, popgen=taj, whichcol='mid',newname='nbh_cm_dwnsmpl')$nbh_cm_dwnsmpl
#taj$elr_cm_dwnsmpl <- conv_popstat(dwnsmpl_ELR, popgen=taj, whichcol='mid',newname='elr_cm_dwnsmpl')$elr_cm_dwnsmpl
################################################################################

pi <- file.path(mpath, 'piper.txt.ncbi.lifted')
pi <- read.table(pi, sep = "\t", header = T)
pi$mid <- pi$start + (abs(pi$end - pi$start) * .5)
pi$Scaffold <- gsub('chr',"",pi$Scaffold)
colnames(pi)[1] <- 'chr'
pi <- pi[!is.na(as.numeric(pi$chr)),]

pi$F.NBH <- pi$NBH - pi$F
pi$BI.NBH <- pi$NBH - pi$BI
pi$ER.KC <- pi$ER - pi$KC
pi$ER.SH <- pi$ER - pi$SH

pi$nbh_cm <- conv_popstat(cross_NBH, popgen=pi, whichcol='mid',newname='nbh_mp')$nbh_cm
pi$elr_cm <- conv_popstat(cross_ELR, popgen=pi, whichcol='mid',newname='elr_mp')$elr_cm

#pi$nbh_cm_dwnsmpl <- conv_popstat(dwnsmpl_NBH, popgen=pi, whichcol='mid',newname='nbh_mp')$nbh_cm_dwnsmpl
#pi$elr_cm_dwnsmpl <- conv_popstat(dwnsmpl_ELR, popgen=pi, whichcol='mid',newname='elr_mp')$elr_cm_dwnsmpl
################################################################################
################################################################################

################################################################################
################################################################################
################################################################################
## Correlate lod and segregation distortion

elr_c2eff <- lapply(c(1:4,6:24),function(X) { scan1coef(elr_pr[ ,as.character(X)], elr$pheno[,"bin"]) })
elr_c2eff <- do.call(rbind,elr_c2eff)
elr_seg <- geno.table(cross_ELR)[rownames(elr_c2eff),'P.value']

nbh_c2eff <- lapply(1:24,function(X) {
 scan1coef(nbh_pr[,as.character(X)], nbh$pheno[,"bin"])
})
nbh_c2eff <- do.call(rbind,nbh_c2eff)
nbh_seg <- geno.table(cross_NBH)[rownames(nbh_c2eff),'P.value']

################################################################################
save.image(file.path(mpath,'supplemental_plot_env.rsave'))
################################################################################
