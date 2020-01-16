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

load(file.path(mpath,'08_phys_plots_pos.rsave'))


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
####################################################################################

####################################################################################
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
####################################################################################

################################################
### ggplot popgen locations
nbh.popgen <- read.table(file.path(mpath,"outliersNBH.txt.ncbi.lifted"), sep = "\t", header = T)
new.popgen <- read.table(file.path(mpath,"outliersNYC.txt.ncbi.lifted"), sep = "\t", header = T)
elr.popgen <- read.table(file.path(mpath,"outliersER.txt.ncbi.lifted"), sep = "\t", header = T)
brp.popgen <- read.table(file.path(mpath,"outliersBP.txt.ncbi.lifted"), sep = "\t", header = T)
################################################

################################################
### Use nbh coords but elr and new popgen
#new.rank <- cnv.popgen(cross.nbh, new.popgen, top = 50)
nbh.rank <- cnv.popgen(cross_NBH, nbh.popgen, top = 50)
elr.rank <- cnv.popgen(cross_ELR, elr.popgen, top = 50)
#brp.rank <- cnv.popgen(cross.nbh, brp.popgen, top = 50)
################################################

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


## get_genes_cm(chr=1, start = 20,stop = 30,models = nbh_gene_models, colm = 'cm_mid')
######## Plot phys pos x marker order ##########################################

png("/home/jmiller1/public_html/ELR_NBH_physpo_filt.png", width=1500, height=1500)
par(mfrow=c(4,6))

for (i in 1:24){

 Y <- c(0, as.numeric(gsub(".*:","",markernames(cross_NBH,i))))
 X <- 1:length(Y)

 A <- c(0, as.numeric(gsub(".*:","",markernames(cross_ELR,i))))
 B <- 1:length(A)

 ymax <- max(c(max(Y),max(A))
 xmax <- max(c(length(Y),length(A))

 plot(c(0,xmax),c(0,ymax), type="n", xlab=paste('chr',i), ylab='physical position')
 points(X,Y, col='blue')
 points(A,B, col='yellow')

 }
dev.off()
################################################################################

######## Plot phys pos x marker order ##########################################

png("/home/jmiller1/public_html/NBH_physpo_filt.png", width=1500, height=1500)
par(mfrow=c(4,6))

for (i in 1:24){

 Y <- c(0, as.numeric(gsub(".*:","",markernames(cross,i))))
 X <- 1:length(Y)

 plot(X,Y, xlab=paste('chr',i), ylab='physical position',cex.axis=2,pch=21,cex.main=2)


 }
dev.off()
################################################################################

#### SEG DISTORTION X PHYSICAL POSITION
gt <- geno.table(cross_nbh)

png("/home/jmiller1/public_html/NBH_segdist.png", width=1500, height=1500)
par(mfrow=c(4,6))

for (i in 1:24){

 ind <- which(gt$chr == i)
 X <- c(as.numeric(gsub(".*:","",rownames(gt)[ind])))
 Y <- c(-log10(gt[ind,'P.value']))

 plot(X,Y,ylim=c(0,4),xlim=c(0,max(X)), xlab=paste('chr',i), ylab='-log10 pval', cex.axis=2,pch=21,cex.main=2)


 }
dev.off()
################################################################################


#### ELR and NBH seg dist
png("/home/jmiller1/public_html/NBH_ELR_segdist.png", width=1750, height=1750)
par(mfrow=c(4,6))

for (i in 1:24){


 ind.nbh <- which(gt.nbh$chr == i)
 X <- c(as.numeric(gsub(".*:","",rownames(gt.nbh)[ind.nbh])))
 Y <- c(-log10(gt.nbh[ind.nbh,'P.value']))


 ind.elr <- which(gt.elr$chr == i)
 A <- c(as.numeric(gsub(".*:","",rownames(gt.elr)[ind.elr])))
 B <- c(-log10(gt.elr[ind.elr,'P.value']))

 ymax <- max(c(max(Y),max(B)))
 xmax <- max(c(max(X),max(A)))

 plot(c(0,xmax),c(0,ymax),ylim=c(0,4), type="n", xlab=paste('chr',i,'phys position'), ylab='-log10 p.value',cex.lab=2,cex.sub=2, cex.axis=2,cex.main=2)
 points(X,Y, col='blue',pch=19)
 points(A,B, col='red',pch=19)

 }
dev.off()


################################################################################
cross_NBH$pheno$pheno_norm <- round(nqrank(cross_NBH$pheno$Pheno),5)
cross_ELR$pheno$pheno_norm <- round(nqrank(cross_ELR$pheno$Pheno),5)

cross_NBH <- sim.geno(cross_NBH, step=1, error.prob=0.025, off.end=5, map.function="kosambi", n.draws=160)
cross_ELR <- sim.geno(cross_ELR, step=1, error.prob=0.025, off.end=5, map.function="kosambi", n.draws=160)
cross_NBH <- calc.genoprob(cross_NBH, step=1, error.prob=0.025, off.end=5, map.function="kosambi")
cross_ELR <- calc.genoprob(cross_ELR, step=1, error.prob=0.025, off.end=5, map.function="kosambi")
scan_nbh <- scanone(cross_NBH, method = "mr", model = "binary", pheno.col = 4)
scan_elr <- scanone(cross_ELR, method = "mr", model = "binary", pheno.col = 4)
scan_nbh <- scanone(cross_NBH, method = "imp", model = "binary", pheno.col = 4)
scan_elr <- scanone(cross_ELR, method = "imp", model = "binary", pheno.col = 4)



png("/home/jmiller1/public_html/NBH_ELR_segdist_binary.png", width=1750, height=1750)
par(mfrow=c(4,6))

for (i in 1:24){

 ind.nbh <- which(gt.nbh$chr == i)
 Ylod <- scan_nbh$lod[ind.nbh]
 X <- c(as.numeric(gsub(".*:","",rownames(gt.nbh)[ind.nbh])))
 Y <- c(-log10(gt.nbh[ind.nbh,'P.value']))

 ind.elr <- which(gt.elr$chr == i)
 Blod <- scan_elr$lod[ind.elr]
 A <- c(as.numeric(gsub(".*:","",rownames(gt.elr)[ind.elr])))
 B <- c(-log10(gt.elr[ind.elr,'P.value']))

 ymax <- max(c(max(Y),max(B)))
 xmax <- max(c(max(X),max(A)))


 if(max(Ylod) > 5 ) Ylod <- rescale(scan_nbh$lod[ind.nbh],c(0,5))
 if(max(Blod) > 5 ) Blod <- rescale(scan_elr$lod[ind.elr],c(0,5))

 plot(c(0,xmax),c(0,ymax),ylim=c(0,5), type="n", xlab=paste('chr',i,'phys position'), ylab='-log10 p.value',cex.lab=2,cex.sub=2, cex.axis=2,cex.main=2)
 points(X,Ylod, col='cornflowerblue',pch=21)
 points(A,Blod, col='pink',pch=21)
 points(X,Y, col='blue',pch=19)
 points(A,B, col='red',pch=19)

 }
dev.off()





################################################################################
## Effect plot (uses sim.geno)
cross_NBH <- sim.geno(cross_NBH, step=1, error.prob=0.05, off.end=5, map.function="kosambi", n.draws=160)
cross_ELR <- sim.geno(cross_ELR, step=1, error.prob=0.05, off.end=5, map.function="kosambi", n.draws=160)


png("/home/jmiller1/public_html/NBH_effect_scan.png", width=1500, height=500)
par(mfrow=c(2,1))
effectscan(cross_NBH, pheno.col=5, chr=c(1:24), get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-2,2), main= 'normal')
effectscan(cross_NBH, pheno.col=4, chr=c(1:24), get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'binary')
dev.off()


png("/home/jmiller1/public_html/ELR_effect_scan.png", width=1500, height=500)
par(mfrow=c(2,1))
effectscan(cross_ELR, pheno.col=5, chr=c(1:24), get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-2,2), main= 'normal')
effectscan(cross_ELR, pheno.col=4, chr=c(1:24), get.se=T, draw=TRUE, gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5), main= 'binary')
dev.off()

################################################################################


################################################################################
################################################################################
col <- c("slateblue", "violetred", "green3")

nbh <- convert2cross2(cross_NBH)
nbh_map <- insert_pseudomarkers(nbh$gmap, step=1)
nbh_pr <- calc_genoprob(nbh, nbh_map, error_prob=0.025, cores=4)

elr <- convert2cross2(cross_ELR)
elr_map <- insert_pseudomarkers(elr$gmap, step=1)
elr_pr <- calc_genoprob(elr, elr_map, error_prob=0.025, cores=4)

cands <- c("AHR1","aip","ARNT","ARNT2","ahrr","ahr1b","AHR2b")

ahr_nbh <- nbh.gens[which(nbh.gens$gene %in% cands),]
ahr_elr <- nbh.gens[which(elr.gens$gene %in% cands),]

plot_ef <- function(crs,map,pr,ahr,popgen){

 for (chr in 1:24){

  c2eff <- scan1coef(pr[,as.character(chr)], crs$pheno[,"pheno_norm"])

  plot(c2eff, map[as.character(chr)], columns=1:3, col=col, ylim=c(0,5), cex.axis = 2)

    if(any( chr %in% ahr$chr )) {
      indx <- which(ahr$chr %in% chr)
      abline(v=as.numeric(ahr[indx,'pos1']), col='red')
    }

    if(any( chr %in% popgen$chr )) {
      indx <- which(popgen$chr %in% chr)
      abline(v=as.numeric(popgen[indx,'pos1']), col='red')
    }


  last_coef <- unclass(c2eff)[nrow(c2eff),] # pull out last coefficients

  for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
  }


 dev.off()
}

png("/home/jmiller1/public_html/NBH_effectplot.png", width=1500, height=1000)
par(mfrow=c(4,6))
plot_ef(crs = nbh, map = nbh_map, pr = nbh_pr, ahr = ahr_nbh, popgen = nbh.rank)

png("/home/jmiller1/public_html/ELR_effectplot.png", width=1500, height=1000)
par(mfrow=c(4,6))
plot_ef(crs = elr, map = elr_map, pr = elr_pr , ahr = ahr_elr, popgen = elr.rank)
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



png("/home/jmiller1/public_html/lodxdist.png", width=500, height=500)
plot(-log10(nbh_seg),nbh_c2eff[,'AA'], col = 'blue',pch=19, ylim=c(0,1))
points(-log10(nbh_seg),nbh_c2eff[,'BB'], col = 'blue',pch=19)
points(-log10(nbh_seg),nbh_c2eff[,'AB'], col='yellow',pch=19)

dev.off()

png("/home/jmiller1/public_html/lodxdist_nbh_AABB.png", width=500, height=500)
plot(-log10(nbh_seg),nbh_c2eff[,'AA'], col = 'blue',pch=19, ylim=c(0,1))
points(-log10(nbh_seg),nbh_c2eff[,'BB'], col = 'red',pch=19)
##points(-log10(nbh_seg),nbh_c2eff[,'AB'], col='yellow',pch=19)
dev.off()


png("/home/jmiller1/public_html/lodxdist_elr_AABB.png", width=500, height=500)
plot(-log10(elr_seg),elr_c2eff[,'AA'], col = 'blue',pch=19, ylim=c(0,1))
points(-log10(elr_seg),elr_c2eff[,'BB'], col = 'red',pch=19)
##points(-log10(elr_seg),nbh_c2eff[,'AB'], col='yellow',pch=19)
dev.off()
