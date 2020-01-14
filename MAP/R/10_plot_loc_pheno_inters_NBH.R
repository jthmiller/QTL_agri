#!/bin/R

##pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]

## USE MARKER REGRESSION TO COMPARE ALL LOCI ON BRP AND NEW
## pop <- 'NBH'
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]
library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################
################################################################################
## Read cross
cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)


loc_2 <- find.marker(cross, full.norm.imp$chr[1],full.norm.imp$pos[1])
loc_8 <- find.marker(cross, full.norm.imp$chr[2],full.norm.imp$pos[2])
loc_13 <- find.marker(cross, full.norm.imp$chr[3],full.norm.imp$pos[3])
loc_18 <- find.marker(cross, full.norm.imp$chr[4],full.norm.imp$pos[4])
loc_24a <- find.marker(cross, full.norm.imp$chr[5],full.norm.imp$pos[5])
loc_24b <- find.marker(cross, full.norm.imp$chr[6],full.norm.imp$pos[6])


################################################################################
### NBH ####
nbh_cross <- cross
nbh_sens <- as.character(nbh_cross$pheno$ID[which(nbh_cross$pheno$ID %in% pheno_ind(nbh_cross,1))])
nbh_tol <- as.character(nbh_cross$pheno$ID[which(nbh_cross$pheno$ID %in% pheno_ind(nbh_cross,0))])

nbh_gts <- pull.geno(nbh_cross)[,c(loc_2,loc_8,loc_13,loc_18,loc_24a,loc_24b)]
rownames(nbh_gts) <- as.character(nbh_cross$pheno$ID)

################################################################################
## single locus
nbh_sens_2 <- factor(nbh_gts[nbh_sens,c(loc_2)],levels=c(NA,1,2,3),exclude = NULL)
nbh_sens_8 <- factor(nbh_gts[nbh_sens,c(loc_8)],levels=c(NA,1,2,3),exclude = NULL)
nbh_sens_13 <- factor(nbh_gts[nbh_sens,c(loc_13)],levels=c(NA,1,2,3),exclude = NULL)
nbh_sens_18 <- factor(nbh_gts[nbh_sens,c(loc_18)],levels=c(NA,1,2,3),exclude = NULL)
nbh_sens_24 <- factor(nbh_gts[nbh_sens,c(loc_24a)],levels=c(NA,1,2,3),exclude = NULL)

nbh_tol_2 <- factor(nbh_gts[nbh_tol,c(loc_2)],levels=c(NA,1,2,3),exclude = NULL)
nbh_tol_8 <- factor(nbh_gts[nbh_tol,c(loc_8)],levels=c(NA,1,2,3),exclude = NULL)
nbh_tol_13 <- factor(nbh_gts[nbh_tol,c(loc_13)],levels=c(NA,1,2,3),exclude = NULL)
nbh_tol_18 <- factor(nbh_gts[nbh_tol,c(loc_18)],levels=c(NA,1,2,3),exclude = NULL)
nbh_tol_24 <- factor(nbh_gts[nbh_tol,c(loc_24a)],levels=c(NA,1,2,3),exclude = NULL)

nbh_2 <- table(nbh_sens_2) + table(nbh_tol_2)
nbh_8 <- table(nbh_sens_8) + table(nbh_tol_8)
nbh_13 <- table(nbh_sens_13) + table(nbh_tol_13)
nbh_18 <- table(nbh_sens_18) + table(nbh_tol_18)
nbh_24 <- table(nbh_sens_24) + table(nbh_tol_24)

nbh_2_total <- table(nbh_sens_2) + table(nbh_tol_2)
nbh_8_total <- table(nbh_sens_8) + table(nbh_tol_8)
nbh_13_total <- table(nbh_sens_13) + table(nbh_tol_13)
nbh_18_total <- table(nbh_sens_18) + table(nbh_tol_18)
nbh_24_total <- table(nbh_sens_24) + table(nbh_tol_24)

nbh_tol_2 <- table(nbh_sens_2) / nbh_2_total
nbh_tol_8 <- table(nbh_sens_8) / nbh_8_total
nbh_tol_13 <- table(nbh_sens_13) / nbh_13_total
nbh_tol_18 <- table(nbh_sens_18) / nbh_18_total
nbh_tol_24 <- table(nbh_sens_24) / nbh_24_total

################################################################################
## single plot 2
cex_single <- c(.25,.5,.25) * 94
xdir <- c(1,2,3)


single <- function(ydir, pop_chr, mainL, chrL, pdfL){

pdf(paste0("/home/jmiller1/public_html/", pdfL,".pdf"), width=3.5)
 plot(c(0.65,3.35), c(-0.1,1.1),
  xaxs="i", xaxt="n", xlab="",
  yaxs="i", yaxt="n", ylab="",
  type="n", main=mainL,
  cex.lab=1.5, cex.main=2)

  lines(xdir, ydir[as.character(c(1,2,3))],col='black',lwd=5)

  points(xdir, ydir[as.character(c(1,2,3))],col=c('black','darkblue','cornflowerblue'), pch=21,bg=c('black','darkblue','cornflowerblue'),
   cex= (12 * (nbh[as.character(c(1,2,3))] / cex_single)))

  text(xdir, ydir[as.character(c(1,2,3))], labels=nbh[as.character(c(1,2,3))],col='white',font=2, cex=2)

  mtext(side=1, line=3, chrL, col="black", font=2,cex=1.5)
  mtext(side=2, line=3, "Proportion Deformed", col="black", font=2, cex=1.5)

  axis(side=1, at=c(1,2,3), labels=c('AA','AB','BB'),font=2, cex.axis=2)
  axis(side=2, at=c(0,0.5,1.0), labels=c('0','0.5','1.0'),font=2, cex.lab=2, cex.axis=2)
dev.off()
}


single(nbh_tol_2, pop_chr = nbh_2, mainL = 'QTL2', chrL = 'CHR 2', pdfL = "nbh_2")
single(nbh_tol_8, pop_chr = nbh_8, mainL = 'QTL8', chrL = 'CHR 8', pdfL = "nbh_8")
single(nbh_tol_13, pop_chr = nbh_13, mainL = 'QTL13', chrL = 'CHR 13', pdfL = "nbh_13")
single(nbh_tol_18, pop_chr = nbh_18, mainL = 'QTL18', chrL = 'CHR 18', pdfL = "nbh_18")
single(nbh_tol_24, pop_chr = nbh_24, mainL = 'QTL24', chrL = 'CHR 24', pdfL = "nbh_24")

#

################################################################################
## single plot AIP

cex_single <- c(.25,.5,.25) * 94
xdir <- c(1,2,3)
ydir <- nbh_tol_18

pdf(paste0("/home/jmiller1/public_html/",pdfL,".pdf"), width=3.5)
 plot(c(0.65,3.35), c(-0.1,1.1),
  xaxs="i", xaxt="n", xlab="",
  yaxs="i", yaxt="n", ylab="",
  type="n", main='AHR',
  cex.lab=1.5, cex.main=2)

  lines(xdir, ydir[as.character(c(1,2,3))],col='black',lwd=5)

  points(xdir, ydir[as.character(c(1,2,3))],col=c('black','darkblue','cornflowerblue'), pch=21,bg=c('black','darkblue','cornflowerblue'),
   cex= (12 * (nbh_ahr[as.character(c(1,2,3))] / cex_single)))

  text(xdir, ydir[as.character(c(1,2,3))], labels=nbh_ahr[as.character(c(1,2,3))],col='white',font=2, cex=2)

  mtext(side=1, line=3, 'AHR', col="black", font=2,cex=1.5)
  mtext(side=2, line=3, "Proportion Deformed", col="black", font=2, cex=1.5)

  axis(side=1, at=c(1,2,3), labels=c('AA','AB','BB'),font=2, cex.axis=2)
  axis(side=2, at=c(0,0.5,1.0), labels=c('0','0.5','1.0'),font=2, cex.lab=2, cex.axis=2)
dev.off()
}
################################################################################



intxs <- function(loc_a,loc_b,popchr,locbN,main){
 ## interactions
 nbh_AA <- table(factor(nbh_gts[names(which(nbh_gts[nbh_sens ,loc_a] == 1)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 nbh_AB <- table(factor(nbh_gts[names(which(nbh_gts[nbh_sens ,loc_a] == 2)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 nbh_BB <- table(factor(nbh_gts[names(which(nbh_gts[nbh_sens ,loc_a] == 3)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 nbh_sensit <- rbind(nbh_AA, nbh_AB, nbh_BB)
 colnames(nbh_sensit) <- c('NA','AA','AB','BB')

 nbh_AA <- table(factor(nbh_gts[names(which(nbh_gts[nbh_tol  ,loc_a] == 1)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 nbh_AB <- table(factor(nbh_gts[names(which(nbh_gts[nbh_tol  ,loc_a] == 2)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 nbh_BB <- table(factor(nbh_gts[names(which(nbh_gts[nbh_tol  ,loc_a] == 3)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 nbh_resist <- rbind(nbh_AA,nbh_AB,nbh_BB)
 colnames(nbh_resist) <- c('NA','AA','AB','BB')

 nbh_AA <- table(factor(nbh_gts[names(which(nbh_gts[c(nbh_sens,nbh_tol) ,loc_a] == 1)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 nbh_AB <- table(factor(nbh_gts[names(which(nbh_gts[c(nbh_sens,nbh_tol) ,loc_a] == 2)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 nbh_BB <- table(factor(nbh_gts[names(which(nbh_gts[c(nbh_sens,nbh_tol) ,loc_a] == 3)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 nbh_total <- rbind(nbh_AA, nbh_AB, nbh_BB)
 colnames(nbh_total) <- c('NA','AA','AB','BB')
 rownames(nbh_total) <- rownames(nbh_resist)

 ## Plot
 xdir <- c(1,2,3)
 ydir <- nbh_sensit/nbh_total

cexs_hom <- c(0.25^2,0.5^2,0.25^2)
cexs_hom <- cexs_hom * 94
cexs_het <- c(0.25*0.5,0.5^2,0.25*0.5)
cexs_het <- cexs_het * 94

pdf(paste0("/home/jmiller1/public_html/",popchr,".pdf"), width=10)
 plot(c(0.65,3.35), c(-0.1,1.1),
  xaxs="i", xaxt="n", xlab="",
  yaxs="i", yaxt="n", ylab="",
  type="n", main=main,
  cex.lab=1.5, cex.main=2)

  rect(1.5, -0.1, 2.5, 1.1,col='lightgrey',border = 'transparent')

  lines(xdir+0.28, ydir[rownames(nbh_total),'BB'],col='cornflowerblue',lwd=5)
  lines(xdir, ydir[rownames(nbh_total),'AB'],col='darkblue',lwd=5)
  lines(xdir-0.28, ydir[rownames(nbh_total),'AA'],col='black',lwd=5)

  points(xdir, ydir[rownames(nbh_total),'AB'],col='darkblue', pch=21, bg='darkblue',
   cex= (12 * (nbh_total[rownames(nbh_total),'AB'] / cexs_het))+2)
  points(xdir+0.28, ydir[rownames(nbh_total),'BB'], col='cornflowerblue', pch=21, bg='cornflowerblue',
   cex= (12 * (nbh_total[rownames(nbh_total),'BB'] / cexs_hom))+2)
  points(xdir-0.28, ydir[rownames(nbh_total),'AA'],col='black', pch=21,bg='black',
   cex= (12 * (nbh_total[rownames(nbh_total),'AA'] / cexs_hom))+2)


  text(xdir+0.28, ydir[rownames(nbh_total),'BB'], labels=nbh_total[rownames(nbh_total),'BB'],col='white',font=2, cex=2)
  text(xdir, ydir[rownames(nbh_total),'AB'], labels=nbh_total[rownames(nbh_total),'AB'],col='white',font=2, cex=2)
  text(xdir-0.28, ydir[rownames(nbh_total),'AA'], labels=nbh_total[rownames(nbh_total),'AA'],col='white',font=2, cex=2)

  mtext(side=1, line=3, locbN , col="black", font=2,cex=1.5)
  mtext(side=2, line=3, "Proportion Deformed", col="black", font=2, cex=1.5)

  axis(side=1, at=c(1,2,3), labels=c('AA','AB','BB'),font=2, cex.axis=2)
  axis(side=2, at=c(0,0.5,1.0), labels=c('0','0.5','1.0'),font=2, cex.lab=2, cex.axis=2)
dev.off()
}
################################################################################

intxs(loc_2,loc_18, popchr = 'popchr',locbN = 'locbN', main = 'main')
