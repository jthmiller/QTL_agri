#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]

## USE MARKER REGRESSION TO COMPARE ALL LOCI ON BRP AND NEW

pop <- 'ELR'
library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################
load(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))
################################################################################

#cross <- switchAlleles(cross, markernames(cross,c(2,13)))

################################################################################

markers <- markernames(cross)
bin.em <- scanone(cross, pheno.col=4, model='binary', method = "em")

sens <- as.character(cross$pheno$ID[which(cross$pheno$ID %in% pheno_ind(cross,1))])
tol <- as.character(cross$pheno$ID[which(cross$pheno$ID %in% pheno_ind(cross,0))])

gts <- pull.geno(cross)[,markers]
rownames(gts) <- as.character(cross$pheno$ID)

################################################################################

## genotypes for plotting
gt_sens <- lapply(markers, function(mark){
 factor(gts[sens,mark],levels=c(NA,1,2,3),exclude = NULL)
 }
)
names(gt_sens) <- markers

gt_tol <- lapply(markers, function(mark){
 factor(gts[tol,mark],levels=c(NA,1,2,3),exclude = NULL)
 }
)
names(gt_tol) <- markers

################################################################################

## list of proportions for plotting
prop_total <- mapply(function(sens,tol){ return(table(sens) + table(tol)) }, gt_sens, gt_tol, SIMPLIFY = F, USE.NAMES = FALSE)
names(prop_total) <- markers

prop_sens <- mapply( function(prop,total){ return(table(prop) / total) }, gt_sens, prop_total, SIMPLIFY = F,USE.NAMES = FALSE)
names(prop_sens) <- markers

#### NAMES FOR PLOTS
pop_chr <- as.list(paste('CHR',gsub(":.*","",names(gt_sens)),'QTL genotype'))
names(pop_chr) <- markers

#loc_name <- as.list(c('AIP (CHR 2)','ARNT (CHR 8)','ARNT (CHR 13)','AHRb (CHR 18)','QTLa (CHR 24)','QTLb (CHR 24)'))
#names(loc_name) <- markers

chr <- names(gt_sens)
chr <- ifelse(duplicated(chr),paste0(chr,'.',sum(duplicated(chr))),chr)
names(chr) <- markers

################################################################################

single <- function(mark, pop = 'ELR'){

 ydir <- prop_sens[[mark]]
 ytot <- prop_total[[mark]]
 chromo <- pop_chr[[mark]]
 mainL <- pop_chr[[mark]]
 chrL <- pop_chr[[mark]]

 pdfL <- paste0("/home/jmiller1/public_html/", pop, chr[mark],".pdf")

 pdf(pdfL, width=3.5)

 cex_single <- c(.25,.5,.25) * nind(cross)
 xdir <- c(1,2,3)

 plot(c(0.65,3.35), c(-0.1,1.1),
  xaxs="i", xaxt="n", xlab="",
  yaxs="i", yaxt="n", ylab="",
  type="n", main=mainL,
  cex.lab=1.5, cex.main=2)

  lines(xdir, ydir[as.character(c(1,2,3))],col='black',lwd=5)

  points(xdir, ydir[as.character(c(1,2,3))],col=c('black','darkblue','cornflowerblue'), pch=21,bg=c('black','darkblue','cornflowerblue'),
   cex= (12 * (ytot[as.character(c(1,2,3))] / cex_single)))

  text(xdir, ydir[as.character(c(1,2,3))], labels=ytot[as.character(c(1,2,3))],col='white',font=2, cex=2)

  mtext(side=1, line=3, chrL, col="black", font=2,cex=1.5)
  mtext(side=2, line=3, "Proportion Deformed", col="black", font=2, cex=1.5)

  axis(side=1, at=c(1,2,3), labels=c('AA','AB','BB'),font=2, cex.axis=2)
  axis(side=2, at=c(0,0.5,1.0), labels=c('0','0.5','1.0'),font=2, cex.lab=2, cex.axis=2)
dev.off()
}

## sapply(amarker,single,pop = 'ELR')
################################################################################

intxs.bin <- function(loc_a, loc_b, popchr, locbN, main){

 ## interactions
 AA <- table(factor(gts[names(which(gts[sens ,loc_a] == 1)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 AB <- table(factor(gts[names(which(gts[sens ,loc_a] == 2)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 BB <- table(factor(gts[names(which(gts[sens ,loc_a] == 3)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 sensit <- rbind(AA, AB, BB)
 colnames(sensit) <- c('NA','AA','AB','BB')

 AA <- table(factor(gts[names(which(gts[tol  ,loc_a] == 1)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 AB <- table(factor(gts[names(which(gts[tol  ,loc_a] == 2)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 BB <- table(factor(gts[names(which(gts[tol  ,loc_a] == 3)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 resist <- rbind(AA,AB,BB)
 colnames(resist) <- c('NA','AA','AB','BB')

 AA <- table(factor(gts[names(which(gts[c(sens,tol) ,loc_a] == 1)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 AB <- table(factor(gts[names(which(gts[c(sens,tol) ,loc_a] == 2)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 BB <- table(factor(gts[names(which(gts[c(sens,tol) ,loc_a] == 3)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
 total <- rbind(AA, AB, BB)
 colnames(total) <- c('NA','AA','AB','BB')
 rownames(total) <- rownames(resist)

 ## Plot
 xdir <- c(1,2,3)
 ydir <- sensit/total

cexs_hom <- c(0.25^2,0.5^2,0.25^2)
cexs_hom <- cexs_hom * nind(cross)
cexs_het <- c(0.25*0.5,0.5^2,0.25*0.5)
cexs_het <- cexs_het * nind(cross)

pdf(paste0("/home/jmiller1/public_html/",popchr,".pdf"), width=10)
 plot(c(0.65,3.35), c(-0.1,1.1),
  xaxs="i", xaxt="n", xlab="",
  yaxs="i", yaxt="n", ylab="",
  type="n", main=main,
  cex.lab=1.5, cex.main=2)

  rect(1.5, -0.1, 2.5, 1.1,col='lightgrey',border = 'transparent')

  lines(xdir+0.28, ydir[rownames(total),'BB'],col='cornflowerblue',lwd=5)
  lines(xdir, ydir[rownames(total),'AB'],col='darkblue',lwd=5)
  lines(xdir-0.28, ydir[rownames(total),'AA'],col='black',lwd=5)

  points(xdir, ydir[rownames(total),'AB'],col='darkblue', pch=21, bg='darkblue',
   cex= (12 * (total[rownames(total),'AB'] / cexs_het))+2)
  points(xdir+0.28, ydir[rownames(total),'BB'], col='cornflowerblue', pch=21, bg='cornflowerblue',
   cex= (12 * (total[rownames(total),'BB'] / cexs_hom))+2)
  points(xdir-0.28, ydir[rownames(total),'AA'],col='black', pch=21,bg='black',
   cex= (12 * (total[rownames(total),'AA'] / cexs_hom))+2)


  text(xdir+0.28, ydir[rownames(total),'BB'], labels=total[rownames(total),'BB'],col='white',font=2, cex=2)
  text(xdir, ydir[rownames(total),'AB'], labels=total[rownames(total),'AB'],col='white',font=2, cex=2)
  text(xdir-0.28, ydir[rownames(total),'AA'], labels=total[rownames(total),'AA'],col='white',font=2, cex=2)

  mtext(side=1, line=3, locbN , col="black", font=2,cex=1.5)
  mtext(side=2, line=3, "Proportion Deformed", col="black", font=2, cex=1.5)

  axis(side=1, at=c(1,2,3), labels=c('AA','AB','BB'),font=2, cex.axis=2)
  axis(side=2, at=c(0,0.5,1.0), labels=c('0','0.5','1.0'),font=2, cex.lab=2, cex.axis=2)
dev.off()
}
################################################################################

elr_qtl18 <- find.marker(cross,18,74.7)
elr_qtl15 <- find.marker(cross,15,8.82)

intxs.bin(elr_qtl18,elr_qtl15,  popchr = "18v15", locbN = 'AHR Genotype', main = 'AHR x CYP')
single(elr_qtl18, pop = 'ELR')
single(elr_qtl15, pop = 'ELR')
