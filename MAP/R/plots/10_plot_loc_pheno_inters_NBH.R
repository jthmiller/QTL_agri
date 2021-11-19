#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]

## USE MARKER REGRESSION TO COMPARE ALL LOCI ON BRP AND NEW

pop <- 'NBH'
library('qtl')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################
load(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))
################################################################################
##cross <- switchAlleles(cross, markernames(cross,c(13)))
################################################################################
## Read cross

#cross2 <- fill.geno(cross, method="argmax", error.prob=0.001, map.function="kosambi", min.prob=0.95)
#cross2 <- argmax.geno(cross2, step=1, off.end=1, error.prob=0.001, map.function="kosambi", stepwidth="fixed")
#markers <- colnames(pull.argmaxgeno(cross2))

################################################################################

##cross <- fill.geno(cross, method="no_dbl_XO", error.prob=0.002, map.function="kosambi", min.prob=0.95)
##
##
##markers <- mapply(
##  function(crs=cross,X,Y){ find.marker(crs,X,Y) },
##   X = bin.em$chr, Y = bin.em$pos)

##markers <- mapply(
##  function(crs=cross,X,Y){ find.marker(crs,X,Y) },
##   X = full.norm.imp$chr, Y = full.norm.imp$pos)

################################################################################

markers <- markernames(cross)
bin.em <- scanone(cross, pheno.col=4, model='binary', method = "em")
bin.hk <- scanone(cross, pheno.col=4, model='binary', method = "hk")

sens <- as.character(cross$pheno$ID[which(cross$pheno$ID %in% pheno_ind(cross,1))])
tol <- as.character(cross$pheno$ID[which(cross$pheno$ID %in% pheno_ind(cross,0))])

gts <- pull.geno(cross)[,markers]
rownames(gts) <- as.character(cross$pheno$ID)


#################################################################################
#gt_sens <- factor(gts[sens,mark],levels=c(NA,1,2,3),exclude = NULL)
#gt_tol <- factor(gts[tol,mark],levels=c(NA,1,2,3),exclude = NULL)
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

single <- function(mark, pop = 'NBH'){

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

nbh_qtl2 <- find.marker(cross,2, 39.87)
nbh_qtl8 <- find.marker(cross,8, 59.22)
nbh_qtl13 <- find.marker(cross,13,46.602)
nbh_qtl24 <- find.marker(cross,24,10.455051)
nbh_qtl1 <- find.marker(cross,1,60.559873)

geno.crosstab(cross, nbh_qtl24,nbh_qtl1)

intxs.bin(nbh_qtl24, nbh_qtl1, popchr = "1v24", locbN = 'test', main = 'test2')

intxs.bin(nbh_qtl1,nbh_qtl24,  popchr = "1v24", locbN = 'test', main = 'test2')

intxs.bin(nbh_qtl2,nbh_qtl13,  popchr = "2v13", locbN = 'AIP Genotype', main = 'AIP x ARNT')

single('2:27373969', pop = 'NBH')

single(nbh_qtl13, pop = 'NBH')
