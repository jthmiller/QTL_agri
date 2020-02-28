#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]

## USE MARKER REGRESSION TO COMPARE ALL LOCI ON BRP AND NEW

##pop <- 'NBH'
library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'

fl <- file.path(mpath,fl)
load(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))

################################################################################

load(file.path(mpath,paste0(pop,1,'_scan_perms_bin_em.rsave')))
perms <- get(paste0('bin.em.perms.2.',1))

for (i in 2:100){
arraynum <- i
load(file.path(mpath,paste0(pop,arraynum,'_scan_perms_bin_em.rsave')))
perms <- c(perms,get(paste0('bin.em.perms.2.',i)))
}

################################################################################

cross <- fill.geno(cross, method="argmax", error.prob=0.002, map.function="kosambi", min.prob=0.95)

cross <- argmax.geno(cross, step=1, off.end=1, error.prob=0.002, map.function="kosambi", stepwidth="fixed")

cross <- calc.genoprob(cross, step=1, off.end=1, error.prob=0.002, map.function="kosambi", stepwidth="fixed")

markers <- colnames(pull.argmaxgeno(cross))

bin.em <- scanone(cross, pheno.col=4, model='binary', method = "em")

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

sens <- as.character(cross$pheno$ID[which(cross$pheno$ID %in% pheno_ind(cross,1))])
tol <- as.character(cross$pheno$ID[which(cross$pheno$ID %in% pheno_ind(cross,0))])

gts <- pull.argmaxgeno(cross)[,markers]
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
##pop_chr <- as.list(paste('CHR',gsub(":.*","",names(gt_sens)),'QTL genotype'))
##names(pop_chr) <- markers

#loc_name <- as.list(c('AIP (CHR 2)','ARNT (CHR 8)','ARNT (CHR 13)','AHRb (CHR 18)','QTLa (CHR 24)','QTLb (CHR 24)'))
#names(loc_name) <- markers

#chr <- names(gt_sens)
#chr <- ifelse(duplicated(chr),paste0(chr,'.',sum(duplicated(chr))),chr)
#names(chr) <- markers
################################################################################

single <- function(mark, pop = 'NBH'){

 ydir <- prop_sens[[mark]]
 ytot <- prop_total[[mark]]
 chromo <- pop_chr[[mark]]
 mainL <- loc_name[[mark]]
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


sapply(markers,single,pop = 'NBH')

################################################################################
## single plot AIP

cex_single <- c(.25,.5,.25) * nind(cross)
xdir <- c(1,2,3)
ydir <- tol_18

pdf(paste0("/home/jmiller1/public_html/",pdfL,".pdf"), width=3.5)
 plot(c(0.65,3.35), c(-0.1,1.1),
  xaxs="i", xaxt="n", xlab="",
  yaxs="i", yaxt="n", ylab="",
  type="n", main='AHR',
  cex.lab=1.5, cex.main=2)

  lines(xdir, ydir[as.character(c(1,2,3))],col='black',lwd=5)

  points(xdir, ydir[as.character(c(1,2,3))],col=c('black','darkblue','cornflowerblue'), pch=21,bg=c('black','darkblue','cornflowerblue'),
   cex= (12 * (ahr[as.character(c(1,2,3))] / cex_single)))

  text(xdir, ydir[as.character(c(1,2,3))], labels=ahr[as.character(c(1,2,3))],col='white',font=2, cex=2)

  mtext(side=1, line=3, 'AHR', col="black", font=2,cex=1.5)
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
cexs_hom <- cexs_hom *  nind(cross)
cexs_het <- c(0.25*0.5,0.5^2,0.25*0.5)
cexs_het <- cexs_het *  nind(cross)

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

summary(bin.em.2,perms=perms, alphas=0.15, pvalues=F, what='int')

summary(bin.em.2, thresholds=c(0, Inf,8, Inf, Inf), what="int")

        pos1  pos2 lod.full lod.fv1 lod.int lod.add lod.av1
c2 :cX  72.12  4.61     9.17    7.42    8.52   0.648  -1.098
c2 :c15 61.12  9.94    10.61    8.87    9.95   0.659  -1.087
c8 :c15 13.43 26.94     9.59    7.15    8.40   1.196  -1.246

c4 :c12 63.24 36.00    11.71    8.84    8.97   2.742  -0.133
cX :c16  7.61  3.18    10.52    7.83    8.25   2.268  -0.421
c6 :c10 34.00 74.64    10.72    7.36    8.21   2.516  -0.850
c6 :c17 48.00  7.59    11.38    8.15    9.45   1.931  -1.296




c12:c14 35.00 14.60    11.45    8.52    8.20   3.247   0.313
c13:c14 77.93 19.60    10.63    4.38    8.53   2.106  -4.152

c13:c19 62.93 35.38    10.49    4.23    9.03   1.463  -4.795

c11:c22 16.00 29.33     9.65    7.09    8.29   1.361  -1.196
c17:c22 32.59 15.33    11.08    7.85    8.35   2.730  -0.498

c3 :c23  4.26 70.00     9.51    6.64    8.30   1.209  -1.661








loc_a <- find.pseudomarker(cross,13,62.93, where = 'prob')
loc_b <- find.pseudomarker(cross,19,35.38, where = 'prob')
intxs.bin(loc_a,loc_b, popchr = 'popchr',locbN = 'locbN', main = 'main')


loc_a <- find.marker(cross,13,62.93)
loc_b <- find.marker(cross,19,35.38)



loc_a <- find.pseudomarker(cross,2,110, where = 'prob')
loc_b <- find.pseudomarker(cross,21,80, where = 'prob')
intxs.bin(loc_a,loc_b, popchr = 'popchr',locbN = 'locbN', main = 'main')

loc_a <- find.pseudomarker(cross,2,110, where = 'prob')
loc_b <- find.pseudomarker(cross,22,56.4, where = 'prob')
intxs.bin(loc_a,loc_b, popchr = 'popchr',locbN = 'locbN', main = 'main')

loc_a <- find.pseudomarker(cross,8, 60.7, where = 'prob')
loc_b <- find.pseudomarker(cross,19, 22.3, where = 'prob')
intxs.bin(loc_a,loc_b, popchr = 'popchr',locbN = 'locbN', main = 'main')

loc_a <- find.pseudomarker(cross,8,32.7, where = 'prob')
loc_b <- find.pseudomarker(cross,21,45, where = 'prob')
intxs.bin(loc_a,loc_b, popchr = 'popchr',locbN = 'locbN', main = 'main')

loc_a <- find.pseudomarker(cross,17,51.5, where = 'prob')
loc_b <- find.pseudomarker(cross,24,72.1, where = 'prob')
intxs.bin(loc_a,loc_b, popchr = 'popchr',locbN = 'locbN', main = 'main')

         pos1 pos2 lod.full lod.fv1 lod.int lod.add lod.av1
c2 :c21 110.0 80.0     17.6 -0.0194   17.57    0.00 -17.586
c2 :c22 110.0 56.4     16.8 -0.7959    8.53    8.26  -9.325
c8 :c19  60.7 22.3     14.5  8.3141    9.83    4.69  -1.521
c8 :c21  32.7 45.0     14.9  8.7378    8.46    6.48   0.274
c17:c24  51.5 72.1     15.1  7.5631    8.41    6.70  -0.846



loc_a <- find.pseudomarker(cross,2,110, where = 'prob')
loc_b <- find.pseudomarker(cross,21,79.97, where = 'prob')
intxs.bin(loc_a,loc_b, popchr = 'popchr',locbN = 'locbN', main = 'main')

loc_a <- find.pseudomarker(cross,2,110, where = 'prob')
loc_b <- find.pseudomarker(cross,22,56.44, where = 'prob')
intxs.bin(loc_a,loc_b, popchr = 'popchr',locbN = 'locbN', main = 'main')


loc_a <- find.pseudomarker(cross,2,108, where = 'prob')
loc_b <- find.pseudomarker(cross,18,38.06, where = 'prob')
intxs.bin(loc_a,loc_b, popchr = 'popchr',locbN = 'locbN', main = 'main')

###########
