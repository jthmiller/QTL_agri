#!/bin/R

## pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]

pop <- 'ELR'
library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

## Read cross
cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

cross <- sim.geno(cross,step=0,off.end=5, error.prob=0.025,map.function="kosambi")
cross <- calc.genoprob(cross,step=0,error.prob=0.025,off.end=5)

norm.add <- stepwiseqtl(cross, incl.markers=T, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=5)
norm.add.qtls <- summary(norm.add)
norm.add.qtls <- makeqtl(cross, chr=as.character(norm.add.qtls$chr), pos=as.numeric(norm.add.qtls$pos), what="draws")

loc_a <- find.marker(cross, norm.add.qtls$chr[2],norm.add.qtls$pos[2])
loc_b <- find.marker(cross, norm.add.qtls$chr[3],norm.add.qtls$pos[3])


scan.norm.imp <- scanone(cross, method = "imp", model = "normal", pheno.col = 5)
plotpub('scan.norm.imp.elr')
plot(scan.norm.imp)
dev.off()

################################################################################
### ELR ####
################################################################################
################################################################################
################################################################################
elr_cross <- cross

elr_sens <- as.character(elr_cross$pheno$ID[which(elr_cross$pheno$ID %in% pheno_ind(elr_cross,1))])
elr_tol <- as.character(elr_cross$pheno$ID[which(elr_cross$pheno$ID %in% pheno_ind(elr_cross,0))])

#em <- rownames(summary(sbe_elr)[c(13,18),])
#mr <- rownames(summary(sbm_elr)[c(13,18),])

####
## what is the parent genotypes at 13?
## "13:1009168" to check nearby locus

##prefilt_par_gts <- pull.markers(cross,unique(c(em,mr,"13:1009168")))
##postfilt_off_gt <- pull.markers(elr_cross,unique(c(em,mr,"13:1009168")))
##
##prefilt_par_gts <- cbind(as.character(prefilt_par_gts$pheno$ID), prefilt_par_gts$geno$'13'$data)
##rownames(prefilt_par_gts) <- prefilt_par_gts[,1]
##prefilt_par_gts[c('ELR_ER1124F','BLI_BI1124M'),]
##
##postfilt_off_gt <- cbind(as.character(postfilt_off_gt$pheno$ID), postfilt_off_gt$geno$'13'$data)
##rownames(postfilt_off_gt) <- postfilt_off_gt[,1]
##
##is_phase_same <- prefilt_par_gts[rownames(postfilt_off_gt),2] == postfilt_off_gt[rownames(postfilt_off_gt),2]
##
##### Sensitive is BB, indicating that a fixed resistant allele are all deformed


elr_gts <- pull.geno(elr_cross)[,c(loc_a,loc_b)]
rownames(elr_gts) <- as.character(elr_cross$pheno$ID)
################################################################################
## single locus
elr_sens_13 <- factor(elr_gts[elr_sens,c(loc_a)],levels=c(NA,1,2,3),exclude = NULL)
elr_sens_18 <- factor(elr_gts[elr_sens,c(loc_b)],levels=c(NA,1,2,3),exclude = NULL)

elr_tol_13 <- factor(elr_gts[elr_tol,c(loc_a)],levels=c(NA,1,2,3),exclude = NULL)
elr_tol_18 <- factor(elr_gts[elr_tol,c(loc_b)],levels=c(NA,1,2,3),exclude = NULL)

elr_aip <- table(elr_sens_13) + table(elr_tol_13)
elr_ahr <- table(elr_sens_18) + table(elr_tol_18)

elr_13_total <- table(elr_sens_13) + table(elr_tol_13)
elr_18_total <- table(elr_sens_18) + table(elr_tol_18)

elr_tol_13 <- table(elr_sens_13) / elr_13_total
elr_tol_18 <- table(elr_sens_18) / elr_18_total
################################################################################
## single plot 13

cex_single <- c(.25,.5,.25) * 79
xdir <- c(1,2,3)
ydir <- elr_tol_13

pdf("/home/jmiller1/public_html/elr_13.pdf", width=3.5)
 plot(c(0.65,3.35), c(-0.1,1.1),
  xaxs="i", xaxt="n", xlab="",
  yaxs="i", yaxt="n", ylab="",
  type="n", main='QTL13',
  cex.lab=1.5, cex.main=2)

  lines(xdir, ydir[as.character(c(1,2,3))],col='black',lwd=5)

  points(xdir, ydir[as.character(c(1,2,3))],col=c('black','darkblue','cornflowerblue'), pch=21,bg=c('black','darkblue','cornflowerblue'),
   cex= (12 * (elr_aip[as.character(c(1,2,3))] / cex_single)))

  text(xdir, ydir[as.character(c(1,2,3))], labels=elr_aip[as.character(c(1,2,3))],col='white',font=2, cex=2)

  ##mtext(side=1, line=3, 'QTL13', col="black", font=2,cex=1.5)
  mtext(side=2, line=3, "Proportion Deformed", col="black", font=2, cex=1.5)

  axis(side=1, at=c(1,2,3), labels=c('AA','AB','BB'),font=2, cex.axis=2)
  axis(side=2, at=c(0,0.5,1.0), labels=c('0','0.5','1.0'),font=2, cex.lab=2, cex.axis=2)
dev.off()
################################################################################
################################################################################
## single plot AHR

cex_single <- c(.25,.5,.25) * 79
xdir <- c(1,2,3)
ydir <- elr_tol_18

pdf("/home/jmiller1/public_html/elr_18.pdf", width=3.5)
 plot(c(0.65,3.35), c(-0.1,1.1),
  xaxs="i", xaxt="n", xlab="",
  yaxs="i", yaxt="n", ylab="",
  type="n", main='AHR',
  cex.lab=1.5, cex.main=2)

  lines(xdir, ydir[as.character(c(1,2,3))],col='black',lwd=5)

  points(xdir, ydir[as.character(c(1,2,3))],col=c('black','darkblue','cornflowerblue'), pch=21,bg=c('black','darkblue','cornflowerblue'),
   cex= (12 * (elr_aip[as.character(c(1,2,3))] / cex_single)))

  text(xdir, ydir[as.character(c(1,2,3))], labels=elr_aip[as.character(c(1,2,3))],col='white',font=2, cex=2)

  ##mtext(side=1, line=3, 'AHR', col="black", font=2,cex=1.5)
  mtext(side=2, line=3, "Proportion Deformed", col="black", font=2, cex=1.5)

  axis(side=1, at=c(1,2,3), labels=c('AA','AB','BB'),font=2, cex.axis=2)
  axis(side=2, at=c(0,0.5,1.0), labels=c('0','0.5','1.0'),font=2, cex.lab=2, cex.axis=2)
dev.off()
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
