#!/bin/R
##pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('new','BRP','NEW','ELR')]

pop <- 'NEW'

library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)


load(file.path(mpath,'single_scans.new.rsave'))

new_cross <- cross
new_cross <- switchAlleles(new_cross, markers = markernames(new_cross,chr=2))
new_cross <- switchAlleles(new_cross, markers = markernames(new_cross,chr=18))


new_cross  <- sim.geno(new_cross)
new_cross  <- calc.genoprob(new_cross ,step=1,error.prob=0.01,off.end=5)

scan.bin.em <- scanone(new_cross, method = "em", model = "binary", pheno.col = 4)
scan.bin.imp <- scanone(new_cross, method = "imp", model = "binary", pheno.col = 4)
scan.bin.mr <- scanone(new_cross, method = "mr", model = "binary", pheno.col = 4)
full.norm.add_only <- stepwiseqtl(new_cross, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=8)

sw_new <- full.norm.add_only
pni_new <- perms.norm.imp
pbe_new <- perms.bin.em
sbe_new <- scan.bin.em
sbm_new <- scan.bin.mr

get_phenos <- function(crs,pheno){
 index <- as.character(crs$pheno$ID[which(crs$pheno$bin == pheno)])
 subset(crs,ind=index)
}

pheno_ind <- function(crs,pheno){
 as.character(crs$pheno$ID[which(crs$pheno$bin == pheno)])
}

################################################################################
### new ####
new_sens <- as.character(new_cross$pheno$ID[which(new_cross$pheno$ID %in% pheno_ind(new_cross,1))])
new_tol <- as.character(new_cross$pheno$ID[which(new_cross$pheno$ID %in% pheno_ind(new_cross,0))])

em <- rownames(summary(sbe_new)[c(2,18),])
mr <- rownames(summary(sbm_new)[c(2,18),])

## From MR
loc_a <- mr[1]
loc_b <- mr[2]

## EM
loc_a <- em[1]
loc_b <- em[2]


################################################################################
new_gts <- pull.geno(new_cross)[,c(loc_a,loc_b)]
rownames(new_gts) <- as.character(new_cross$pheno$ID)
################################################################################
## single locus
new_sens_2 <- factor(new_gts[new_sens,c(loc_a)],levels=c(NA,1,2,3),exclude = NULL)
new_sens_18 <- factor(new_gts[new_sens,c(loc_b)],levels=c(NA,1,2,3),exclude = NULL)

new_tol_2 <- factor(new_gts[new_tol,c(loc_a)],levels=c(NA,1,2,3),exclude = NULL)
new_tol_18 <- factor(new_gts[new_tol,c(loc_b)],levels=c(NA,1,2,3),exclude = NULL)

new_aip <- table(new_sens_2) + table(new_tol_2)
new_ahr <- table(new_sens_18) + table(new_tol_18)

new_2_total <- table(new_sens_2) + table(new_tol_2)
new_18_total <- table(new_sens_18) + table(new_tol_18)

new_tol_2 <- table(new_sens_2) / new_2_total
new_tol_18 <- table(new_sens_18) / new_18_total
################################################################################
## single plot AIP

cex_single <- c(.25,.5,.25) * 96
xdir <- c(1,2,3)
ydir <- new_tol_2

pdf("/home/jmiller1/public_html/new_2.pdf", width=3.5)
 plot(c(0.65,3.35), c(-0.1,1.1),
  xaxs="i", xaxt="n", xlab="",
  yaxs="i", yaxt="n", ylab="",
  type="n", main='AIP',
  cex.lab=1.5, cex.main=2)

  lines(xdir, ydir[as.character(c(1,2,3))],col='black',lwd=5)

  points(xdir, ydir[as.character(c(1,2,3))],col=c('black','darkblue','cornflowerblue'), pch=21,bg=c('black','darkblue','cornflowerblue'),
   cex= (12 * (new_aip[as.character(c(1,2,3))] / cex_single)))

  text(xdir, ydir[as.character(c(1,2,3))], labels=new_aip[as.character(c(1,2,3))],col='white',font=2, cex=2)

  mtext(side=1, line=3, 'AIP', col="black", font=2,cex=1.5)
  mtext(side=2, line=3, "Proportion Deformed", col="black", font=2, cex=1.5)

  axis(side=1, at=c(1,2,3), labels=c('AA','AB','BB'),font=2, cex.axis=2)
  axis(side=2, at=c(0,0.5,1.0), labels=c('0','0.5','1.0'),font=2, cex.lab=2, cex.axis=2)
dev.off()
################################################################################
################################################################################
## single plot AIP

cex_single <- c(.25,.5,.25) * 96
xdir <- c(1,2,3)
ydir <- new_tol_18

pdf("/home/jmiller1/public_html/new_18.pdf", width=3.5)
 plot(c(0.65,3.35), c(-0.1,1.1),
  xaxs="i", xaxt="n", xlab="",
  yaxs="i", yaxt="n", ylab="",
  type="n", main='AHR',
  cex.lab=1.5, cex.main=2)

  lines(xdir, ydir[as.character(c(1,2,3))],col='black',lwd=5)

  points(xdir, ydir[as.character(c(1,2,3))],col=c('black','darkblue','cornflowerblue'), pch=21,bg=c('black','darkblue','cornflowerblue'),
   cex= (12 * (new_ahr[as.character(c(1,2,3))] / cex_single)))

  text(xdir, ydir[as.character(c(1,2,3))], labels=new_ahr[as.character(c(1,2,3))],col='white',font=2, cex=2)

  mtext(side=1, line=3, 'AHR', col="black", font=2,cex=1.5)
  mtext(side=2, line=3, "Proportion Deformed", col="black", font=2, cex=1.5)

  axis(side=1, at=c(1,2,3), labels=c('AA','AB','BB'),font=2, cex.axis=2)
  axis(side=2, at=c(0,0.5,1.0), labels=c('0','0.5','1.0'),font=2, cex.lab=2, cex.axis=2)
dev.off()
################################################################################

## interactions
new_2_AA <- table(factor(new_gts[names(which(new_gts[new_sens ,loc_a] == 1)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
new_2_AB <- table(factor(new_gts[names(which(new_gts[new_sens ,loc_a] == 2)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
new_2_BB <- table(factor(new_gts[names(which(new_gts[new_sens ,loc_a] == 3)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
new_sensit <- rbind(new_2_AA,new_2_AB,new_2_BB)
colnames(new_sensit) <- c('NA','AA','AB','BB')

new_2_AA <- table(factor(new_gts[names(which(new_gts[new_tol  ,loc_a] == 1)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
new_2_AB <- table(factor(new_gts[names(which(new_gts[new_tol  ,loc_a] == 2)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
new_2_BB <- table(factor(new_gts[names(which(new_gts[new_tol  ,loc_a] == 3)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
new_resist <- rbind(new_2_AA,new_2_AB,new_2_BB)
colnames(new_resist) <- c('NA','AA','AB','BB')

new_2_AA <- table(factor(new_gts[names(which(new_gts[c(new_sens,new_tol) ,loc_a] == 1)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
new_2_AB <- table(factor(new_gts[names(which(new_gts[c(new_sens,new_tol) ,loc_a] == 2)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
new_2_BB <- table(factor(new_gts[names(which(new_gts[c(new_sens,new_tol) ,loc_a] == 3)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
new_total <- rbind(new_2_AA, new_2_AB, new_2_BB)
colnames(new_total) <- c('NA','AA','AB','BB')
rownames(new_total) <- rownames(new_resist)

## Plot
xdir <- c(1,2,3)
ydir <- new_sensit/new_total

#colnames(ydir) <- c('NA','AA','AB','BB')
#rownames(ydir) <- c('AA','AB','BB')
#rownames(new_total) <- c('AA','AB','BB')

cexs_hom <- c(0.25^2,0.5^2,0.25^2)
cexs_hom <- cexs_hom * 96
cexs_het <- c(0.25*0.5,0.5^2,0.25*0.5)
cexs_het <- cexs_het * 96

pdf("/home/jmiller1/public_html/new_2_18.pdf", width=10)
 plot(c(0.65,3.35), c(-0.1,1.1),
  xaxs="i", xaxt="n", xlab="",
  yaxs="i", yaxt="n", ylab="",
  type="n", main='AIP x AHR',
  cex.lab=1.5, cex.main=2)

  rect(1.5, -0.1, 2.5, 1.1,col='lightgrey',border = 'transparent')

  lines(xdir-0.28, ydir[rownames(new_total),'BB'],col='cornflowerblue',lwd=5)
  lines(xdir, ydir[rownames(new_total),'AB'],col='darkblue',lwd=5)
  lines(xdir+0.28, ydir[rownames(new_total),'AA'],col='black',lwd=5)

  points(xdir, ydir[rownames(new_total),'AB'],col='darkblue', pch=21, bg='darkblue',
   cex= (12 * (new_total[rownames(new_total),'AB'] / cexs_het))+2)
  points(xdir-0.28, ydir[rownames(new_total),'BB'], col='cornflowerblue', pch=21, bg='cornflowerblue',
   cex= (12 * (new_total[rownames(new_total),'BB'] / cexs_hom))+2)
  points(xdir+0.28, ydir[rownames(new_total),'AA'],col='black', pch=21,bg='black',
   cex= (12 * (new_total[rownames(new_total),'AA'] / cexs_hom))+2)


  text(xdir-0.28, ydir[rownames(new_total),'BB'], labels=new_total[rownames(new_total),'BB'],col='white',font=2, cex=2)
  text(xdir, ydir[rownames(new_total),'AB'], labels=new_total[rownames(new_total),'AB'],col='white',font=2, cex=2)
  text(xdir+0.28, ydir[rownames(new_total),'AA'], labels=new_total[rownames(new_total),'AA'],col='white',font=2, cex=2)

  mtext(side=1, line=3, 'AHR Genotype', col="black", font=2,cex=1.5)
  mtext(side=2, line=3, "Proportion Deformed", col="black", font=2, cex=1.5)

  axis(side=1, at=c(1,2,3), labels=c('AA','AB','BB'),font=2, cex.axis=2)
  axis(side=2, at=c(0,0.5,1.0), labels=c('0','0.5','1.0'),font=2, cex.lab=2, cex.axis=2)
dev.off()
################################################################################
