#!/bin/R

pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR')]

pop <- 'ELR'

library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

# single scans
load(file.path(mpath,'single_scans.elr.rsave'))
sw_elr <- full.norm.add_only
pni_elr <- perms.norm.imp
pbe_elr <- perms.bin.em
sbe_elr <- scan.bin.em
sbm_elr <- scan.bin.mr
elr_cross <- cross

load(file.path(mpath,'single_scans.nbh.rsave'))
sw_nbh <- full.norm.add_only
pni_nbh <- perms.norm.imp
pbe_nbh <- perms.bin.em
sbe_nbh <- scan.bin.em
sbm_nbh <- scan.bin.mr
nbh_cross <- cross

load(file.path(mpath,'single_scans.new.rsave'))
sw_new <- full.norm.add_only
pni_new <- perms.norm.imp
pbe_new <- perms.bin.em
sbe_new <- scan.bin.em
sbm_new <- scan.bin.mr
new_cross <- cross


load(file.path(mpath,'single_scans.brp.rsave'))
sw_brp <- full.norm.add_only
pni_brp <- perms.norm.imp
pbe_brp <- perms.bin.em
sbe_brp <- scan.bin.em
sbm_brp <- scan.bin.mr
brp_cross <- cross

get_phenos <- function(crs,pheno){
 index <- as.character(crs$pheno$ID[which(crs$pheno$bin == pheno)])
 subset(crs,ind=index)
}

pheno_ind <- function(crs,pheno){
 as.character(crs$pheno$ID[which(crs$pheno$bin == pheno)])
}

################################################################################
### NBH ####
nbh_sens <- as.character(nbh_cross$pheno$ID[which(nbh_cross$pheno$ID %in% pheno_ind(nbh_cross,1))])
nbh_tol <- as.character(nbh_cross$pheno$ID[which(nbh_cross$pheno$ID %in% pheno_ind(nbh_cross,0))])

## From MR
loc_a <- '2:35718468'
loc_b <- '18:17874376'

## EM
loc_a <- '2:35740383'
loc_b <- '18:17665127'

nbh_gts <- pull.geno(nbh_cross)[,c(loc_a,loc_b)]
rownames(nbh_gts) <- as.character(nbh_cross$pheno$ID)
################################################################################
## single locus
nbh_sens_2 <- factor(nbh_gts[nbh_sens,c(loc_a)],levels=c(NA,1,2,3),exclude = NULL)
nbh_sens_18 <- factor(nbh_gts[nbh_sens,c(loc_b)],levels=c(NA,1,2,3),exclude = NULL)

nbh_tol_2 <- factor(nbh_gts[nbh_tol,c(loc_a)],levels=c(NA,1,2,3),exclude = NULL)
nbh_tol_18 <- factor(nbh_gts[nbh_tol,c(loc_b)],levels=c(NA,1,2,3),exclude = NULL)

nbh_aip <- table(nbh_sens_2) + table(nbh_tol_2)
nbh_ahr <- table(nbh_sens_18) + table(nbh_tol_18)

nbh_2_total <- table(nbh_sens_2) + table(nbh_tol_2)
nbh_18_total <- table(nbh_sens_18) + table(nbh_tol_18)

nbh_tol_2 <- table(nbh_sens_2) / nbh_2_total
nbh_tol_18 <- table(nbh_sens_18) / nbh_18_total
################################################################################
## single plot AIP

cex_single <- c(.25,.5,.25) * 94
xdir <- c(1,2,3)
ydir <- nbh_tol_2

pdf("/home/jmiller1/public_html/nbh_2.pdf", width=3.5)
 plot(c(0.65,3.35), c(-0.1,1.1),
  xaxs="i", xaxt="n", xlab="",
  yaxs="i", yaxt="n", ylab="",
  type="n", main='AIP',
  cex.lab=1.5, cex.main=2)

  lines(xdir, ydir[as.character(c(1,2,3))],col='black',lwd=5)

  points(xdir, ydir[as.character(c(1,2,3))],col=c('black','darkblue','cornflowerblue'), pch=21,bg=c('black','darkblue','cornflowerblue'),
   cex= (12 * (nbh_aip[as.character(c(1,2,3))] / cex_single)))

  text(xdir, ydir[as.character(c(1,2,3))], labels=nbh_aip[as.character(c(1,2,3))],col='white',font=2, cex=2)

  mtext(side=1, line=3, 'AIP', col="black", font=2,cex=1.5)
  mtext(side=2, line=3, "Proportion Deformed", col="black", font=2, cex=1.5)

  axis(side=1, at=c(1,2,3), labels=c('AA','AB','BB'),font=2, cex.axis=2)
  axis(side=2, at=c(0,0.5,1.0), labels=c('0','0.5','1.0'),font=2, cex.lab=2, cex.axis=2)
dev.off()
################################################################################
################################################################################
## single plot AIP

cex_single <- c(.25,.5,.25) * 94
xdir <- c(1,2,3)
ydir <- nbh_tol_18

pdf("/home/jmiller1/public_html/nbh_18.pdf", width=3.5)
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
################################################################################

## interactions
nbh_2_AA <- table(factor(nbh_gts[names(which(nbh_gts[nbh_sens ,loc_a] == 1)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
nbh_2_AB <- table(factor(nbh_gts[names(which(nbh_gts[nbh_sens ,loc_a] == 2)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
nbh_2_BB <- table(factor(nbh_gts[names(which(nbh_gts[nbh_sens ,loc_a] == 3)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
nbh_sensit <- rbind(nbh_2_AA,nbh_2_AB,nbh_2_BB)
colnames(nbh_sensit) <- c('NA','AA','AB','BB')

nbh_2_AA <- table(factor(nbh_gts[names(which(nbh_gts[nbh_tol  ,loc_a] == 1)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
nbh_2_AB <- table(factor(nbh_gts[names(which(nbh_gts[nbh_tol  ,loc_a] == 2)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
nbh_2_BB <- table(factor(nbh_gts[names(which(nbh_gts[nbh_tol  ,loc_a] == 3)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
nbh_resist <- rbind(nbh_2_AA,nbh_2_AB,nbh_2_BB)
colnames(nbh_resist) <- c('NA','AA','AB','BB')

nbh_2_AA <- table(factor(nbh_gts[names(which(nbh_gts[c(nbh_sens,nbh_tol) ,loc_a] == 1)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
nbh_2_AB <- table(factor(nbh_gts[names(which(nbh_gts[c(nbh_sens,nbh_tol) ,loc_a] == 2)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
nbh_2_BB <- table(factor(nbh_gts[names(which(nbh_gts[c(nbh_sens,nbh_tol) ,loc_a] == 3)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
nbh_total <- rbind(nbh_2_AA, nbh_2_AB, nbh_2_BB)
colnames(nbh_total) <- c('NA','AA','AB','BB')
rownames(nbh_total) <- rownames(nbh_resist)

## Plot
xdir <- c(1,2,3)
ydir <- nbh_sensit/nbh_total

#colnames(ydir) <- c('NA','AA','AB','BB')
#rownames(ydir) <- c('AA','AB','BB')
#rownames(nbh_total) <- c('AA','AB','BB')

cexs_hom <- c(0.25^2,0.5^2,0.25^2)
cexs_hom <- cexs_hom * 94
cexs_het <- c(0.25*0.5,0.5^2,0.25*0.5)
cexs_het <- cexs_het * 94

pdf("/home/jmiller1/public_html/nbh_2_18.pdf", width=10)
 plot(c(0.65,3.35), c(-0.1,1.1),
  xaxs="i", xaxt="n", xlab="",
  yaxs="i", yaxt="n", ylab="",
  type="n", main='AIP x AHR (NBH)',
  cex.lab=1.5, cex.main=2)

  rect(1.5, -0.1, 2.5, 1.1,col='lightgrey',border = 'transparent')

  lines(xdir-0.28, ydir[rownames(nbh_total),'BB'],col='cornflowerblue',lwd=5)
  lines(xdir, ydir[rownames(nbh_total),'AB'],col='darkblue',lwd=5)
  lines(xdir+0.28, ydir[rownames(nbh_total),'AA'],col='black',lwd=5)

  points(xdir, ydir[rownames(nbh_total),'AB'],col='darkblue', pch=21, bg='darkblue',
   cex= (12 * (nbh_total[rownames(nbh_total),'AB'] / cexs_het))+2)
  points(xdir-0.28, ydir[rownames(nbh_total),'BB'], col='cornflowerblue', pch=21, bg='cornflowerblue',
   cex= (12 * (nbh_total[rownames(nbh_total),'BB'] / cexs_hom))+2)
  points(xdir+0.28, ydir[rownames(nbh_total),'AA'],col='black', pch=21,bg='black',
   cex= (12 * (nbh_total[rownames(nbh_total),'AA'] / cexs_hom))+2)


  text(xdir-0.28, ydir[rownames(nbh_total),'BB'], labels=nbh_total[rownames(nbh_total),'BB'],col='white',font=2, cex=2)
  text(xdir, ydir[rownames(nbh_total),'AB'], labels=nbh_total[rownames(nbh_total),'AB'],col='white',font=2, cex=2)
  text(xdir+0.28, ydir[rownames(nbh_total),'AA'], labels=nbh_total[rownames(nbh_total),'AA'],col='white',font=2, cex=2)

  mtext(side=1, line=3, 'AHR Genotype', col="black", font=2,cex=1.5)
  mtext(side=2, line=3, "Proportion Deformed", col="black", font=2, cex=1.5)

  axis(side=1, at=c(1,2,3), labels=c('AA','AB','BB'),font=2, cex.axis=2)
  axis(side=2, at=c(0,0.5,1.0), labels=c('0','0.5','1.0'),font=2, cex.lab=2, cex.axis=2)
dev.off()
################################################################################





################################################################################
### ELR ####
################################################################################
################################################################################
################################################################################
elr_sens <- as.character(elr_cross$pheno$ID[which(elr_cross$pheno$ID %in% pheno_ind(elr_cross,1))])
elr_tol <- as.character(elr_cross$pheno$ID[which(elr_cross$pheno$ID %in% pheno_ind(elr_cross,0))])

em <- rownames(summary(sbe_elr)[c(13,18),])
mr <- rownames(summary(sbm_elr)[c(13,18),])

####
## what is the parent genotypes at 13?
## "13:1009168" to check nearby locus

prefilt_par_gts <- pull.markers(cross,unique(c(em,mr,"13:1009168")))
postfilt_off_gt <- pull.markers(elr_cross,unique(c(em,mr,"13:1009168")))

prefilt_par_gts <- cbind(as.character(prefilt_par_gts$pheno$ID), prefilt_par_gts$geno$'13'$data)
rownames(prefilt_par_gts) <- prefilt_par_gts[,1]
prefilt_par_gts[c('ELR_ER1124F','BLI_BI1124M'),]

postfilt_off_gt <- cbind(as.character(postfilt_off_gt$pheno$ID), postfilt_off_gt$geno$'13'$data)
rownames(postfilt_off_gt) <- postfilt_off_gt[,1]

is_phase_same <- prefilt_par_gts[rownames(postfilt_off_gt),2] == postfilt_off_gt[rownames(postfilt_off_gt),2]

### Sensitive is BB, indicating that a fixed resistant allele are all deformed


## From MR
loc_a <- mr[1]
loc_b <- mr[2]

## EM
loc_a <- em[1]
loc_b <- em[2]

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

cex_single <- c(.25,.5,.25) * 94
xdir <- c(1,2,3)
ydir <- elr_tol_13

png("/home/jmiller1/public_html/elr_13.png", width=250)
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

cex_single <- c(.25,.5,.25) * 94
xdir <- c(1,2,3)
ydir <- elr_tol_18

png("/home/jmiller1/public_html/elr_18.png", width=250)
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


## interactions
elr_2_AA <- table(factor(elr_gts[names(which(elr_gts[elr_sens ,loc_a] == 1)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
elr_2_AB <- table(factor(elr_gts[names(which(elr_gts[elr_sens ,loc_a] == 2)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
elr_2_BB <- table(factor(elr_gts[names(which(elr_gts[elr_sens ,loc_a] == 3)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
elr_sensit <- rbind(elr_2_AA,elr_2_AB,elr_2_BB)
colnames(elr_sensit) <- c('NA','AA','AB','BB')

elr_2_AA <- table(factor(elr_gts[names(which(elr_gts[elr_tol  ,loc_a] == 1)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
elr_2_AB <- table(factor(elr_gts[names(which(elr_gts[elr_tol  ,loc_a] == 2)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
elr_2_BB <- table(factor(elr_gts[names(which(elr_gts[elr_tol  ,loc_a] == 3)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
elr_resist <- rbind(elr_2_AA,elr_2_AB,elr_2_BB)
colnames(elr_resist) <- c('NA','AA','AB','BB')

elr_2_AA <- table(factor(elr_gts[names(which(elr_gts[c(elr_sens,elr_tol) ,loc_a] == 1)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
elr_2_AB <- table(factor(elr_gts[names(which(elr_gts[c(elr_sens,elr_tol) ,loc_a] == 2)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
elr_2_BB <- table(factor(elr_gts[names(which(elr_gts[c(elr_sens,elr_tol) ,loc_a] == 3)),loc_b],levels=c(NA,1,2,3),exclude = NULL))
elr_total <- rbind(elr_2_AA, elr_2_AB, elr_2_BB)
colnames(elr_total) <- c('NA','AA','AB','BB')
rownames(elr_total) <- rownames(elr_resist)

## Plot
xdir <- c(1,2,3)
ydir <- elr_sensit/elr_total

#colnames(ydir) <- c('NA','AA','AB','BB')
#rownames(ydir) <- c('AA','AB','BB')
#rownames(elr_total) <- c('AA','AB','BB')

cexs_hom <- c(0.25^2,0.5^2,0.25^2)
cexs_hom <- cexs_hom * sum(elr_total)
cexs_het <- c(0.25*0.5,0.5^2,0.25*0.5)
cexs_het <- cexs_het * sum(elr_total)

png("/home/jmiller1/public_html/elr_13_18.png", width=750)
 plot(c(0.65,3.35), c(-0.1,1.1),
  xaxs="i", xaxt="n", xlab="",
  yaxs="i", yaxt="n", ylab="",
  type="n", main='QTL13 x AHR',
  cex.lab=1.5, cex.main=2)

  rect(1.5, -0.1, 2.5, 1.1,col='lightgrey',border = 'transparent')

  lines(xdir-0.28, ydir[rownames(elr_total),'BB'],col='cornflowerblue',lwd=5)
  lines(xdir, ydir[rownames(elr_total),'AB'],col='darkblue',lwd=5)
  lines(xdir+0.28, ydir[rownames(elr_total),'AA'],col='black',lwd=5)

  points(xdir, ydir[rownames(elr_total),'AB'],col='darkblue', pch=21, bg='darkblue',
   cex= (12 * (elr_total[rownames(elr_total),'AB'] / cexs_het))+1)

  points(xdir-0.28, ydir[rownames(elr_total),'BB'], col='cornflowerblue', pch=21, bg='cornflowerblue',
   cex= (12 * (elr_total[rownames(elr_total),'BB'] / cexs_hom))+1)

  points(xdir+0.28, ydir[rownames(elr_total),'AA'],col='black', pch=21,bg='black',
   cex= (12 * (elr_total[rownames(elr_total),'AA'] / cexs_hom))+1)


  text(xdir-0.28, ydir[rownames(elr_total),'BB'], labels=elr_total[rownames(elr_total),'BB'],col='white',font=2, cex=2)
  text(xdir, ydir[rownames(elr_total),'AB'], labels=elr_total[rownames(elr_total),'AB'],col='white',font=2, cex=2)
  text(xdir+0.28, ydir[rownames(elr_total),'AA'], labels=elr_total[rownames(elr_total),'AA'],col='white',font=2, cex=2)

  mtext(side=1, line=3, 'QTL13 Genotype', col="black", font=2,cex=1.5)
  mtext(side=2, line=3, "Proportion Deformed", col="black", font=2, cex=1.5)

  axis(side=1, at=c(1,2,3), labels=c('AA','AB','BB'),font=2, cex.axis=2)
  axis(side=2, at=c(0,0.5,1.0), labels=c('0','0.5','1.0'),font=2, cex.lab=2, cex.axis=2)
dev.off()
