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

png("/home/jmiller1/public_html/nbh_2.png", width=250)
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

png("/home/jmiller1/public_html/nbh_18.png", width=250)
 plot(c(0.65,3.35), c(-0.1,1.1),
  xaxs="i", xaxt="n", xlab="",
  yaxs="i", yaxt="n", ylab="",
  type="n", main='AHR',
  cex.lab=1.5, cex.main=2)

  lines(xdir, ydir[as.character(c(1,2,3))],col='black',lwd=5)

  points(xdir, ydir[as.character(c(1,2,3))],col=c('black','darkblue','cornflowerblue'), pch=21,bg=c('black','darkblue','cornflowerblue'),
   cex= (12 * (nbh_aip[as.character(c(1,2,3))] / cex_single)))

  text(xdir, ydir[as.character(c(1,2,3))], labels=nbh_aip[as.character(c(1,2,3))],col='white',font=2, cex=2)

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

png("/home/jmiller1/public_html/nbh_2_18.png", width=750)
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

### ELR ####













nbh_gts[which(nbh_gts[nbh_sens ,loc_a] == 1),'18:17874376']



nbh_sens_2 <- factor(,levels=c(NA,1,2,3),exclude = NULL)
nbh_sens_18 <- factor(pull.geno(nbh_cross)[nbh_sens,c()],levels=c(NA,1,2,3),exclude = NULL)

nbh_tol_2 <- factor(pull.geno(nbh_cross)[nbh_tol,c('2:35718468')],levels=c(NA,1,2,3),exclude = NULL)
nbh_tol_18 <- factor(pull.geno(nbh_cross)[nbh_tol,c('18:17874376')],levels=c(NA,1,2,3),exclude = NULL)


nbh_total <- table(nbh_sens_2) + table(nbh_tol_2)
nbh_tol_2 <- table(nbh_tol_2) / nbh_total

nbh_total <- table(nbh_sens_18) + table(nbh_tol_18)
nbh_tol_18 <- table(nbh_tol_18) / nbh_total








nbh_tol <- pull.geno(nbh_cross)[nbh_tol,c('2:35718468','18:17874376')]

nbh_sens <- apply(nbh_sens,1,factor,levels=c(NA,1,2,3),exclude = NULL)

nbh_sens <- factor(nbh_sens,levels=c(NA,1,2,3),exclude = NULL)
nbh_tol <- factor(,levels=c(NA,1,2,3),exclude = NULL)



pull.geno(nbh_cross)[,'2:35718468','18:17874376']

nbh_cross$


nbh_sens <- geno.crosstab(get_phenos(nbh_cross,1), '2:35718468', '18:17874376', eliminate.zeros=TRUE)
nbh_tol <- geno.crosstab(get_phenos(nbh_cross,0), '2:35718468', '18:17874376', eliminate.zeros=TRUE)


nbh_tol['AA',gts]
nbh_sens['AA',gts]
