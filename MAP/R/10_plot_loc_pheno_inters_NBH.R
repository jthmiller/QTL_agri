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

cross <- sim.geno(cross,step=0,off.end=5, error.prob=0.025,map.function="kosambi")
cross <- calc.genoprob(cross,step=0,error.prob=0.025,off.end=5)

bin.add.imp <- stepwiseqtl(cross, incl.markers=T, additive.only = T, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=3)
bin.add.imp.qtls <- summary(bin.add.imp)
bin.add.imp.qtls <- makeqtl(gg_step2, chr=as.character(bin.add.imp.qtls$chr), pos=as.numeric(bin.add.imp.qtls$pos), what="draws")


norm.add <- stepwiseqtl(cross, incl.markers=T, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=5)
norm.add.qtls <- summary(norm.add)
norm.add.qtls <- makeqtl(cross, chr=as.character(norm.add.qtls$chr), pos=as.numeric(norm.add.qtls$pos), what="draws")

loc_a <- find.marker(cross, norm.add.qtls$chr[1],norm.add.qtls$pos[1])
loc_b <- find.marker(cross, norm.add.qtls$chr[2],norm.add.qtls$pos[2])


get_phenos <- function(crs,pheno){
 index <- as.character(crs$pheno$ID[which(crs$pheno$bin == pheno)])
 subset(crs,ind=index)
}

pheno_ind <- function(crs,pheno){
 as.character(crs$pheno$ID[which(crs$pheno$bin == pheno)])
}

################################################################################
### NBH ####
nbh_cross <- cross
nbh_sens <- as.character(nbh_cross$pheno$ID[which(nbh_cross$pheno$ID %in% pheno_ind(nbh_cross,1))])
nbh_tol <- as.character(nbh_cross$pheno$ID[which(nbh_cross$pheno$ID %in% pheno_ind(nbh_cross,0))])


#### From MR
##loc_a <- '2:35718468'
##loc_b <- '18:17874376'
##
#### EM
##loc_a <- '2:35740383'
##loc_b <- '18:17665127'

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

  mtext(side=1, line=3, 'AHR Genotype', col="black", font=2,cex=1.5)
  mtext(side=2, line=3, "Proportion Deformed", col="black", font=2, cex=1.5)

  axis(side=1, at=c(1,2,3), labels=c('AA','AB','BB'),font=2, cex.axis=2)
  axis(side=2, at=c(0,0.5,1.0), labels=c('0','0.5','1.0'),font=2, cex.lab=2, cex.axis=2)
dev.off()
################################################################################
