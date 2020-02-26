#!/bin/R
### first run combine pops for multi-pop cross objects
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
library("ggridges")
library("plyr")
library("scales")
library("ggrepel")
library('qtl')
library('RColorBrewer')
library('parallel')
mpath <- '/home/jmiller1/QTL_agri/data'
setwd(mpath)

##cores <- detectCores() - 2
##
##fl <- file.path(mpath,paste0(pop'.mapped.tsp.csv'))
##
##cross <- read.cross(
## file = fl,
## format = "csv", genotypes=c("1","2","3"),
## estimate.map = FALSE
##)
################################################################################
################################################################################

#load(file.path(mpath,'supplemental_plot_env.rsave'))
#cross <- cross_NBH
#cross <- cross_ELR
#cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno),5)

################################################################################
################################################################################

erp <- 0.0025

cross <- jittermap(cross)

cross <- sim.geno(cross, step=0.5, error.prob=erp, off.end=5,
 map.function="kosambi", n.draws=100, stepwidth="fixed")

cross <- calc.genoprob(cross, step=0.5, error.prob=erp, off.end=5,
 map.function="kosambi", stepwidth="fixed")

norm <- scanone(cross, method = "imp", model = "normal", pheno.col = 5)
bin <- scanone(cross, method = "em", model = "binary", pheno.col = 4)
gt <- geno.table(cross)
map <- pull.map(cross)
map_sum <- summary(pull.map(cross))

################################################################################

cross_grid <- reduce2grid(cross)

bin_grid <- scanone(cross_grid, method = "em", model = "binary", pheno.col = 4)
norm_grid <- scanone(cross_grid, method = "imp", model = "normal", pheno.col = 5)
gt_grid <- geno.table(cross_grid)
map_grid <- pull.map(cross_grid)
map_grid_sum <- summary(pull.map(cross_grid))

################################################################################
################################################################################

col <- c("slateblue", "violetred", "green3")
cross2 <- convert2cross2(cross)
map2 <- insert_pseudomarkers(cross2$gmap, step=0.5)
pr2 <- calc_genoprob(cross2, map2, error_prob=0.0025, cores=cores)
bin2 <- scan1(pr2, pheno=cross2$pheno[,'bin'] , model = "binary", cores = cores)

################################################################################
if(pop == 'NBH') { rank <- nbh.rank ; ahr <- ahr_nbh ;  pbsname <- 'NBH'; pfstNSname <- 'F.NBH' ; piname <- pfstname <- 'BI.NBH'}
if(pop == 'ELR') { rank <- elr.rank ; ahr <- ahr_elr ;  pbsname <- 'ER'; pfstNSname <- 'ER.SH' ; piname <- pfstname <- 'ER.KC'}
################################################################################
#bottom, left, top and right
#location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks
stat <- 'pfst'
for (ch in 1:24){

 pdf(paste0("/home/jmiller1/public_html/",pop,"_all",ch,"segdist.pdf"), width=4.5,height=6)
 mat<-matrix(c(1:8),8,1, byrow=T)
 layout(mat, widths=1, heights= c(0.1, 0.1, 0.1, 0.025, 0.1, 0.1, 0.1, 0.05))
 par(mar=c(0.5,2.25,0,1)+0.1, oma = c(1, 1, 0, 1))

 len <- map_sum[ch,'length']
 ind <- which(gt$chr == ch)
 segX <- map[[ch]][rownames(gt)[ind]]
 segY <- -log10(gt[ind,'P.value'])
 out <- which(!gt$chr == ch)
 segB <- -log10(gt[out,'P.value'])
 segA <- rescale(unlist(map, use.names=F),to = c(0,len))[out]
 ab <- ahr[which(ahr$chr==ch),'pos1']

 ################################################################################
 plot_pgen(
  crs = cross, chrs=ch, stat = pfst, map = 'mid', mgp = c(1.25, 0.5, 0),
  ahr = ahr, ahr_clm= 'stp',  colnm = pfstNSname, popgen = rank,
  rank_clm='end', ylimo=c(-0.025,0.8), stat_name='pfst NxS', pch=16,
  cex=0.15, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75, xaxt="n"
 )

 if(stat == 'pfst'){
  plot_pgen(
   crs = cross, chrs=ch, stat = pfst, map = 'mid', mgp = c(1.25, 0.5, 0),
   ahr = ahr, ahr_clm= 'stp',  colnm = pfstname, popgen = rank,
   rank_clm='end', ylimo=c(-0.025,0.6), stat_name='pfst', pch=16,
   cex=0.15, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75, xaxt="n"
  )
 } else {
  plot_pgen(
   crs = cross, chrs=ch, stat = pbs, map = 'mid', mgp = c(1.25, 0.5, 0),
   ahr = ahr, ahr_clm= 'stp',  colnm = pbsname, popgen = rank,
   rank_clm='end', ylimo=c(-0.08,0.5), stat_name='pbs', pch=16,
   cex=0.15, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75, xaxt="n"
  )
 }

 plot_pgen(
  crs = cross, chrs=ch, stat = pi, map = 'mid', mgp = c(1.25, 0.5, 0),
  ahr = ahr, ahr_clm= 'stp', colnm = piname, popgen = rank,
  rank_clm = 'end', ylimo=c(-0.03,0.03), stat_name='delta pi', pch=16,
  cex=0.15, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75
 )

 ## spacer empty plot ## ## ## ## ##
 plot(0,type='n',axes=FALSE,ann=FALSE)
 ## ## ## ## ## ## ## ## ## ## ## ##

 ## CM
 plot(
  bin2, map,chr=ch, lodcolumn=1, ylim=c(0,8),
  xaxt="n", mgp = c(1.25, 0.5, 0)
 )
 abline(v=ab, col='red',lwd=0.5)

 plot_ef(
  crs = cross2, map = map, pr = pr2 , ahr = ahr,
  popgen = rank, chs=ch, main=NULL, model="bin", xaxt="n",
  cex.axis = 0.5, cex.lab= 0.75, cex.main = 0.75
 )

 effectscan(
  cross_grid, pheno.col=4, chr=ch, get.se=T, draw=TRUE,
  gap=25, mtick="line",add.legend=F, alternate.chrid=T, ylim=c(-0.5,0.5),
  xlim=c(0,len), main=NULL,ylab='Effect Est. (+/- 1 SE)', xaxs="i",, xaxt="n",
  cex.axis = 0.75, cex.lab= 0.75, cex.main = 0.75, mgp = c(1.5, 0.5, 0)
 )
 abline(v=ab,col='red',lwd=0.5)

 plot(segA,segB, ylim=c(0,4), xlim=c(0,max(segX)), pch=16, col='lightgrey',
  xlab="",ylab='-log10 pval',xaxs="i", mgp = c(1.25, 0.5, 0),
  cex=0.25, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75)

 points(segX,segY,pch=16,cex=0.25)

 abline(v=ab,col='red',lwd=0.5)
 dev.off()
}
################################################################################
################################################################################
