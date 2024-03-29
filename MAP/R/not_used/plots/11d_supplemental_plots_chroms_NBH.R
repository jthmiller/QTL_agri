#!/bin/R
### first run combine pops for multi-pop cross objects
pop <- 'NBH'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
library("ggridges")
library("plyr")
library("scales")
library("ggrepel")
library('qtl')
library('RColorBrewer')

mpath <- '/home/jmiller1/QTL_agri/data'
setwd(mpath)

cores <- detectCores() - 2

fl <- file.path(mpath,paste0(pop'.mapped.tsp.csv'))

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno),5)
################################################################################
################################################################################

#load(file.path(mpath,'supplemental_plot_env.rsave'))
#cross <- cross_NBH
#cross <- cross_ELR

################################################################################
################################################################################

norm <- scanone(cross, pheno.col = 5, model='normal', method = "imp")
bin <- scanone(cross, pheno.col = 4, method = "em", model = "binary")
gt <- geno.table(cross)
map <- pull.map(cross)
map_sum <- summary(pull.map(cross))

################################################################################
cross_grid <- subset(cross,chr=c(1:4,6:24))
cross_grid <- reduce2grid(cross_grid)

bin_grid <- scanone(cross, method = "em", model = "binary", pheno.col = 4)
norm_grid <- scanone(cross, method = "imp", model = "normal", pheno.col = 5)
gt_grid <- geno.table(cross)
map_grid <- pull.map(cross)
map_grid_sum <- summary(pull.map(cross))

################################################################################
################################################################################

col <- c("slateblue", "violetred", "green3")
cross2 <- convert2cross2(cross)
map <- insert_pseudomarkers(cross2$gmap, step=0.5)
pr <- calc_genoprob(cross2, map2, error_prob=0.0025, cores=cores)
bin2<- scan1(pr2, pheno=cross2$pheno[,'bin'] , model = "binary", cores = cores)

################################################################################
################################################################################
#bottom, left, top and right
#location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks
ch <- 2
len <- map_sum[ch,'length']
ind <- which(gt$chr == ch)
segX <- map[[ch]][rownames(gt)[ind]]
#segX <- c(as.numeric(gsub(".*:","",rownames(gt)[ind])))/1000000
segY <- c(-log10(gt[ind,'P.value']))

if(pop == 'NBH') { rank <- nbh.rank ; ahr <- ahr_nbh ; statname <- popname <- 'NBH' }
if(pop == 'ELR') { rank <- elr.rank ; ahr <- ahr_elr ; statname <- popname <- 'ER' }
ab <- ahr[which(ahr$chr==ch),'pos1']
################################################################################
################################################################################
stat <- 'pfst' ; statname <- 'BI.NBH'
stat <- 'pbst'
stat <- 'pfst' ; statname <- 'ER.KC'
################################################################################
################################################################################
pdf(paste0("/home/jmiller1/public_html/",pop,"_all",ch,"segdist2.pdf"), width=4.5,height=6)
mat < -matrix(c(1:7),7,1, byrow=T)

layout(mat, widths=1, heights= c(0.1, 0.1, 0.025, 0.1, 0.1, 0.1, 0.05))
par(mar=c(0.5,2.25,0,1)+0.1, oma = c(1, 1, 0, 1))

if(stat == 'pfst'){
 plot_pgen(
  crs = cross, chrs=ch, stat = pfst, map = 'mid', mgp = c(1.25, 0.5, 0),
  ahr = ahr, ahr_clm= 'stp',  colnm = statname, popgen = rank,
  rank_clm='end', ylimo=c(-0.01,0.6), stat_name='pfst NxS', pch=16,
  cex=0.15, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75, xaxt="n"
 )
} else {
 plot_pgen(
  crs = cross, chrs=ch, stat = pbs, map = 'mid', mgp = c(1.25, 0.5, 0),
  ahr = ahr, ahr_clm= 'stp',  colnm = popname, popgen = rank,
  rank_clm='end', ylimo=c(-0.08,0.5), stat_name='pbs', pch=16,
  cex=0.15, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75, xaxt="n"
 )
}

plot_pgen(
 crs = cross, chrs=ch, stat = pi, map = 'mid', mgp = c(1.25, 0.5, 0),
 ahr = ahr, ahr_clm= 'stp', colnm = statname, popgen = rank,
 rank_clm = 'end', ylimo=c(-0.02,0.015), stat_name='delta pi', pch=16,
 cex=0.15, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75
)

## spacer empty plot ## ## ## ## ##
plot(0,type='n',axes=FALSE,ann=FALSE)
## ## ## ## ## ## ## ## ## ## ## ##

## CM
plot(
 bin2, map,chr=ch, lodcolumn=1, ylim=c(0,10),
 xaxt="n", mgp = c(1.25, 0.5, 0)
)
abline(v=ab, col='red',lwd=0.5)

plot_ef(
 crs = cross2, map = map, pr = pr , ahr = ahr,
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

plot(
 segX,segY,ylim=c(0,4),xlim=c(0,max(segX)), xlab="",
 ylab='-log10 pval', pch=16, xaxs="i", mgp = c(1.25, 0.5, 0),
 cex=0.25, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75
)
abline(v=ab,col='red',lwd=0.5)
dev.off()
################################################################################
################################################################################
