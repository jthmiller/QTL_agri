#!/bin/R
### first run combine pops for multi-pop cross objects
pop <- 'NBH'
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
library("ggridges")
library("plyr")
library("scales")
library("ggrepel")
library('qtl')
library('RColorBrewer')

mpath <- '/home/jmiller1/QTL_agri/data'
setwd(mpath)

#load(file.path(mpath,'supplemental_plot_env.rsave'))
#cross <- cross_NBH
cross <- cross_ELR

fl <- file.path(mpath,paste0(pop'.mapped.tsp.csv'))

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)
################################################################################
################################################################################

erp <- 0.0025
cross <- sim.geno(cross, step=1, error.prob=erp, off.end=5, map.function="kosambi", n.draws=100, stepwidth="fixed")
cross <- calc.genoprob(cross, step=1, error.prob=erp, off.end=5, map.function="kosambi", stepwidth="fixed")
cross <- jittermap(cross)
gt <- geno.table(cross)
gt.map <- pull.map(cross)

cross_grid <- reduce2grid(cross)
scan <- scanone(cross_grid, method = "em", model = "binary", pheno.col = 4)
mpl <- summary(pull.map(cross_grid))

cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno),5)
norm <- scanone(cross_grid, method = "imp", model = "normal", pheno.col = 5)
bin <- scanone(cross_grid, method = "em", model = "binary", pheno.col = 4)

col <- c("slateblue", "violetred", "green3")

cross2 <- convert2cross2(cross)
map <- insert_pseudomarkers(cross2$gmap, step=1)
pr <- calc_genoprob(cross2, map, error_prob=0.0025, cores=4)
bin <- scan1(pr, pheno=cross2$pheno[,'bin'] , model = "binary", cores = 8)

################################################################################
################################################################################
#bottom, left, top and right
#location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks
ch <- 13
len <- mpl[ch,'length']
ind <- which(gt$chr == ch)
segX <- gt.map[[ch]][rownames(gt)[ind]]
#segX <- c(as.numeric(gsub(".*:","",rownames(gt)[ind])))/1000000
segY <- c(-log10(gt[ind,'P.value']))
ahr <- ahr_elr
ab <- ahr[which(ahr$chr==ch),'pos1']
################################################################################
################################################################################

################################################################################
################################################################################
pdf(paste0("/home/jmiller1/public_html/NBH_all",ch,"segdist.pdf"), width=4.5,height=7)
mat<-matrix(c(1:7),7,1, byrow=T)

layout(mat, widths=1, heights= c(0.1, 0.1, 0.025, 0.15, 0.15, 0.1, 0.1))
par(mar=c(0.5,2,0,1)+0.1, oma = c(1, 1, 0, 1))

plot_pgen(crs = cross, chrs=ch, stat = pbs, map = 'mid' ,
 ahr = ahr, ahr_clm= 'stp',  colnm = pop, popgen = nbh.rank,
 rank_clm='end', ylimo=c(-0.08,0.5), stat_name='pbs', pch=16, mgp = c(1.25, 0.5, 0),
 cex=0.15, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75, xaxt="n")

plot_pgen(crs = cross, chrs=ch, stat = taj, map = 'mid',
 ahr = ahr, ahr_clm= 'stp', colnm = pop, popgen = nbh.rank,
 rank_clm = 'end', ylimo=c(-4,4), stat_name='taj', pch=16, mgp = c(1.25, 0.5, 0),
 cex=0.15, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75)

#plot(0,type='n',cex=0.25, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75, xaxt="n")
plot(0,type='n',axes=FALSE,ann=FALSE)
## CM
plot(nbh_bin, nbh_map,chr=ch, lodcolumn=1, ylim=c(0,10), xaxt="n", mgp = c(1.25, 0.5, 0),)
abline(v=nbh_ab, col='red')

plot_ef(crs = nbh, map = nbh_map, pr = nbh_pr , ahr = ahr,
 popgen = nbh.rank, chs=ch, main=NULL, model="bin", xaxt="n",
 cex.axis = 0.5, cex.lab= 0.75, cex.main = 0.75)

effectscan(cross_grid, pheno.col=4, chr=ch, get.se=T, draw=TRUE,
 gap=25, mtick="line",add.legend=F, alternate.chrid=T,ylim=c(-0.5,0.5),
 xlim=c(0,len), main=NULL,ylab='Effect Est. (+/- 1 SE)', xaxs="i",, xaxt="n",
 cex.axis = 0.75, cex.lab= 0.75, cex.main = 0.75, mgp = c(1.5, 0.5, 0),)
abline(v=nbh_ab,col='red',lwd=0.5)

plot(segX,segY,ylim=c(0,4),xlim=c(0,max(segX)), xlab="",
 ylab='-log10 pval', pch=16, xaxs="i", mgp = c(1.25, 0.5, 0),
 cex=0.25, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75)
abline(v=nbh_ab,col='red',lwd=0.5)
dev.off()
################################################################################
################################################################################
