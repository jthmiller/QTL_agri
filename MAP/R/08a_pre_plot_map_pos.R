#!/bin/R
### first run combine pops for multi-pop cross objects

pop <- 'ELR'

source("/home/jmiller1/QTL_agri/MAP/control_file.R")
library("ggridges")
library("plyr")
library("scales")
library("ggrepel")
library('qtl')
library('RColorBrewer')

mpath <- '/home/jmiller1/QTL_agri/data'
setwd(mpath)

## FUNCTIONS
plot_stat_sep <- function(Z,ch,poplot){

  ind <- which(Z[,1] == ch)

  pops <- names(poplot)

  ymx_mn <- c(
    quantile(as.matrix(Z[ind,pops]), probs = 0.00001, na.rm = T),
    quantile(as.matrix(Z[ind,pops]), probs = 0.99999, na.rm = T))

  x_mx_mn <- c(min(Z[ind,'mid_midpo'],na.rm=T),max(Z[ind,'mid_midpo'],na.rm=T))

  X <- Z[ind,'mid_midpo']

  Y <- as.list(Z[ind,pops])
  names(Y) <- pops

  par(mfrow=c(length(pops),1),mar = c(1, 1, 1, 1),oma = c(1.5, 1.5, 1.5, 1.5))

  sapply(pops,plot_pop_sep,X,Y,poplot,x_mx_mn,ymx_mn)

  axis(side=1)

}

plot_pop_sep <- function(stat,X,Y,poplot,x_mx_mn,ymx_mn){
 plot(x_mx_mn, ymx_mn, type="n",xaxs="i", yaxs="i",main=NULL,xaxt="n",bty='n')
 points(X, Y[[stat]], pch=20, col=poplot[stat])
}

plot_stat <- function(Z,ch,poplot){

  ind <- which(Z[,1] == ch)

  pops <- names(poplot)

  ymx_mn <- c(
    quantile(as.matrix(Z[ind,pops]), probs = 0.00001, na.rm = T),
    quantile(as.matrix(Z[ind,pops]), probs = 0.99999, na.rm = T))

  x_mx_mn <- c(min(Z[ind,'mid_midpo'],na.rm=T),max(Z[ind,'mid_midpo'],na.rm=T))

  X <- Z[ind,'mid_midpo']

  Y <- as.list(Z[ind,pops])
  names(Y) <- pops

  plot(x_mx_mn, ymx_mn, type="n")
  sapply(pops,plot_pnts,X,Y,poplot)

}

plot_pnts <- function(stat,X,Y,poplot){ points(X, Y[[stat]], pch=20, col=poplot[stat]) }

##keeping colors consistent####################
all.pops <- c("NBH", "BRP", "ELR", "NEW")
popcol <- brewer.pal(8, "Paired")[c(2, 4, 6, 8)]
names(popcol) <- all.pops

popgen <- popcol
names(popgen) <- c('NBH','BP','ER','NYC')

popout <- c(popgen,'grey')
names(popout) <- c('NBH','BP','ER','NYC','BI')

### Color for stat comparisons
statcol <- popcol
names(statcol) <- c('BI.NBH','ER.KC','BP.F','NYC.SH')
################################################

pbs <- file.path(mpath, 'pbstat.txt.ncbi.lifted')
pbs <- read.table(pbs, sep = "\t", header = T)
pbs$mid <- pbs$V2 + (abs(pbs$V3 - pbs$V2) * .5)
pbs$V1 <- gsub('chr',"",pbs$V1)

pfst <- file.path(mpath, 'pfst.txt.ncbi.lifted')
pfst <- read.table(pfst, sep = "\t", header = T)
pfst$mid <- pfst$start + (abs(pfst$end - pfst$start) * .5)
pfst$Scaffold <- gsub('chr',"",pfst$Scaffold)

taj <- file.path(mpath, 'tajstat.txt.ncbi.lifted')
taj <- read.table(taj, sep = "\t", header = T)
taj$mid <- taj$start + (abs(taj$end - taj$start) * .5)
taj$Scaffold <- gsub('chr',"",taj$Scaffold)

pi <- file.path(mpath, 'piper.txt.ncbi.lifted')
pi <- read.table(pi, sep = "\t", header = T)
pi$mid <- pi$start + (abs(pi$end - pi$start) * .5)
pi$Scaffold <- gsub('chr',"",pi$Scaffold)

#### AHRs #####
AHR.bed <- read.table("lift_AHR_genes.bed", stringsAsFactors = F, header = F)
colnames(AHR.bed) <- c("chrom", "str", "stp", "gene")
AHR.bed$chrom <- as.numeric(gsub("chr", "", AHR.bed$chrom))
AHR.bed$str <- as.numeric(AHR.bed$str)
AHR.bed$stp <- as.numeric(AHR.bed$stp)
AHR.notmap <- AHR.bed[is.na(AHR.bed$chrom), ]
AHR.bed <- AHR.bed[!is.na(AHR.bed$chrom), ]
AHR.bed$gene <- gsub(":158640", "", AHR.bed$gene)
# add arnts (forgot to scan for them)
################################################

## Phenotypes
################################################
cross.BRP <- read.cross(format = "csv", dir = mpath, file = 'BRP.mapped.tsp.csv', genotypes=c("1","2","3"), estimate.map = FALSE)
cross.ELR <- read.cross(format = "csv", dir = mpath, file = 'ELR.mapped.tsp.csv', genotypes=c("1","2","3"), estimate.map = FALSE)
cross.NBH <- read.cross(format = "csv", dir = mpath, file = 'NBH.mapped.tsp.csv', genotypes=c("1","2","3"), estimate.map = FALSE)
cross.NEW <- read.cross(format = "csv", dir = mpath, file = 'NEW.mapped.tsp.csv', genotypes=c("1","2","3"), estimate.map = FALSE)
################################################

################################################
cor_nbh <- get_cor(cross.NBH)
cor_elr <- get_cor(cross.ELR)
cor_brp <- get_cor(cross.BRP)
cor_new <- get_cor(cross.NEW)

cross.BRP <- flip.order(cross.BRP, names(cor_brp)[which(cor_brp < 0)])
cross.NBH <- flip.order(cross.NBH, names(cor_nbh)[which(cor_nbh < 0)])
cross.NEW <- flip.order(cross.NEW, names(cor_new)[which(cor_new < 0)])
cross.ELR <- flip.order(cross.ELR, names(cor_elr)[which(cor_elr < 0)])
################################################

################################################
cross.nbh <- sim.geno(cross.NBH, n.draws = 500, step = 5, off.end = 10, error.prob = 0.025,
  map.function = "kosambi", stepwidth = "fixed")
cross.new <- sim.geno(cross.NEW, n.draws = 500, step = 5, off.end = 10, error.prob = 0.025,
  map.function = "kosambi", stepwidth = "fixed")
cross.elr <- sim.geno(cross.ELR, n.draws = 500, step = 5, off.end = 10, error.prob = 0.025,
  map.function = "kosambi", stepwidth = "fixed")
cross.brp <- sim.geno(cross.BRP, n.draws = 500, step = 5, off.end = 10, error.prob = 0.025,
  map.function = "kosambi", stepwidth = "fixed")
################################################

################################################
cross.nbh <- reduce2grid(cross.nbh)
cross.new <- reduce2grid(cross.new)
cross.elr <- reduce2grid(cross.elr)
cross.brp <- reduce2grid(cross.brp)
################################################

################################################
scan.norm.imp.NBH <- scanone(cross.nbh, method = "imp", model = "normal", pheno.col = 5)
scan.bin.imp.NBH <-  scanone(cross.nbh, method = "em", model = "binary", pheno.col = 4)
scan.norm.imp.ELR <- scanone(cross.elr, method = "imp", model = "normal", pheno.col = 5)
scan.bin.imp.ELR <-  scanone(cross.elr, method = "em", model = "binary", pheno.col = 4)
scan.norm.imp.NEW <- scanone(cross.new, method = "imp", model = "normal", pheno.col = 5)
scan.bin.imp.NEW <-  scanone(cross.new, method = "em", model = "binary", pheno.col = 4)
scan.norm.imp.BRP <- scanone(cross.brp, method = "imp", model = "normal", pheno.col = 5)
scan.bin.imp.BRP <-  scanone(cross.brp, method = "em", model = "binary", pheno.col = 4)
################################################
### marker regression plots
scan.bin.mr.NBH <-  scanone(cross.nbh, method = "mr", model = "binary", pheno.col = 4)
scan.bin.mr.ELR <-  scanone(cross.elr, method = "mr", model = "binary", pheno.col = 4)
scan.bin.mr.BRP <-  scanone(cross.brp, method = "mr", model = "binary", pheno.col = 4)
scan.bin.mr.NEW <-  scanone(cross.new, method = "mr", model = "binary", pheno.col = 4)
################################################

### use scanone for plots
themelt.nbh <- scan.bin.imp.NBH
themelt.new <- scan.bin.imp.NEW
themelt.elr <- scan.bin.imp.ELR
themelt.brp <- scan.bin.imp.BRP

themelt.nbh$pop <- "NBH"
themelt.new$pop <- "NEW"
themelt.elr$pop <- "ELR"
themelt.brp$pop <- "BRP"

themelt.brp.mr <- scan.bin.mr.BRP
themelt.elr.mr <- scan.bin.mr.ELR
themelt.nbh.mr <- scan.bin.mr.NBH
themelt.new.mr <- scan.bin.mr.NEW

themelt.brp.mr$pop <- 'BRP'
themelt.elr.mr$pop <- 'ELR'
themelt.nbh.mr$pop <- 'NBH'
themelt.new.mr$pop <- 'NEW'

save.image('08_phys_plots_pos.rsave')
################################################
################################################
### FROM MAP MAPPING
plot_stat <- function(Z,ch,poplot){

  ind <- which(Z[,1] == ch)

  pops <- names(poplot)

  ymx_mn <- c(
    quantile(as.matrix(Z[ind,pops]), probs = 0.00001, na.rm = T),
    quantile(as.matrix(Z[ind,pops]), probs = 0.99999, na.rm = T))

  x_mx_mn <- c(min(Z[ind,'mid'],na.rm=T),max(Z[ind,'mid'],na.rm=T))

  X <- Z[ind,'mid']

  Y <- as.list(Z[ind,pops])
  names(Y) <- pops

  plot(x_mx_mn, ymx_mn, type="n")
  sapply(pops,plot_pnts,X,Y,poplot)

}
plot_pnts <- function(stat,X,Y,poplot){ points(X, Y[[stat]], pch=20, col=poplot[stat]) }
plot_stat_sep <- function(Z,ch,poplot){

  ind <- which(Z[,1] == ch)

  pops <- names(poplot)

  ymx_mn <- c(
    quantile(as.matrix(Z[ind,pops]), probs = 0.00001, na.rm = T),
    quantile(as.matrix(Z[ind,pops]), probs = 0.99999, na.rm = T))

  x_mx_mn <- c(min(Z[ind,'mid'],na.rm=T),max(Z[ind,'mid'],na.rm=T))

  X <- Z[ind,'mid']

  Y <- as.list(Z[ind,pops])
  names(Y) <- pops

  par(mfrow=c(length(pops),1),mar = c(1, 1, 1, 1),oma = c(1.5, 1.5, 1.5, 1.5))

  sapply(pops,plot_pop_sep,X,Y,poplot,x_mx_mn,ymx_mn)

  axis(side=1)

}
plot_pop_sep <- function(stat,X,Y,poplot,x_mx_mn,ymx_mn){
 plot(x_mx_mn, ymx_mn, type="n",xaxs="i", yaxs="i",main=NULL,xaxt="n",bty='n')
 points(X, Y[[stat]], pch=20, col=poplot[stat])
}

### FROM PHYS MAPPING
## FUNCTIONS
plot_stat_sep <- function(Z,ch,poplot){

  ind <- which(Z[,1] == ch)

  pops <- names(poplot)

  ymx_mn <- c(
    quantile(as.matrix(Z[ind,pops]), probs = 0.00001, na.rm = T),
    quantile(as.matrix(Z[ind,pops]), probs = 0.99999, na.rm = T))

  x_mx_mn <- c(min(Z[ind,'mid_midpo'],na.rm=T),max(Z[ind,'mid_midpo'],na.rm=T))

  X <- Z[ind,'mid_midpo']

  Y <- as.list(Z[ind,pops])
  names(Y) <- pops

  par(mfrow=c(length(pops),1),mar = c(1, 1, 1, 1),oma = c(1.5, 1.5, 1.5, 1.5))

  sapply(pops,plot_pop_sep,X,Y,poplot,x_mx_mn,ymx_mn)

  axis(side=1)

}
plot_pop_sep <- function(stat,X,Y,poplot,x_mx_mn,ymx_mn){
 plot(x_mx_mn, ymx_mn, type="n",xaxs="i", yaxs="i",main=NULL,xaxt="n",bty='n')
 points(X, Y[[stat]], pch=20, col=poplot[stat])
}
plot_stat <- function(Z,ch,poplot){

  ind <- which(Z[,1] == ch)

  pops <- names(poplot)

  ymx_mn <- c(
    quantile(as.matrix(Z[ind,pops]), probs = 0.00001, na.rm = T),
    quantile(as.matrix(Z[ind,pops]), probs = 0.99999, na.rm = T))

  x_mx_mn <- c(min(Z[ind,'mid_midpo'],na.rm=T),max(Z[ind,'mid_midpo'],na.rm=T))

  X <- Z[ind,'mid_midpo']

  Y <- as.list(Z[ind,pops])
  names(Y) <- pops

  plot(x_mx_mn, ymx_mn, type="n")
  sapply(pops,plot_pnts,X,Y,poplot)

}
plot_pnts <- function(stat,X,Y,poplot){ points(X, Y[[stat]], pch=20, col=poplot[stat]) }
