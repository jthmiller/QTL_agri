#!/bin/R
### first run combine pops for multi-pop cross objects
### plot with ER and NBH only

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

load('08_phys_plots_pos.rsave')

png("/home/jmiller1/public_html/pfst_nbh_18.png", width = 1000)
plot_stat_midpo(pfst,ch=18,poplot=statcol)
dev.off()

png("/home/jmiller1/public_html/pfst_elr_18.png", width = 1000)
plot_stat_midpo(pfst,ch=18,poplot=statcol)
dev.off()

png("/home/jmiller1/public_html/pfst_18.png", width = 1000)
plot_stat_mid(pfst,ch=18,poplot=statcol)
dev.off()

png("/home/jmiller1/public_html/pfst_18.png", width = 1000)
plot_stat_mid(pfst,ch=18,poplot=statcol)
dev.off()

plot_stat_midpo(pfst_conv,ch=8,poplot=statcol)
dev.off()

png("/home/jmiller1/public_html/pfst_8.png", width = 1000)
plot_stat_mid(pfst_conv,ch=8,poplot=statcol)
dev.off()




png("/home/jmiller1/public_html/pfst_18.png", width = 3000)
plot_stat_midpo(pfst_conv,ch=18,poplot=statcol)
dev.off()


png("/home/jmiller1/public_html/pfst_18.png", width = 3000)
plot_stat_sep(pfst,ch=18,poplot=statcol)
dev.off()

png("/home/jmiller1/public_html/pfst_18.png", width = 3000)
plot_stat_sep(pfst,ch=18,poplot=statcol)
dev.off()

png("/home/jmiller1/public_html/pbs.png", width = 3000)
plot_stat(pbs,ch=2,poplot=popgen)
dev.off()

png("/home/jmiller1/public_html/taj.png", width = 3000)
plot_stat(taj,ch=2,poplot=popout)
dev.off()

png("/home/jmiller1/public_html/pfst_conv.png", width = 1000)
plot_stat_sep(pfst_conv,ch=18,poplot=statcol)
dev.off()

################################################
### get positions of genes
nbh.gens <- cnv.ahrs(cross.nbh, AHRdf = AHR.bed, EXP = F)
new.gens <- cnv.ahrs(cross.new, AHRdf = AHR.bed, EXP = F)
elr.gens <- cnv.ahrs(cross.elr, AHRdf = AHR.bed, EXP = F)
brp.gens <- cnv.ahrs(cross.brp, AHRdf = AHR.bed, EXP = F)

qtl.gens <- nbh.gens[which(nbh.gens$chr %in% c(1, 2, 5, 8, 10, 12, 13, 18, 24)),]
minor.gens <- nbh.gens[which(nbh.gens$chr %in% c(8, 13, 23, 24)), ]
incompat.gens <- nbh.gens[which(nbh.gens$chr %in% c(8, 13)), ]
qtl_pg <- c(2,8, 13, 18, 24)
ol.gens <- nbh.gens[which(nbh.gens$chr %in% qtl_pg), ]
################################################

################################################
### ggplot popgen locations
nbh.popgen <- read.table("outliersNBH.txt.ncbi.lifted", sep = "\t", header = T)
new.popgen <- read.table("outliersNYC.txt.ncbi.lifted", sep = "\t", header = T)
elr.popgen <- read.table("outliersER.txt.ncbi.lifted", sep = "\t", header = T)
brp.popgen <- read.table("outliersBP.txt.ncbi.lifted", sep = "\t", header = T)
################################################

################################################
### Use nbh coords but elr and new popgen
new.rank <- cnv.popgen(cross.nbh, new.popgen, top = 50)
nbh.rank <- cnv.popgen(cross.nbh, nbh.popgen, top = 50)
elr.rank <- cnv.popgen(cross.nbh, elr.popgen, top = 50)
brp.rank <- cnv.popgen(cross.nbh, brp.popgen, top = 50)

nbh.rank$pop <- "NBH"
new.rank$pop <- "NEW"
elr.rank$pop <- "ELR"
brp.rank$pop <- "BRP"

all.rank <- rbind(new.rank, nbh.rank, elr.rank, brp.rank)
all.rank$pop <- factor(all.rank$pop, levels = c("NBH", "BRP", "NEW", "ELR"))
qtl.rank <- all.rank[which(all.rank$chr %in% c(1,2,5,8,10,12,13,18,23,24)),]
minor.rank <- all.rank[which(all.rank$chr %in% c(8, 13, 23, 24)), ]
incompat.rank <- all.rank[which(all.rank$chr %in% c(8, 13)), ]

qtl_pg <- c(2,8, 13, 18, 24)
ol.rank <- all.rank[which(all.rank$chr %in% qtl_pg), ]
################################################

### GGriges plot
scan_NBH <- scan.norm.imp.NBH
scan_ELR <- scan.norm.imp.ELR
scan_NEW <- scan.norm.imp.NEW
scan_BRP <- scan.norm.imp.BRP

scan_NBH <- scan.bin.mr.NBH
scan_ELR <- scan.bin.mr.ELR
scan_NEW <- scan.bin.mr.NEW
scan_BRP <- scan.bin.mr.BRP

scan_NBH <- scan.bin.imp.NBH
scan_ELR <- scan.bin.imp.ELR
scan_NEW <- scan.bin.imp.NEW
scan_BRP <- scan.bin.imp.BRP

png("/home/jmiller1/public_html/nbh.png", width = 3000)
plot(scan.bin.imp.NBH)
dev.off()

png("/home/jmiller1/public_html/BRP.png", width = 3000)
plot(scan.bin.imp.BRP)
dev.off()


melted.nbh <- data.frame(pop = "NBH", chr = scan_NBH$chr, pos = scan_NBH$pos,
  lod = scan_NBH$lod)
melted.new <- data.frame(pop = "NEW", chr = scan_NEW$chr, pos = scan_NEW$pos,
  lod = scan_NEW$lod)
melted.elr <- data.frame(pop = "ELR", chr = scan_ELR$chr, pos = scan_ELR$pos,
  lod = scan_ELR$lod)
melted.brp <- data.frame(pop = "BRP", chr = scan_BRP$chr, pos = scan_BRP$pos,
  lod = scan_BRP$lod)

melted <- rbind(melted.nbh, melted.new, melted.elr, melted.brp)
melted$pop <- factor(melted$pop, levels = rev(c("NBH", "BRP", "NEW", "ELR")))

## Total CM length of NBH. Rescale to NBH

mxes <- sapply(1:24,get_mxes,Y=themelt.nbh)
mxes <- sapply(1:24,get_mxes,Y=themelt.elr)

nbh.rescale <- melso(themelt.nbh)
new.rescale <- melso(themelt.new)
brp.rescale <- melso(themelt.brp)
elr.rescale <- melso(themelt.elr)

allmelt <- rbind(themelt.elr, new.rescale, nbh.rescale, brp.rescale)
allmelt$pop <- factor(allmelt$pop, levels = c("NBH", "BRP", "NEW", "ELR"))
qtlmelt <- allmelt[which(allmelt$chr %in% c(1,2,5,8,10,12,13,18,19,23,24)),
  ]
qtlminor <- allmelt[which(allmelt$chr %in% c(8,13,19,23,24)), ]
incompat <- allmelt[which(allmelt$chr %in% c(8,13)), ]

qtl_pg <- c(2,8, 13, 18, 24)
ol.melt <- allmelt[which(allmelt$chr %in% qtl_pg), ]

##save.image('/home/jmiller1/public_html/QTL_plot.Rsave')
################################################
qtlmelt <- allmelt[which(allmelt$chr %in% c( 2, 13, 18)),]
qtl.rank <- qtl.rank[which(qtl.rank$rank <21),]

qtl.rank <- all.rank[which(all.rank$chr %in% c(2, 13, 18)),]
qtl.gens <- nbh.gens[which(nbh.gens$chr %in% c( 2, 13, 18)),]
qtl.gens <- qtl.gens[as.character(c(1,3,4,5,6,7,8,24,26,45,46,56,70,96,94,57)),]
qtl.gens <- qtl.gens[c(1,3,4,5,6,7),]
qtl.gens <- qtl.gens[-5,]
################################################

qtl_pg <- c(1,5, 10, 12, 23)
ol.rank <- all.rank[which(all.rank$chr %in% qtl_pg), ]
ol.melt <- allmelt[which(allmelt$chr %in% qtl_pg), ]
ol.gens <- nbh.gens[which(nbh.gens$chr %in% qtl_pg), ]

names(popcol)[2] <- 'BRP'

################################################
################################################
