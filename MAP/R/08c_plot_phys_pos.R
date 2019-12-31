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

#############################################
### pfst
statcol_mod <- popcol
names(statcol_mod) <- c('BI.NBH','F.NBH','ER.KC','ER.SH')
statcol_mod[c(2,3)] <- c('darkblue','red4')
ab1 <- which.min(abs(pfst[which(pfst$chr==13),'mid'] - 1009168))
ab2 <- which.min(abs(pfst[which(pfst$chr==13),'mid'] - 3595485))
ab1 <- pfst[which(pfst$chr==13),][ab1,'mid']
ab2 <- pfst[which(pfst$chr==13),][ab2,'mid']

for (ch in c(1,2,8,13,9,1,18,24)){
 png(paste0("/home/jmiller1/public_html/pfst",ch,".png"), width = 500)
 plot_stat(pfst,ch=ch,poplot=statcol_mod,colnm='mid')

if (ch==13){
 abline(v=ab1)
 abline(v=ab2)
}

dev.off()
}
###################################################
#############################################
## pbs

for (ch in c(1,2,8,13,9,1,18,24)){
 png(paste0("/home/jmiller1/public_html/pbs",ch,".png"), width = 500)
 plot_stat(pbs,ch=ch,poplot=popgen,colnm='mid')

 if (ch==13){
  abline(v=ab1)
  abline(v=ab2)
 }

 dev.off()
}

################################################]
## tajimas d
for (ch in c(1,2,8,13,9,1,18,24)){
 png(paste0("/home/jmiller1/public_html/taj",ch,".png"), width = 500)
 plot_stat(taj,ch=ch,poplot=popgen,colnm='mid')

 if (ch==13){
  abline(v=ab1)
  abline(v=ab2)
 }

 dev.off()
}
################################################
## pi
for (ch in c(1,2,8,13,9,1,18,24)){
 png(paste0("/home/jmiller1/public_html/pi",ch,".png"), width = 500)
 plot_stat(pi,ch=ch,poplot=popgen,colnm='mid')

 if (ch==13){
  abline(v=ab1)
  abline(v=ab2)
 }

 dev.off()
}


load('08_phys_plots_pos.rsave')
################################################

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

##TRY WITH ELR AND NBH COORDS

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

################################################
### GGriges plot
scan_NBH <- scan.norm.imp.NBH
scan_NBH <- scan.norm.imp.NBH
scan_NBH <- scan.norm.imp.NBH
scan_NBH <- scan.norm.imp.NBH

scan_NBH <- scan.bin.mr.NBH
scan_ELR <- scan.bin.mr.ELR
scan_NEW <- scan.bin.mr.NEW
scan_BRP <- scan.bin.mr.BRP

melted.nbh <- data.frame(pop = "NBH", chr = scan_NBH$chr, pos = scan_NBH$pos,
  lod = scan_NBH$lod)
melted.new <- data.frame(pop = "NEW", chr = scan_NEW$chr, pos = scan_NEW$pos,
  lod = scan_NEW$lod)
melted.elr <- data.frame(pop = "ELR", chr = scan_ELR$chr, pos = scan_ELR$pos,
  lod = scan_ELR$lod)
melted.brp <- data.frame(pop = "BRP", chr = scan_BRP$chr, pos = scan_BRP$pos,
  lod = scan_BRP$lod)


nbh.pos <- data.frame(pop = "NBH", chr = scan_NBH$chr, pos = rownames(scan_NBH), lod = scan_NBH$lod)
elr.pos <- data.frame(pop = "ELR", chr = scan_ELR$chr, pos = rownames(scan_ELR), lod = scan_ELR$lod)
new.pos <- data.frame(pop = "NEW", chr = scan_NEW$chr, pos = rownames(scan_NEW), lod = scan_NEW$lod)
brp.pos <- data.frame(pop = "BRP", chr = scan_BRP$chr, pos = rownames(scan_BRP), lod = scan_BRP$lod)

nbh.pos$pos <- as.numeric(gsub("*.:","",nbh.pos$pos))
elr.pos$pos <- as.numeric(gsub("*.:","",elr.pos$pos))
new.pos$pos <- as.numeric(gsub("*.:","",new.pos$pos))
brp.pos$pos <- as.numeric(gsub("*.:","",brp.pos$pos))

maxmin <- function(ch){
 a <- rbind(nbh.pos,elr.pos,new.pos,brp.pos)
 max(a[which(a$chr == ch),'pos'],na.rm=T)
}

ydir <- function(data,chr){
 data[which(data$chr == chr),'lod']
}

xdir <- function(data,chr){
 data[which(data$chr == chr),'pos']
}






plot_R_stat <- function(ch){

 X <- c(0,maxmin(ch))
 Y <- c(0,20)
 plot(X,Y, type="n")
 points(xdir(nbh.pos,ch), ydir(nbh.pos,ch), pch=20, col=popcol['NBH'])
 points(xdir(elr.pos,ch), ydir(elr.pos,ch), pch=20, col=popcol['ELR'])
 points(xdir(new.pos,ch), ydir(new.pos,ch), pch=20, col=popcol['NEW'])
 points(xdir(brp.pos,ch), ydir(brp.pos,ch), pch=20, col=popcol['BRP'])

}

png("/home/jmiller1/public_html/log_mr.png", width = 1000)
plot_R_stat(ch=18)
dev.off()


#### base R ##########
plot_R_stat <- function(Z,ch,poplot){

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






##### ggplot ########
melted <- rbind(melted.nbh, melted.new, melted.elr, melted.brp)
melted$pop <- factor(melted$pop, levels = rev(c("NBH", "BRP", "NEW", "ELR")))

## Total CM length of NBH. Rescale to NBH
mxes <- sapply(1:24, function(X) {
  max(themelt.nbh.mr$pos[which(themelt.nbh.mr$chr == X)])
})

melso <- function(tomelt){
 ts <- tomelt[which(tomelt$chr == 1), ]
 ts$pos <- rescale(ts$pos, to = c(-10, mxes[1]))
 the_rescale <- ts
 for (i in 2:24) {
   ts <- tomelt[which(tomelt$chr == i), ]
   ts$pos <- rescale(ts$pos, to = c(-10, mxes[i]))
   the_rescale <- rbind(the_rescale, ts)
  }
 the_rescale
}


elr.rescale.mr <- melso(themelt.elr.mr)
brp.rescale.mr <- melso(themelt.brp.mr)
new.rescale.mr <- melso(themelt.new.mr)

new.rescale <- melso(themelt.new)
brp.rescale <- melso(themelt.brp)
elr.rescale <- melso(themelt.elr)

##ts <- themelt.new[which(themelt.new$chr == 1), ]
##ts$pos <- rescale(ts$pos, to = c(-10, mxes[1]))
##new.rescale <- ts
##for (i in 2:24) {
##  ts <- themelt.new[which(themelt.new$chr == i), ]
##  ts$pos <- rescale(ts$pos, to = c(-10, mxes[i]))
##  new.rescale <- rbind(new.rescale, ts)
##}

##ts <- themelt.elr[which(themelt.elr$chr == 1), ]
##ts$pos <- rescale(ts$pos, to = c(-10, mxes[1]))
##elr.rescale <- ts
##for (i in 2:24) {
##  ts <- themelt.elr[which(themelt.elr$chr == i), ]
##  ts$pos <- rescale(ts$pos, to = c(-10, mxes[i]))
##  elr.rescale <- rbind(elr.rescale, ts)
##}
##
##ts <- themelt.brp[which(themelt.brp$chr == 1), ]
##ts$pos <- rescale(ts$pos, to = c(-10, mxes[1]))
##brp.rescale <- ts
##for (i in 2:24) {
##  ts <- themelt.brp[which(themelt.brp$chr == i), ]
##  ts$pos <- rescale(ts$pos, to = c(-10, mxes[i]))
##  brp.rescale <- rbind(brp.rescale, ts)
##}

allmelt <- rbind(themelt.nbh.mr, new.rescale.mr, elr.rescale.mr, brp.rescale.mr)
allmelt$pop <- factor(allmelt$pop, levels = c("NBH", "BRP", "NEW", "ELR"))
qtlmelt <- allmelt[which(allmelt$chr %in% c(1,2,5,8,10,12,13,18,19,23,24)),
  ]
qtlminor <- allmelt[which(allmelt$chr %in% c(8,13,19,23,24)), ]
incompat <- allmelt[which(allmelt$chr %in% c(8,13)), ]

qtl_pg <- c(2,8, 13, 18, 24)
ol.melt <- allmelt[which(allmelt$chr %in% qtl_pg), ]


################################################
p <- ggplot(themelt.nbh.mr, aes(x = pos, y = lod))
png("/home/jmiller1/public_html/all_popgen_rank_scaled.qtl.png", width = 3000)
p + facet_wrap(~chr, scales = "free_x", nrow = 1, ncol = 24) + scale_y_continuous(limits = c(-12,
  23)) + geom_line(size = 2, alpha = 0.6) + theme_minimal() + theme(axis.text = element_text(size = 10)) +
  labs(x = "Chromosome", y = "LOD", linetype = "") + geom_label_repel(aes(x = pos,
  y = -0.1, label = gene), box.padding = unit(0.25, "lines"), parse = T, point.padding = unit(0.2,
  "lines"), force = 10, label.padding = unit(0.2, "lines"), ylim = c(0, -12), segment.size = 1,
  max.iter = 6000, data = nbh.gens, direction = "y", size = 4, seed = 666, nudge_y = -0.01,
  vjust = 3) + geom_label_repel(aes(x = pos, y = -0.1, label = gene), box.padding = unit(0.25,
  "lines"), parse = T, point.padding = unit(0.2, "lines"), force = 10, label.padding = unit(0.2,
  "lines"), ylim = c(0, -12), segment.size = 0, max.iter = 6000, data = nbh.gens,
  direction = "y", size = 4, seed = 666, nudge_y = -0.01, vjust = 3) + geom_point(aes(size = rank),
  data = nbh.rank) + geom_label_repel(aes(x = pos, y = lod, label = rank), data = nbh.rank,
  size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines"),
  vjust = 1)
dev.off()

p <- ggplot(themelt.new, aes(x = pos, y = lod))
png("/home/jmiller1/public_html/new_genes_below_rank_above.qtl.png", width = 3000)
p + facet_wrap(~chr, scales = "free_x", nrow = 1, ncol = 24) + scale_y_continuous(limits = c(-12,
  23)) + geom_line(size = 2, alpha = 0.6) + theme_minimal() + theme(axis.text = element_text(size = 10)) +
  labs(x = "Chromosome", y = "LOD", linetype = "") + geom_label_repel(aes(x = pos,
  y = -0.1, label = gene), box.padding = unit(0.25, "lines"), parse = T, point.padding = unit(0.2,
  "lines"), force = 10, label.padding = unit(0.2, "lines"), ylim = c(0, -12), segment.size = 1,
  max.iter = 6000, data = new.gens, direction = "y", size = 4, seed = 666, nudge_y = -0.01,
  vjust = 3) + geom_label_repel(aes(x = pos, y = -0.1, label = gene), box.padding = unit(0.25,
  "lines"), parse = T, point.padding = unit(0.2, "lines"), force = 10, label.padding = unit(0.2,
  "lines"), ylim = c(0, -12), segment.size = 0, max.iter = 6000, data = new.gens,
  direction = "y", size = 4, seed = 666, nudge_y = -0.01, vjust = 3) + geom_label_repel(aes(x = pos,
  y = lod, label = rank), data = nbh.rank, size = 4, box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"), vjust = 1)
dev.off()

p <- ggplot(themelt.elr, aes(x = pos, y = lod, color = pop))
png("/home/jmiller1/public_html/rank_top_scaled_random.qtl.png", width = 3000)
p + facet_wrap(~chr, scales = "free_x", nrow = 1, ncol = 24) + scale_y_continuous(limits = c(-12,
  23)) + scale_color_manual(values = popcol) + geom_line(size = 2, alpha = 0.6) +
  theme_minimal() + theme(axis.text = element_text(size = 10)) + labs(x = "Chromosome",
  y = "LOD", linetype = "") + geom_label_repel(aes(x = pos, y = -0.1, label = gene),
  box.padding = unit(0.25, "lines"), parse = T, point.padding = unit(0.2, "lines"),
  force = 10, label.padding = unit(0.2, "lines"), ylim = c(0, -12), segment.size = 1,
  max.iter = 6000, data = elr.gens, direction = "y", size = 4, seed = 666, nudge_y = -0.01,
  vjust = 3) + geom_label_repel(aes(x = pos, y = -0.1, label = gene), box.padding = unit(0.25,
  "lines"), parse = T, point.padding = unit(0.2, "lines"), force = 10, label.padding = unit(0.2,
  "lines"), ylim = c(0, -12), segment.size = 0, max.iter = 6000, data = elr.gens,
  direction = "y", size = 4, seed = 666, nudge_y = -0.01, vjust = 3) + geom_label_repel(aes(x = pos,
  y = lod, label = rank), data = elr.rank, size = 4, box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines"), vjust = 1)
dev.off()

p <- ggplot(allmelt, aes(x = pos, y = lod, colour = pop))
png("/home/jmiller1/public_html/no_annot_.qtl.png", width = 2000)
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 24) + scale_y_continuous(limits = c(0,
  22)) + scale_color_manual(values = popcol) + geom_line(size = 2, alpha = 0.85) +
  theme_minimal() + theme(axis.text = element_text(size = 10)) + labs(x = "Chromosome",
  y = "LOD", linetype = "")
dev.off()

p <- ggplot(allmelt, aes(x = pos, y = lod, colour = pop))
png("/home/jmiller1/public_html/ahr_genes_only.qtl.png", width = 2000)
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 24) + scale_y_continuous(limits = c(0,
  22)) + scale_color_manual(values = popcol) + geom_line(size = 2, alpha = 0.6) +
  theme_minimal() + theme(axis.text = element_text(size = 10)) + geom_label_repel(aes(x = pos,
  y = 0.1, label = gene), box.padding = unit(0.5, "lines"), parse = T, point.padding = unit(0.8,
  "lines"), force = 1, label.padding = unit(0.2, "lines"), ylim = c(5, 20), segment.size = 1,
  max.iter = 6000, data = nbh.gens, direction = "y", size = 3, seed = 666, nudge_y = 0.01,
  vjust = 0.1) + geom_label_repel(aes(x = pos, y = 0.1, label = gene), box.padding = unit(0.5,
  "lines"), parse = T, point.padding = unit(0.8, "lines"), force = 1, label.padding = unit(0.2,
  "lines"), ylim = c(5, 20), segment.size = 0, max.iter = 6000, data = nbh.gens,
  direction = "y", size = 3, seed = 666, nudge_y = 0.01, vjust = 0.1) + labs(x = "Chromosome",
  y = "LOD", linetype = "")
dev.off()

p <- ggplot(allmelt, aes(x = pos, y = lod, colour = pop))
png("/home/jmiller1/public_html/All_chr.qtl.png", width = 3000)
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 24) + scale_color_manual(values = popcol) +
  scale_y_continuous(limits = c(-5, 22)) + geom_label_repel(aes(x = pos, y = 0.1,
  label = gene), color = "black", segment.size = 1, label.padding = unit(0.2, "lines"),
  data = nbh.gens, direction = "y", size = 3, seed = 666, ylim = c(5, 20), segment.alpha = 0.5) +
  geom_line(size = 2, alpha = 0.5) + geom_label_repel(aes(x = pos, y = 0.1, label = gene),
  color = "black", label.padding = unit(0.2, "lines"), ylim = c(5, 20), segment.size = 0,
  max.iter = 6000, nudge_y = 2, data = nbh.gens, direction = "y", size = 3.5, seed = 666) +
  geom_label_repel(aes(x = pos, y = 0, fill = pop, label = rank), data = all.rank,
    size = 3.5, segment.size = 1, force = 4, min.segment.length = 0.1, point.padding = unit(0.4,
      "lines"), direction = "both", ylim = c(0, -8), seed = 666, box.padding = 0.1) +
  geom_label_repel(aes(x = pos, y = 0, label = rank, ), fontface = "bold", data = all.rank,
    size = 3.5, segment.size = 0, force = 4, min.segment.length = 0.1, point.padding = unit(0.4,
      "lines"), direction = "both", ylim = c(0, -8), seed = 666, box.padding = 0.1) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = "none") + labs(x = "Chromosome", y = "LOD", linetype = "")
dev.off()


qtlmelt <- allmelt[which(allmelt$chr %in% c( 2, 13, 18)),]
qtl.rank <- qtl.rank[which(qtl.rank$rank <21),]

qtl.rank <- all.rank[which(all.rank$chr %in% c(2, 13, 18)),]
qtl.gens <- nbh.gens[which(nbh.gens$chr %in% c( 2, 13, 18)),]
qtl.gens <- qtl.gens[as.character(c(1,3,4,5,6,7,8,24,26,45,46,56,70,96,94,57)),]
qtl.gens <- qtl.gens[c(1,3,4,5,6,7),]
qtl.gens <- qtl.gens[-5,]

p <- ggplot(qtlmelt, aes(x = pos, y = lod, colour = pop))
pdf("/home/jmiller1/public_html/all_pop_qtl_only.pdf",width=15,height=11)
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 10) +
  scale_color_manual(values = popcol) +
  scale_y_continuous(limits = c(-5, 22)) +
  ### gene line
  geom_label_repel(aes(x = pos, y = 0.1, label = gene),fill='white', color = "black", segment.size = 1,
    label.padding = unit(0.2, "lines"), data = qtl.gens, direction = "y", size = 5,
    seed = 666, nudge_y = 2, ylim = c(5, 20), segment.alpha = 0.5) + ## plot lines

  geom_line(size = 2, alpha = 0.75) + ### gene label

  geom_label_repel(aes(x = pos, y = 0.1, label = gene), color = "black", segment.size = 0,
    label.padding = unit(0.2, "lines"), data = qtl.gens, direction = "y", size = 5,
    seed = 666, nudge_y = 2, ylim = c(5, 20), max.iter = 6000, ) + ### Rank line

  geom_label_repel(aes(x = pos, y = 0, color = pop, label = rank), data = qtl.rank,
    size = 5, segment.size = 1, force = 4, min.segment.length = 0.1, point.padding = unit(0.4,
    "lines"), direction = "both", ylim = c(0, -8), seed = 666, box.padding = 0.1) +
    ### Rank label

  geom_label_repel(aes(x = pos, y = 0, color = pop,label = rank), fontface = "bold", data = qtl.rank,
    size = 5, segment.size = 0, force = 4, min.segment.length = 0.1, point.padding = unit(0.4,
    "lines"), direction = "both", ylim = c(0, -8), seed = 666, box.padding = 0.1) +

  theme(axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", color = "black", size = 16),
    axis.text.y=element_text(face = "bold", color = "black", size = 16),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(face = "bold", color = "black",size=16),
    legend.position = "none") +
  labs(x = "Chromosome", y = "LOD", linetype = "")
dev.off()
## Why would chromosomes with top ranked loci have fewest signatures?

p <- ggplot(qtlminor, aes(x = pos, y = lod, colour = pop))
png("/home/jmiller1/public_html/minor_qtl_only.png", width = 2000)
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 9) + scale_color_manual(values = popcol) +
  scale_y_continuous(limits = c(-1, 8)) + theme_minimal() + geom_label_repel(aes(x = pos,
  y = 0.1, label = gene), color = "black", segment.size = 1, label.padding = unit(0.2,
  "lines"), data = minor.gens, direction = "y", size = 5, seed = 666, nudge_y = 2,
  ylim = c(5, 20), segment.alpha = 0.5) + geom_line(size = 2, alpha = 0.5) + geom_label_repel(aes(x = pos,
  y = 0.1, label = gene), color = "black", label.padding = unit(0.2, "lines"),
  ylim = c(5, 20), segment.size = 0, max.iter = 6000, nudge_y = 2, data = minor.gens,
  direction = "y", size = 5, seed = 666) + geom_label_repel(aes(x = pos, y = 0,
  color = pop, label = rank), data = minor.rank, size = 5, segment.size = 1, force = 4,
  min.segment.length = 0.1, point.padding = unit(0.4, "lines"), direction = "both",
  ylim = c(0, -8), seed = 666, box.padding = 0.1) + geom_label_repel(aes(x = pos,
  y = 0, label = rank), fontface = "bold", data = minor.rank, size = 5, segment.size = 0,
  force = 4, min.segment.length = 0.1, point.padding = unit(0.4, "lines"), direction = "both",
  ylim = c(0, -8), seed = 666, box.padding = 0.1) + theme(axis.title.x = element_blank(),
  axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  labs(x = "Chromosome", y = "LOD", linetype = "")
dev.off()

p <- ggplot(incompat, aes(x = pos, y = lod, colour = pop))
png("/home/jmiller1/public_html/elr_incompat.png", width = 750)
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 9) +
  theme(strip.text.x = element_text(size = 12)) +
  scale_color_manual(values = popcol) +
  scale_y_continuous(limits = c(-1, 8)) +
  geom_label_repel(aes(x = pos, y = 0.1, label = gene), color = "black", segment.size = 1,
    label.padding = unit(0.2, "lines"), data = incompat.gens, direction = "y",
    size = 5, seed = 666, nudge_y = 2, ylim = c(5, 20), segment.alpha = 0.5) +
  geom_line(size = 2, alpha = 0.5) + geom_label_repel(aes(x = pos, y = 0.1, label = gene),
  color = "black", label.padding = unit(0.2, "lines"), ylim = c(5, 20), segment.size = 0,
  max.iter = 6000, nudge_y = 2, data = incompat.gens, direction = "y", size = 5,
  seed = 666) +
  geom_label_repel(aes(x = pos, y = 0, color = pop, label = rank),
  data = incompat.rank, size = 5, segment.size = 1, force = 4, min.segment.length = 0.1,
  point.padding = unit(0.4, "lines"), direction = "both", ylim = c(0, -8), seed = 666,
  box.padding = 0.1) +
  geom_label_repel(aes(x = pos, y = 0, label = rank), fontface = "bold",
  data = incompat.rank, size = 5, segment.size = 0, force = 4, min.segment.length = 0.1,
  point.padding = unit(0.4, "lines"), direction = "both", ylim = c(0, -8), seed = 666,
  box.padding = 0.1) + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
  axis.ticks.x = element_blank(), legend.position = "none") + labs(x = "Chromosome",
  y = "LOD", linetype = "")
dev.off()

### QTL only

p <- ggplot(ol.melt, aes(x = pos, y = lod, colour = pop))
png("/home/jmiller1/public_html/ol_qtl_only.png", width = 2000)
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 9) +
scale_color_manual(values = popcol) +
  scale_y_continuous(limits = c(-1, 8)) +
  theme_minimal() + geom_label_repel(aes(x = pos,
  y = 0.1, label = gene), color = "black", segment.size = 1, label.padding = unit(0.2,
  "lines"), data = minor.gens, direction = "y", size = 5, seed = 666, nudge_y = 2,
  ylim = c(5, 20), segment.alpha = 0.5) +
  geom_line(size = 2, alpha = 0.5) +
  geom_label_repel(aes(x = pos,y = 0.1, label = gene),
  color = "black", label.padding = unit(0.2, "lines"),
  ylim = c(5, 20), segment.size = 0, max.iter = 6000, nudge_y = 2, data = minor.gens,
  direction = "y", size = 5, seed = 666) +
  geom_label_repel(aes(x = pos, y = 0,color = pop, label = rank),
  data = minor.rank, size = 5, segment.size = 1, force = 4,
  min.segment.length = 0.1, point.padding = unit(0.4, "lines"), direction = "both",
  ylim = c(0, -8), seed = 666, box.padding = 0.1) +
  geom_label_repel(aes(x = pos,y = 0, label = rank),
  fontface = "bold", data = minor.rank, size = 5, segment.size = 0,
  force = 4, min.segment.length = 0.1, point.padding = unit(0.4, "lines"), direction = "both",
  ylim = c(0, -8), seed = 666, box.padding = 0.1) +
  theme(axis.title.x = element_blank(),
  axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  labs(x = "Chromosome", y = "LOD", linetype = "")
dev.off()


names(popcol)[2] <- 'BRP'

p <- ggplot(ol.melt, aes(x = pos, y = lod, colour = pop))
pdf("/home/jmiller1/public_html/OL_chr.qtl.pdf")
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 24) +
  scale_color_manual(values = popcol) +
  scale_y_continuous(limits = c(-5, 22)) +
  geom_label_repel(aes(x = pos, y = 0.1,label = gene),data = ol.gens,
    color = "black", segment.size = 1, label.padding = unit(0.2, "lines"),
    direction = "y", size = 3, seed = 666, ylim = c(5, 20), segment.alpha = 0.5) +
  geom_line(size = 2, alpha = 0.5) +
  geom_label_repel(aes(x = pos, y = 0.1, label = gene), data = ol.gens,
    color = "black", label.padding = unit(0.2, "lines"), ylim = c(5, 20), segment.size = 0,
    max.iter = 6000, nudge_y = 2, direction = "y", size = 3.5, seed = 666) +
  geom_label_repel(aes(x = pos, y = 0, fill = pop, label = rank), data = ol.rank,
    size = 3.5, segment.size = 1, force = 4, min.segment.length = 0.1, point.padding = unit(0.4,
      "lines"), direction = "both", ylim = c(0, -8), seed = 666, box.padding = 0.1) +
  geom_label_repel(aes(x = pos, y = 0, label = rank, ), fontface = "bold", data = ol.rank,
    size = 3.5, segment.size = 0, force = 4, min.segment.length = 0.1, point.padding = unit(0.4,
      "lines"), direction = "both", ylim = c(0, -8), seed = 666, box.padding = 0.1) +
  theme(axis.title.y = element_text(face='bold',size = 12),
    strip.text.x = element_text(face='bold',size = 12),
    axis.text.x = element_text(face='bold', size = 12),
    axis.text.y = element_text(face='bold', size = 12),
    legend.position = "none") +
  labs(x = "Chromosome", y = "LOD", linetype = "")
dev.off()


qtl_pg <- c(1,5, 10, 12, 23)
ol.rank <- all.rank[which(all.rank$chr %in% qtl_pg), ]
ol.melt <- allmelt[which(allmelt$chr %in% qtl_pg), ]
ol.gens <- nbh.gens[which(nbh.gens$chr %in% qtl_pg), ]


p <- ggplot(ol.melt, aes(x = pos, y = lod, colour = pop))
pdf("/home/jmiller1/public_html/OL_noQTL_chr.qtl.pdf")
p + facet_wrap(~chr, nrow = 1, scales = "free_x", ncol = 24) +
  scale_color_manual(values = popcol) +
  scale_y_continuous(limits = c(-5, 22)) +
  geom_label_repel(aes(x = pos, y = 0.1,label = gene),data = ol.gens,
    color = "black", segment.size = 1, label.padding = unit(0.2, "lines"),
    direction = "y", size = 3, seed = 666, ylim = c(5, 20), segment.alpha = 0.5) +
  geom_line(size = 2, alpha = 0.5) +
  geom_label_repel(aes(x = pos, y = 0.1, label = gene), data = ol.gens,
    color = "black", label.padding = unit(0.2, "lines"), ylim = c(5, 20), segment.size = 0,
    max.iter = 6000, nudge_y = 2, direction = "y", size = 3.5, seed = 666) +
  geom_label_repel(aes(x = pos, y = 0, fill = pop, label = rank), data = ol.rank,
    size = 3.5, segment.size = 1, force = 4, min.segment.length = 0.1, point.padding = unit(0.4,
      "lines"), direction = "both", ylim = c(0, -8), seed = 666, box.padding = 0.1) +
  geom_label_repel(aes(x = pos, y = 0, label = rank, ), fontface = "bold", data = ol.rank,
    size = 3.5, segment.size = 0, force = 4, min.segment.length = 0.1, point.padding = unit(0.4,
      "lines"), direction = "both", ylim = c(0, -8), seed = 666, box.padding = 0.1) +
  theme(axis.title.y = element_text(face='bold',size = 12),
    strip.text.x = element_text(face='bold',size = 12),
    axis.text.x = element_text(face='bold', size = 12),
    axis.text.y = element_text(face='bold', size = 12),
    legend.position = "none") +
  labs(x = "Chromosome", y = "LOD", linetype = "")
dev.off()


### Entropy

NBH <- subset(cross.NBH, ind = cross.NBH$pheno$gt == 1)
png("/home/jmiller1/public_html/nbh_entropy.png", width = 3000)
plotInfo(NBH, chr = c(1:24), main = "NBH", method = "both", include.genofreq = T)
dev.off()

NEW <- subset(cross.NEW, ind = cross.NEW$pheno$gt == 1)
png("/home/jmiller1/public_html/new_entropy.png", width = 3000)
plotInfo(NEW, chr = c(1:24), main = "NEW", method = "both", include.genofreq = T)
dev.off()

ELR <- subset(cross.ELR, ind = cross.ELR$pheno$gt == 1)
png("/home/jmiller1/public_html/elr_entropy.png", width = 3000)
plotInfo(ELR, chr = c(1:24), main = "ELR", method = "both", include.genofreq = T)
dev.off()

BRP <- subset(brp.remap, ind = brp.remap$pheno$gt == 1)
png("/home/jmiller1/public_html/brp_entropy.png", width = 3000)
plotInfo(BRP, chr = c(1:24), main = "BRP", method = "both", include.genofreq = T)
dev.off()

#### Density

pheno <- read.table("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/data/pheno.csv",
  stringsAsFactors = T, header = T,sep=',')
pheno$pop_all <- factor(pheno$pop_all, levels = rev(c('NBH','BRP','NEW','ELR')))
pheno$pheno_all <- factor(pheno$pheno_all, levels = c(NA,0:5))
pheno$gtd <- pheno$pheno_all
pheno[which(pheno$GT_NG_ALT=='NG'),6] <- NA
