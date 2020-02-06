#!/bin/R
pop <- 'NBH'
library('qtl')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

################################################################################
## perms.1
## perms.2
## pens
load(file.path(mpath,paste0(pop,'_all_perms_bin_hk.rsave')))
################################################################################
## bin.em.2
load(file.path(mpath,paste0(pop,'_scan2_bin_hk.rsave')))
load(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))

load(file.path(mpath,paste0(pop,'_scan2_bin_em_noCof.rsave')))
################################################################################
#sone.o <- scanone(cross,pheno.col=4, model="binary", method="em")
#sone.a <- scanone(cross,pheno.col=4, model="binary", method="em", addcovar=g[,1])
#sone.i <- scanone(cross,pheno.col=4, model="binary", method="em", addcovar=g[,1],intcovar=g[,1])
#sone.io <- scanone(cross,pheno.col=4, model="binary", method="em", addcovar=g[,1],intcovar=g[,1])
#cbind(summary(sone.o),summary(sone.a)$lod,summary(sone.i)$lod,summary(sone.io)$lod)

rf <- subset(cross, chr = c(1:4,6:24))
rf <- est.rf(rf, maxit=100000, tol=1e-6)

s1 <- scanone(rf,pheno.col=4, model="binary", method="em")
s1l <- matrix(s1$lod, nrow = 1991, ncol = 1991)
s1l[lower.tri(s1l, diag = T)] <- NA

mars <- find.marker(rf, bin.em.2$map$chr, bin.em.2$map$pos)
rf.df <- pull.rf(rf)
rf.df <- rf.df[mars,mars]
rf.df[lower.tri(rf.df, diag = T)] <- NA

lod.df <- pull.rf(rf, what='lod')
lod.df <- lod.df [mars,mars]
lod.df[lower.tri(lod.df, diag = T)] <- NA

lod_phen <- bin.em.2$lod
lod_phen[lower.tri(lod_phen, diag = T)] <- NA

#mars <- unlist(lapply(attr(bin.em.2,"fullmap"),names))
#rownames(lod_phen) <- rownames(s1l) <- mars
#colnames(lod_phen) <- colnames(s1l) <- mars

for (i in unique(bin.em.2$map$chr)){
 ind <- which(bin.em.2$map$chr == i)
 s1l[ind,ind] <- NA
 lod_phen[ind,ind] <- NA
 rf.df[ind,ind] <- NA
 lod.df[ind,ind] <- NA
}

mar.names <- matrix(mars, nrow = dim(rf.df)[1], ncol = dim(rf.df)[2])
mat.names <- matrix(mars, nrow = dim(rf.df)[1], ncol = dim(rf.df)[2])
mat.names <- gsub(":.*","",mat.names)

################################################################################
AHR.bed <- read.table(file.path(mpath,"lift_AHR_genes.bed"), stringsAsFactors = F, header = F)
colnames(AHR.bed) <- c("chrom", "str", "stp", "gene")
AHR.bed$chrom <- as.numeric(gsub("chr", "", AHR.bed$chrom))
AHR.bed$str <- as.numeric(AHR.bed$str)
AHR.bed$stp <- as.numeric(AHR.bed$stp)
AHR.notmap <- AHR.bed[is.na(AHR.bed$chrom), ]
AHR.bed <- AHR.bed[!is.na(AHR.bed$chrom), ]
AHR.bed$gene <- gsub(":158640", "", AHR.bed$gene)
AHR.bed <- AHR.bed[!AHR.bed$chr == 5,]
nbh_gens <- cnv.ahrs(rf, AHRdf = AHR.bed, EXP = F)
################################################################################

################################################################################
plot_test('rf', width=1250,height=1250)
 plot(lod_phen[!is.na(lod.df)],lod.df[!is.na(lod.df)], col = 'grey', pch=19, xlim=c(0,26))
dev.off()

plot_test('phen_rf', width=1250,height=1250)
 plot(lod_phen[!is.na(lod_phen)],rf.df[!is.na(lod_phen)], col = 'grey', pch=19, xlim=c(0,26))
dev.off()

plot_test('sc1_lod_link', width=1250,height=1250)
 plot(s1l[!is.na(rf.df)],rf.df[!is.na(rf.df)], col = 'grey', pch=19, xlim=c(0,26))
dev.off()
################################################################################

hrf <- quantile(rf.df, 0.999,na.rm=T)
lrf <- quantile(rf.df, 0.001,na.rm=T)

rf.cut <- which(rf.df < lrf | rf.df > hrf)

lod.cut <- which(lod_phen > 6)
ind <- intersect(rf.cut,lod.cut)
hi_lod <- sort(table(t(mat.names)[ind]))

lod.cut <- which(lod_phen > 3 & lod_phen < 6 )
ind <- intersect(rf.cut,lod.cut)
med_lod <- sort(table(t(mat.names)[ind]))

lod.cut <- which(lod_phen < 3)
ind <- intersect(rf.cut,lod.cut)
low_lod <- sort(table(t(mat.names)[ind]))

hi_lod
med_lod
low_lod

high <- c(10,2,1,18)
med <- c(3,8,24,19)
low <- c(13,11,12,17,15,23)




3,16,

################################################################################
################################################################################

chr1 <- gsub(":.*","",colnames(rf.df)) %in% c(1)
chr18 <- gsub(":.*","",colnames(rf.df)) %in% c(18)
chr2 <- gsub(":.*","",colnames(rf.df)) %in% c(2)

col <- matrix('grey', nrow = 1991, ncol = 1991)

col[,chr1] <- 'black'
col[,chr1] <- 'black'

col[ind1,] <- 'black'
col[,ind2] <- 'black'
col[ind2,] <- 'black'
col[ind1,ind2] <- 'red'
col[ind2,ind1] <- 'red'

plot_test('rf', width=1250,height=1250)
 plot(lod_phen[col == 'grey'],rf.df[col == 'grey'], col = 'grey', pch=19, xlim=c(0,26))
 #plot(lod_phen[col == 'black'],rf.df[col == 'black'], mat.names[col == 'black'], col = 'black')
 text(lod_phen[col == 'black'] ,rf.df[col == 'black'], mat.names[col == 'black'], col = col[col == 'black'])
 text(lod_phen[col == 'red'],rf.df[col == 'red'], mat.names[col == 'red'], col = col[col == 'red'])
dev.off()


1, 18
2, 8


y <- sort(lod.df)[!is.na(lod.df)]
x <- 1:length(lod.df[!is.na(lod.df)])
dat <- data.frame(cbind(x,y))

model <- lm(formula = y ~ x,data=dat)

newdat <- data.frame(y)
newdat <- predict(model, data.frame(x))


plot_test('rf', width=250,height=250)
plot(newdat,y, xlim=c(0,3),ylim=c(0,3))
abline(lm(newdat ~ y))
dev.off()











colnames(lod_phen) <- markernames(rf)
rownames(lod_phen) <- markernames(rf)




plot_test('rf', width=1250, height=1250)
plotRF(rf,zmax=8, col.scheme="redblue")
dev.off()


############################################

############################################
no_qtl_mr <- scanone(cross, pheno.col=4, method="mr", model="binary")
qtl <- summary(no_qtl_mr, 5)

for (i in chrnames(rf)){
 ind <- which(markernames(rf) %in% markernames(rf,i))
 rf$rf[ind,ind] <- NA
}



h <- quantile(rf.df, 0.999,na.rm=T)
l <- quantile(rf.df, 0.001,na.rm=T)
h9 <- quantile(rf.df, 0.95,na.rm=T)
l9 <- quantile(rf.df, 0.05,na.rm=T)

plot_test('nbh_1_2_8_18', width=1250,height=1250)
par(mfrow = c(9,1))

plot(pull.rf(rf), rownames(summary(no_qtl_mr))[17], ylim=c(0.3,0.75))
abline(h=h, col='red')
abline(h=h9, col='grey')
abline(h=0.5, col='black')
abline(h=l9, col='grey')
abline(h=l, col='red')

plot(pull.rf(rf), rownames(summary(no_qtl_mr))[12], ylim=c(0.3,0.75))
abline(h=h, col='red')
abline(h=h9, col='grey')
abline(h=0.5, col='black')
abline(h=l9, col='grey')
abline(h=l, col='red')


plot(pull.rf(rf), rownames(summary(no_qtl_mr))[2], ylim=c(0.3,0.75))
abline(h=h, col='red')
abline(h=h9, col='grey')
abline(h=0.5, col='black')
abline(h=l9, col='grey')
abline(h=l, col='red')


plot(pull.rf(rf), rownames(summary(no_qtl_mr))[7], ylim=c(0.3,0.75))
abline(h=h, col='red')
abline(h=h9, col='grey')
abline(h=0.5, col='black')
abline(h=l9, col='grey')
abline(h=l, col='red')


plot(pull.rf(rf), rownames(summary(no_qtl_mr))[23], ylim=c(0.3,0.75))
abline(h=h, col='red')
abline(h=h9, col='grey')
abline(h=0.5, col='black')
abline(h=l9, col='grey')
abline(h=l, col='red')


plot(pull.rf(rf), find.marker(rf,1,0), ylim=c(0.3,0.75))
abline(h=h, col='red')
abline(h=h9, col='grey')
abline(h=0.5, col='black')
abline(h=l9, col='grey')
abline(h=l, col='red')

plot(pull.rf(rf), find.marker(rf,10,45), ylim=c(0.3,0.75))
abline(h=h, col='red')
abline(h=h9, col='grey')
abline(h=0.5, col='black')
abline(h=l9, col='grey')
abline(h=l, col='red')

plot(pull.rf(rf), find.marker(rf,3,17), ylim=c(0.3,0.75))
abline(h=h, col='red')
abline(h=h9, col='grey')
abline(h=0.5, col='black')
abline(h=l9, col='grey')
abline(h=l, col='red')

plot(pull.rf(rf), find.marker(rf,16,15.48), ylim=c(0.3,0.75))
abline(h=h, col='red')
abline(h=h9, col='grey')
abline(h=0.5, col='black')
abline(h=l9, col='grey')
abline(h=l, col='red')


dev.off()
########################################################################################
summary(no_qtl_mr)
########################################################################################

h <- quantile(pull.rf(rf, what='lod'), 0.99,na.rm=T)
l <- quantile(pull.rf(rf, what='lod'), 0.01,na.rm=T)


plot_test('rf', width=1250)
par(mfrow = c(3,1))
plot(pull.rf(rf,what='lod'), rownames(summary(no_qtl_mr))[12], ylim=c(-0.50,2))
abline(h=h)
abline(h=l)
plot(pull.rf(rf,what='lod'), rownames(summary(no_qtl_mr))[17], ylim=c(-0.50,2))
abline(h=h)
abline(h=l)
plot(pull.rf(rf,what='lod'), rownames(summary(no_qtl_mr))[10], ylim=c(-0.50,2))
abline(h=h)
abline(h=l)
dev.off()



summary(bin.em.2, thresholds=c(0, Inf, 5, Inf, Inf), what="int")
summary(bin.em.2, thresholds=c(16, 0, 0, 0, 0), what="full")


rfm <- matrix(pull.rf(rf), nrow = 1618, ncol = 1618)
#ind <- gsub(":.*","",markernames(rf)) %in% c(1,2,8,18,13,24)

ind1 <- gsub(":.*","",markernames(rf)) %in% c(2)
ind2 <- gsub(":.*","",markernames(rf)) %in% c(18)

col <- matrix('black', nrow = 1618, ncol = 1618)
col[,ind] <- 'grey'
col[ind,] <- 'grey'
col[ind1,ind2] <- 'red'
col[ind2,ind1] <- 'red'

plot_test('rf', width=1250,height=1250)
 plot(lod_phen[col == 'black'],rfm[col == 'black'], pch=19, xlim= c(0,25))
 points(lod_phen[col == 'grey'] ,rfm[col == 'grey'], pch=19, col = col)
 points(lod_phen[col == 'red'],rfm[col == 'red'], pch=19, col = col)
dev.off()




plot_test('rf', width=1250,height=1250)
 plot(lod_phen,rfm, col = NA)
 text(lod_phen[col == 'black'],rfm[col == 'black'], mat.names[col == 'black'], pch=19, col = 'black')
 text(lod_phen[col == 'grey'] ,rfm[col == 'grey'], mat.names[col == 'grey'], pch=19, col = col)
 text(lod_phen[col == 'red'],rfm[col == 'red'], mat.names[col == 'red'], pch=19, col = col)
dev.off()

## correlate recombination frequency and lod w phenotype


mat.names <- matrix(mars, nrow = 1991, ncol = 1991)
mat.names <- gsub(":.*","",mat.names)


ind1 <- gsub(":.*","",colnames(rf.df)) %in% c(2)
ind2 <- gsub(":.*","",colnames(rf.df)) %in% c(18)

col <- matrix('black', nrow = 1991, ncol = 1991)
col[,ind1] <- 'grey'
col[,ind2] <- 'grey'
col[ind1,ind2] <- 'red'
col[ind2,ind1] <- 'red'

plot_test('rf', width=1250,height=1250)
 plot(lod_phen,rf.df, col = NA)
 text(lod_phen[col == 'black'],rf.df[col == 'black'], mat.names[col == 'black'], pch=19, col = 'black')
 text(lod_phen[col == 'grey'] ,rf.df[col == 'grey'], mat.names[col == 'grey'], pch=19, col = col)
 text(lod_phen[col == 'red'],rf.df[col == 'red'], mat.names[col == 'red'], pch=19, col = col)
dev.off()



summary(no_qtl_mr)

plot_test('effect_chr1')
effectplot(rf, pheno.col=4, mname1 = find.marker(rf,1,10), var.flag="pooled")
dev.off()


plot_test('effect_chr1_chr18')
effectplot(rf, pheno.col=1, mname1 = find.marker(rf,1,10), mname2 = find.marker(rf, 18, 46.54), var.flag="group")
dev.off()

effectplot(rf, pheno.col=1, mname1, mark1, geno1, mname2, mark2, geno2, main,
ylim, xlab, ylab, col, add.legend=TRUE, legend.lab, draw=TRUE, var.flag=c("pooled","group"))


###########
qtl <- summary(no_qtl_mr, 5)
Q2 <- makeqtl(rf, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")
Q2 <- refineqtl(rf, pheno.col = 4, qtl=Q2, method = "imp", model='binary',incl.markers=T)

plot_test('lod_prof',width=1000)
plotLodProfile(Q2, incl.markers=TRUE, gap=25, lwd=2, lty=1, col="black", qtl.labels=TRUE,showallchr=T, labelsep=5)
dev.off()

plot_test('pxg1',width=666)
plotPXG(rf, find.marker(rf,1,10), pheno.col=4, jitter=2, infer=TRUE)
dev.off()

plot_test('pxg18',width=666)
plotPXG(rf, find.marker(rf, 18, 46.54), pheno.col=4, jitter=2, infer=TRUE)
dev.off()


geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,1,10), mname2 = find.marker(rf, 18, 46.54))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,1,10), mname2 = find.marker(rf, 18, 46.54))
geno.crosstab(rf,mname1 = find.marker(rf,1,10), mname2 = find.marker(rf, 18, 46.54))

geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,1,10), mname2 = find.marker(rf, 2, 102.7))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,1,10), mname2 = find.marker(rf, 2, 102.7))
geno.crosstab(rf,mname1 = find.marker(rf,1,10), mname2 = find.marker(rf, 2, 102.7))


geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,1,8), mname2 = find.marker(rf, 18, 48))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,1,8), mname2 = find.marker(rf, 18, 48))
geno.crosstab(rf,mname1 = find.marker(rf,1,8), mname2 = find.marker(rf, 18, 48))

geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,2,102), mname2 = find.marker(rf, 18, 48))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,2,102), mname2 = find.marker(rf, 18, 48))
geno.crosstab(rf,mname1 = find.marker(rf,2,102), mname2 = find.marker(rf, 18, 48))

geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,10,45), mname2 = find.marker(rf, 18, 57))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,10,45), mname2 = find.marker(rf, 18, 57))
geno.crosstab(rf,mname1 = find.marker(rf,10,45), mname2 = find.marker(rf, 18, 57))

geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,10,45), mname2 = find.marker(rf, 2, 102))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,10,45), mname2 = find.marker(rf, 2, 102))
geno.crosstab(rf,mname1 = find.marker(rf,10,45), mname2 = find.marker(rf, 2, 102))
