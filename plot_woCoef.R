#load(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

load(file.path(mpath,paste0(pop,'_scan2_bin_em_noCof.rsave')))
## with coefload(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
library(circlize)

###############
rf <- subset(cross, chr = c(1:4,6:24))
rf <- est.rf(rf, maxit=100000, tol=1e-6)
mars <- find.marker(rf, bin.em.2$map$chr, bin.em.2$map$pos)
###############

nbh_gens <- cnv.ahrs(rf, AHRdf = AHR.bed, EXP = F)
elr_gens <- cnv.ahrs(rf, AHRdf = AHR.bed, EXP = F)

###########################################################################
rf.df <- pull.rf(rf)
rf.df <- rf.df[mars,mars]
diag(rf.df) = 0
rf.df[lower.tri(rf.df)] = 0
###########################################################################

###############
mar.names <- matrix(mars, nrow = dim(rf.df)[1], ncol = dim(rf.df)[2])
mat.names <- matrix(mars, nrow = dim(rf.df)[1], ncol = dim(rf.df)[2])
mat.names <- gsub(":.*","",mat.names)
###############

###############
lod.df <- pull.rf(rf, what='lod')
mat = lod.df
mat <- mat[mars,mars]
diag(mat) = 0
mat[lower.tri(mat)] = 0
n = nrow(mat)
rn = rownames(mat)

diag(mat.names) = 0
mat.names[lower.tri(mat.names)] = 0
n = nrow(mat.names)

rownames(mat.names) <- rownames(mat)
colnames(mat.names) <- rownames(mat)
rn = rownames(mat.names)
###############

###############

###############
s1 <- scanone(rf,pheno.col=4, model="binary", method="em")
s1l <- matrix(s1$lod, nrow = dim(s1), ncol = dim(s1))
s1l[lower.tri(s1l, diag = T)] <- NA
###############

###############
lod_phen <- bin.em.2$lod
lod_phen[lower.tri(lod_phen, diag = T)] <- NA
diag(lod_phen) = 0
lod_phen[lower.tri(lod_phen)] = 0

rownames(lod_phen) <- rownames(mat)
colnames(lod_phen) <- rownames(mat)
###############

###############
map <- map2table(pull.map(rf))

for (i in unique(bin.em.2$map$chr)){
 mars <- which(bin.em.2$map$chr == i)
 mat[mars,mars] <- 0
}
###########################################################################

###########################################################################
###############
### test seg dist

probs <- c(0.0625,0.125,0.25)
gts <- c('AA','AB','BB')

homs <- c('AA','BB')
hets <- 'AB'

#homs <- c('1','3')
#hets <- '2'

tr.table <- matrix(NA, ncol=3, nrow=3)
rownames(tr.table) <- colnames(tr.table) <- gts

tr.table[homs,homs] <- 0.0625
tr.table[hets,homs] <- 0.125
tr.table[homs,hets] <- 0.125
tr.table[hets,hets] <- 0.25

gtf <- c('AA','AB','BB')
gt_gt <- cbind(rep(gtf,3),c(rep('AA',3),rep('AB',3),rep('BB',3)))
gt_names <- paste0(gt_gt[,1],gt_gt[,2])
gt_probs <- setNames(tr.table[gt_gt], gt_names)

rf.gts <- pull.geno(rf)

csq <- function(mara, marb) {
 test <- factor(paste0(factor(mara, labels = gtf), factor(marb, labels = gtf)), levels = gt_names)
 chisq.test(table(test), p = gt_probs)$p.value
}

#### long ################################
csq.pval <- apply(rf.gts, 2, function(X){
 apply(rf.gts, 2, csq, marb = X)
})
colnames(csq.pval) <- rownames(csq.pval) <- colnames(rf.gts)

csq.bk <- csq.pval
########################################

## Set within chromosomes to zero #####
for (i in unique(bin.em.2$map$chr)){
 mars <- markernames(rf, i)
 csq.pval[mars,mars] <- NA
}
########################################

maxdist <- lapply(as.character(unique(bin.em.2$map$chr)), function(i) {
  mars <- markernames(rf, i)
  a <- which(csq.pval[mars,] == min(csq.pval[mars,], na.rm = T), arr.ind=T)
  b <- markernames(rf)[a[,'col']]
  a <- rownames(a)
  cbind(a, b, -log10(csq.pval[cbind(a,b)]))[1,]
})
maxdist <- do.call(rbind,maxdist)
maxdist <- maxdist[order(as.numeric(maxdist[,3])),]
maxdist <- data.frame(maxdist, stringsAsFactors = F)
rownames(maxdist)  <- as.character(unique(bin.em.2$map$chr))
###########################################################################
###########################################################################

###########################################################################
###########################################################################
mat.dwn <- which(mat > 0 & mat < 100, arr.ind=T)
mar_b <- colnames(mat)[mat.dwn[,'col']]
mar_a <- rownames(mat.dwn)
lod_ab <- unlist(lod.df)
lod_ab <- mapply(function(x,y){ lod_ab[x,y, drop = T] },mar_a,mar_b)
lod_p <- mapply(function(x,y){ lod_phen[x,y, drop = T] },mar_a,mar_b)
rfs <- mapply(function(x,y){ rf.df[x,y, drop = T] },mar_a,mar_b)
lod_seg.dist <- mapply(function(x,y){ -log10(csq.pval[x,y, drop = T]) },mar_a,mar_b)

chr_a <- map[mar_a,c('chr','pos')]
chr_b <- map[mar_b,c('chr','pos')]
links <- cbind(chr_a,chr_b,lod_ab,lod_p,rfs,lod_seg.dist)
###########################################################################

crstb <- function(X,Y) {
 print(geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = X, mname2 = Y))
 print(geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = X, mname2 = Y))
 print(geno.crosstab(rf, mname1 = X, mname2 = Y))
}
crstb(X = '1:6612852', Y = "19:42466531")

################################################################################

#save.image(file.path(mpath,paste0(pop,'woCoef_circos.rsave')))
load(file.path(mpath,paste0(pop,'_circos_wo_Coef.rsave')))
################################################################################

pc <- function(links){

circos.par("track.height" = 0.1)
circos.initialize(factors = map$chr, x = map$pos)

circos.track(factor = map$chr, y = map$pos,
    panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"),
            CELL_META$sector.index)
        circos.axis(labels.cex = 0.6)
})
for(i in 1:length(links[,1])){
 circos.link(links[i,1], links[i,2],links[i,3], links[i,4], h = 0.4)
 }
}

################################################################################
### PLOTs ######################################################################

h.segdist <- quantile(links$lod_seg.dist, 0.99995,na.rm=T)
l.segdist <- quantile(links$lod_seg.dist, 0.001,na.rm=T)


hp <- quantile(links$lod_p, 0.99,na.rm=T)
hab <- quantile(links$lod_ab, 0.99,na.rm=T)

hrf <- quantile(links$rfs, 0.98,na.rm=T)
lrf <- quantile(links$rfs, 0.01,na.rm=T)


high_2seg <- links[which(links$lod_seg.dist > h.segdist),]
plot_test(paste0(pop,'circ_hi_2segdist'))
pc(high_2seg)
dev.off()

hdis <- quantile(links$lod_seg.dist, 0.999,na.rm=T)
hab <- quantile(links$lod_ab, 0.9999,na.rm=T)

low_p <- links[which(links$lod_seg.dist > hdis & links$lod_ab > hab),]
plot_test(paste0(pop,'high_dis_hi_ab'))
pc(low_p)
dev.off()

hdis <- quantile(links$lod_seg.dist, 0.999,na.rm=T)
hab <- quantile(links$lod_p, 0.99, na.rm=T)

low_p <- links[which(links$lod_seg.dist > hdis & links$lod_p > hab),]
plot_test(paste0(pop,'high_dis_hilo_rf'))
pc(low_p)
dev.off()


lp <- quantile(links$lod_p, 0.85,na.rm=T)
hab <- quantile(links$lod_ab, 0.9999,na.rm=T)
low_p <- links[which(links$lod_p < lp & links$lod_ab > hab),]
plot_test(paste0(pop,'circ_low_p_hi_ab'))
pc(low_p)
dev.off()

hp <- quantile(links$lod_p, 0.98,na.rm=T)
hab <- quantile(links$lod_ab, 0.98,na.rm=T)
hi_p <- links[which(links$lod_p > hp & links$lod_ab > hab),]
plot_test(paste0(pop,'_circ_hi_p_hi_ab_98'))
pc(hi_p)
dev.off()

hrf <- quantile(links$rfs, 0.9999,na.rm=T)
lrf <- quantile(links$rfs, 0.0001,na.rm=T)
hp <- quantile(links$lod_p, 0.85,na.rm=T)
hilo_rf <- links[which(links$lod_p > hp & links$rfs < lrf | links$lod_p > hp & links$rfs > hrf),]
plot_test(paste0(pop,'circ_hilo_rf_med_dp'))
pc(hilo_rf)
dev.off()

mhrf <- quantile(links$rfs, 0.99,na.rm=T)
mlrf <- quantile(links$rfs, 0.01,na.rm=T)
hp <- quantile(links$lod_p, 0.99,na.rm=T)
hi_p <- links[which(links$lod_p > hp & links$rfs > mhrf | links$lod_p > hp & links$rfs < mlrf ),]
plot_test(paste0(pop,'circ_hi_rf_hi_p'))
pc(hi_p)
dev.off()

################################################################################
hrf <- quantile(links$rfs, 0.98,na.rm=T)
lrf <- quantile(links$rfs, 0.02,na.rm=T)

hp <- quantile(links$lod_p, 0.999,na.rm=T)
hab <- quantile(links$lod_ab, 0.995,na.rm=T)

rf_p <- links[which(links$lod_p > hp & links$rfs < lrf | links$lod_p > hp & links$rfs > hrf ),]
rf_p <- rf_p[order(rf_p$rfs),]

ab_p <- links[which(links$lod_p > hp & links$lod_ab > hab ),]
ab_p <- ab_p[order(ab_p$lod_ab),]




################################################################################
plot_test(paste0(pop,'cors_lod'), width= 1500, height=1000)
 par(mfrow=c(1,2))
 plot(links$lod_p, links$lod_ab, pch=19, cex=0.5, col = NA)
 text(links$lod_p, links$lod_ab, links[,1])
 plot(links$lod_p, links$lod_ab, pch=19, cex=0.5, col = NA)
 text(links$lod_p, links$lod_ab, links[,3])
dev.off()

plot_test(paste0(pop,'cors_rf'), width= 1500, height=1000)
 par(mfrow=c(1,2))
 plot(links$lod_p, links$rfs, pch=19, cex=0.5, col = NA)
 text(links$lod_p, links$rfs, links[,1])
 plot(links$lod_p, links$rfs, pch=19, cex=0.5, col = NA)
 text(links$lod_p, links$rfs, links[,3])
dev.off()

#col_fun = colorRamp2(c(0, 4), c("white", "red"), transparency = 0.5)
