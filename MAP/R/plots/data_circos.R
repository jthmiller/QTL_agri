#load(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
library('qtl')
library(circlize)
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'


###########################################################################
#load(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))
load(file.path(mpath,paste0(pop,'_csq_scan.rsave')))
##cross1 <- est.rf(cross)
cross1 <- cross
names(cross1$geno) <- ifelse(names(cross$geno) == "X","5",names(cross$geno))
attr(cross1$geno[["5"]], 'class') <- 'A'

cross1 <- est.rf(cross1)

rf <- pull.rf(cross1)
lod <- pull.rf(cross1, what='lod')
mars <- markernames(cross1)
ahr_genes <- get_AHR(cross1)
###########################################################################

rf <- pull.rf(cross1)
rownames(rf) <- colnames(rf) <- markernames(cross1)

lod <- pull.rf(cross1, what='lod')
rownames(lod) <- colnames(lod) <- markernames(cross1)

lod_hom <- data.matrix(-log10(csq_mod.pval))
rownames(lod_hom) <- colnames(lod_hom) <- markernames(cross1)

lod_inc <- data.matrix(-log10(csq.pval))
rownames(lod_inc) <- colnames(lod_inc) <- markernames(cross1)

s1a <- scanone(cross1, pheno.col=5, model="normal", method="mr")
s1 <- scanone(cross1, pheno.col=5, model="normal", method="mr")
s1 <- s1[mars,'lod']
s1 <- matrix(s1, nrow = length(mars), ncol = length(mars))
rownames(s1) <- colnames(s1) <- mars

###########################################################################
load(file.path(mpath,paste0(pop,'_scan2_normal_mr.rsave')))
lod_phen <- data.matrix(norm.mr.2$lod)
rownames(lod_phen) <- colnames(lod_phen) <- markernames(cross)
###########################################################################

###########################################################################

mars <- grep('NW',markernames(cross1), invert=T,value=T)
mat.names1 <- matrix(mars, nrow = length(mars), ncol = length(mars))
mat.names2 <- t(mat.names1)

rownames(mat.names1) <- colnames(mat.names1) <- mars
rownames(mat.names2) <- colnames(mat.names2) <- mars

###########################################################################

###########################################################################
for (i in chrnames(cross1)){
 mars2 <- markernames(cross1, i)
 mars2 <- mars2[which(mars2 %in% mars)]

 lod[mars2, mars2] <- NA
 lod_phen[mars2, mars2] <- NA
 rf[mars2, mars2] <- NA
 lod_hom[mars2, mars2] <- NA
 lod_inc[mars2, mars2] <- NA
}
###########################################################################

###########################################################################
map <- map2table(pull.map(cross1))
###########################################################################

mat.names1 <- mat.names1 [lower.tri(mat.names1, diag = F)]
mat.names2 <- mat.names2[lower.tri(mat.names2, diag = F)]

rf <- round(rf[cbind(mat.names1,mat.names2)],digits = 3)
lod <- round(lod[cbind(mat.names1,mat.names2)],digits = 3)
lod_hom <- round(lod_hom[cbind(mat.names1,mat.names2)],digits = 3)
lod_inc <- round(lod_inc[cbind(mat.names1,mat.names2)],digits = 3)
lod_phen <- round(lod_phen[cbind(mat.names1,mat.names2)],digits = 3)
s1 <- round(s1[cbind(mat.names1,mat.names2)],digits = 3)

links <- data.frame(cbind(mat.names1,mat.names2,lod,rf,lod_hom,lod_inc,lod_phen,s1),stringsAsFactors=F)

cols <- c('lod','rf','lod_hom','lod_inc','lod_phen','s1')
links_ord <- links[,cols]
links_ord <- apply(links_ord, 2 , rank, na.last = F)

gt <- geno.table(cross1)
pvals <- data.frame(chr = as.character(map[rownames(gt),'chr']),pos = map[rownames(gt),'pos'], pval = -log10(gt$P.value))
rownames(pvals) <- rownames(gt)
################################################################################
################################################################################

################################################################################
################################################################################

facp = pvals[rownames(map),'chr']
xp = pvals[rownames(map),'pos']
yp = pvals[rownames(map),'pval']

facP = s1a[rownames(map),'chr']
xP = s1a[rownames(map),'pos']
yP = s1a[rownames(map),'lod']

################################################################################
################################################################################


################################################################################
################################################################################
cross2 <- convert2cross2(cross1)
map2 <- insert_pseudomarkers(cross2$gmap, step=0.5)
pr2 <- calc_genoprob(cross2, map2, error_prob=0.001, cores=4)
bin2 <- scan1(pr2, pheno=cross2$pheno[,'bin'] , model = "binary", cores = cores)

plot_test(paste0(pop,'_scan2'), width = 2000, height = 250)
plot(bin2, map2, cex=2, lwd=2, col = 'cornflowerblue', cex.axis = 2, ylim = c(0,12))
dev.off()

plot_test(paste0('_chr_scan2'), width = 750, height = 250)
plot(bin2, map2, chr = c(1,2,8,13,18,24), cex=2, lwd=2, col = 'cornflowerblue', cex.axis = 2, ylim = c(0,12))
dev.off()
################################################################################

################################################################################
################################################################################
pc.plot <- function(link, col){
 for (i in 1:length(link[,1])){
 circos.link( map[link[i,'mat.names1'],'chr'], map[link[i,'mat.names1'],'pos'],
   map[link[i,'mat.names2'],'chr'], map[link[i,'mat.names2'],'pos'], col = col,  h = 0.2, border = 1)
 }
}

top_links <- function(links2,n,col){
 tp <- links2[which.max(links2[,col]),col] - n
 return(links[which(links2[,col] > tp),])
}

quant_links <- function(links2, n, col){
 return( links2[ which(col > quantile(col, n, na.rm=T)) ,] )
}
################################################################################
################################################################################

viz_hom <- data.matrix(-log10(csq_mod.pval))
rownames(viz_hom) <- colnames(viz_hom) <- markernames(cross1)

viz_inc <- data.matrix(-log10(csq.pval))
rownames(viz_inc) <- colnames(viz_inc) <- markernames(cross1)

################################################################################

################################################################################
ahr <- ahr_genes[!is.na(ahr_genes$PATH),]

gt <- geno.table(cross1)

gt.count <- rowSums(gt[,c('AA','AB','BB')])

gt.AA <- gt$AA - (gt.count * 0.25)
gt.AB <- gt$AB - (gt.count * 0.5)
gt.BB <- gt$BB - (gt.count * 0.25)
chrs <- as.factor(gt$chr)

ymx <- max(c(gt.AA,gt.AB,gt.BB))
ymn <- min(c(gt.AA,gt.AB,gt.BB))

## vectors of genotype counts
gts <- list(gt.AA,gt.AB,gt.BB)

### Standard deviation
gts.sd <- lapply(gts,sd)
me <- lapply(gts.sd, function(X){  qnorm(.995)*(as.numeric(X)/sqrt(sum(nmar(cross1)))) })
gts.mean <- lapply(gts,mean)

### mean for each vector
gts.H <- mapply(function(X,Y){ X + Y },gts.mean,me)
gts.L <- mapply(function(X,Y){ X - Y },gts.mean,me)

################################################################################
################################################################################
### PREDICTION INTERVALS #######################################################
lm_preds <- lapply(gts,function(X){
 lms <- lm(X~1)
 lms_pred <- predict(lms,  interval="predict")[1,]
 list(lms,lms_pred)
})

outliers <- mapply(function(X,Y) {  which(X < Y[[2]]['lwr'] | X > Y[[2]]['upr'])  },gts,lm_preds)

models <- lapply(gts, function(X){
 lapply(1:24, function(i){
  ind <- which(gt$chr == i)
  x <- map[names(X[ind]),'pos']
  loess(X[ind]~x)
  })
})

### PLOT LINES #################################################################
models.predicted <- lapply(models,function(X){
 lapply(X,predict)
})
################################################################################

################################################################################
save.image(file.path(mpath,paste0(pop,'_circos.rsave')))
################################################################################
