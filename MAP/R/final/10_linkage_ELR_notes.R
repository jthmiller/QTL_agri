#!/bin/R
mpath <- '/home/jmiller1/QTL_agri/data'
pop <- 'ELR'

library('qtl')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")

################################################################################
load(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))
################################################################################

names(cross$geno) <- ifelse(names(cross$geno) == "X","5",names(cross$geno))
attr(cross$geno[["5"]], 'class') <- 'A'

#############################################
ahr_genes <- get_AHR(cross)
gt <- geno.table(cross)
ahr_genes$segdist <- -log10(gt[ahr_genes$close_marker,'P.value'])
ahr_genes_sub <- ahr_genes[!is.na(ahr_genes$PATH),]
#############################################

################ INCOMPATABILITY ##############################################
geno.crosstab(cross,'2:35401205','13:22410721')


#### THESE LOCI IN NBH
c2 :c13 80.59 39.42    27.32    4.68    9.38  17.937 -4.7046






plot_test('sdf')
plot(bin.mr.2)
dev.off()



plot_test('sdf', width = 1000, height = 1000)
 plot(norm.em.2, zmax = c(25,14), col.scheme = "redblue", contours=T, gamma=0.8)
dev.off()

 col.scheme = c("viridis", "redblue","cm","gray","heat","terrain","topo")


plot_test('sdf', width = 1000, height = 1000)
 plot(norm.em.2, zmax = c(25,14), col.scheme = "redblue", contours=T, gamma=0.8)
dev.off()



plot_test('gray', width = 1000, height = 1000)
 plot(norm.em.2, zmax = c(25,14), col.scheme = "gray", contours=T)
dev.off()

plot_test('heat', width = 1000, height = 1000)
 plot(norm.em.2, zmax = c(25,14), col.scheme = "heat", contours=T)
dev.off()

plot_test('terrain', width = 1000, height = 1000)
 plot(norm.em.2, zlim = c(20,10), col.scheme = "terrain", contours=T)
dev.off()

plot_test('topo', width = 1000, height = 1000)
 plot(norm.em.2, zmax = c(25,14), col.scheme = "topo", contours=T)
dev.off()



summary(bin.imp.2, thresholds=c(0, Inf, 5.5, Inf, Inf), what="int")


########################################
### TABLE OF THE HIGHEST LOD SCORES OF LINKAGE FOR EACH CHR
maxdist <- lapply(chrnames(cross), function(i) {
  mars <- markernames(rf, i)
  a <- which(csq.pval[mars,] == min(csq.pval[mars,], na.rm = T), arr.ind=T)
  b <- markernames(rf)[a[,'col']][1]
  a <- rownames(a)[1]
  cbind(a, b, -log10(csq.pval[cbind(a,b)]))
})
maxdist <- do.call(rbind,maxdist)
maxdist <- maxdist[order(as.numeric(maxdist[,3])),]
maxdist <- data.frame(maxdist, stringsAsFactors = F)

##rownames(maxdist)  <- as.character(unique(bin.em.2$map$chr))
######################################################################

##### HEATMAP ######################################################################
csq.pval.hm <- data.matrix(-log10(csq.pval))
plot_test('heatmap_dist_elr',height=3000,width=3000)
heatmap(csq.pval.hm, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column")
dev.off()
######################################################################

######################################################################
### WHITHOUT DISTORTED 17
csq.pval.2 <- csq.pval
mars <- markernames(rf, 17)
csq.pval.2[mars,] <- NA
csq.pval.2[,mars] <- NA
rf.plots.2 <- rf.plots
rf.plots.2$rf[ ,mars] <- NA
rf.plots.2$rf[mars, ] <- NA

### TABLE OF THE HIGHEST LOD SCORES OF LINKAGE FOR EACH CHR
maxdist2 <- lapply(chrnames(cross)[-17], function(i) {
  mars <- markernames(rf, i)
  a <- which(csq.pval.2[mars,] == min(csq.pval.2[mars,], na.rm = T), arr.ind=T)
  b <- markernames(rf)[a[,'col']][1]
  a <- rownames(a)[1]
  cbind(a, b, -log10(csq.pval.2[cbind(a,b)]))
})
maxdist2 <- do.call(rbind,maxdist2)
maxdist2 <- maxdist2[order(as.numeric(maxdist2[,3])),]
maxdist2 <- data.frame(maxdist2, stringsAsFactors = F)

##### HEATMAP ######################################################################
csq.pval.hm <- data.matrix(-log10(csq.pval.2))
plot_test('heatmap_dist_elr',height=3000,width=3000)
heatmap(csq.pval.hm, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column")
dev.off()
######################################################################
######################################################################















ahr_genes_sub <- ahr_genes_sub[which(ahr_genes_sub$PATH == 'AHR'),]
ahr.dist <- csq.pval[ahr_genes_sub$close_marker,]
ahr_col <- as.factor(ahr_genes_sub$chr)

which.max(-log10(csq.pval['22:19528880',]))


plot_test('sfsd', width=5000, height=2000)

plot(1:length(ahr.dist[2,]), -log10(ahr.dist[2,]), pch=19, cex=0.5, ylim=c(0,8))
for(i in 2:11){
#points(1:length(ahr.dist[i,]), -log10(ahr.dist[i,]), pch=19, cex=0.5)
points(1:length(ahr.dist[i,]), -log10(ahr.dist[i,]), col=ahr_col[i], pch=19, cex=0.5)
}
points(1:length(gt[,1]), -log10(gt$P.value), col=as.factor(gt$chr), pch=19, cex=2, ylim=c(0,5))
dev.off()


diag(rf.plots$rf) <- NA

plot_test('CHRxCHR_LOD_scores',height=1000,width=1000)
plotRF(rf.plots,zmax=4,col.scheme="redblue")
dev.off()


a <- which(csq.pval[markernames(cross,22),] == min(csq.pval[markernames(cross,22),],na.rm=T ), arr.ind =T)
markernames(cross)[43:49]
markernames(cross)[196:197]

which.min(csq.pval['20:19045496',])

###############

csq.pval[ahr_genes_sub$close_marker,ahr_genes_sub$close_marker]


crs <- pull.markers(rf.plots, ahr_genes$close_marker)

plot_test('qtlxqtl_SegLOD_scores')
plotRF(crs,zmax=7,col.scheme="redblue")
dev.off()


which(csq.pval[markernames(cross,20),] == min(csq.pval[markernames(cross,20),],na.rm=T ), arr.ind =T)
