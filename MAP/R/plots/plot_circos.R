#load(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
library('qtl')
library(circlize)
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'

library(circlize)

##load(file.path(mpath,paste0(pop,'_scan2_bin_em_noCof.rsave')))
#load(file.path(mpath,paste0(pop,'_csq_scan.rsave')))
load(file.path(mpath,paste0(pop,'_scan2_normal_em.rsave')))

load(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))

###############
ahr_genes <- get_AHR(cross)

nbh_gens <- cnv.ahrs(rf, AHRdf = AHR.bed, EXP = F)
elr_gens <- cnv.ahrs(rf, AHRdf = AHR.bed, EXP = F)

###############
rf <- drop.markers(cross, grep('NW',markernames(cross), value=T))
mars <- markernames(rf)
###############

###############
rf.df <- pull.rf(rf)
rf.df <- rf.df[mars,mars]
diag(rf.df) = 0
rf.df[lower.tri(rf.df)] = 0
###################

###############
mar.names <- matrix(mars, nrow = dim(rf.df)[1], ncol = dim(rf.df)[2])
mat.names <- matrix(mars, nrow = dim(rf.df)[1], ncol = dim(rf.df)[2])
mat.names <- gsub(":.*","",mat.names)
###############

###############
lod.df <- pull.rf(rf, what='lod')
lod_ab <- matrix(lod.df, nrow = sum(nmar(rf)), ncol = sum(nmar(rf)))
rownames(lod_ab) <- markernames(rf)
colnames(lod_ab) <- markernames(rf)

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
s1 <- scanone(rf,pheno.col=4, model="binary", method="mr")
s1l <- matrix(s1$lod, nrow = sum(nmar(rf)), ncol = sum(nmar(rf)))
s1l[lower.tri(s1l, diag = T)] <- NA
###############

###############
lod_phen <- norm.em.2$lod
lod_phen[lower.tri(lod_phen, diag = T)] <- NA
diag(lod_phen) = 0
lod_phen[lower.tri(lod_phen)] = 0

lod_phen <- matrix(norm.em.2$lod, nrow = 7457, ncol = 7457)
rownames(lod_phen) <- markernames(cross)
colnames(lod_phen) <- markernames(cross)

###############

###########################################################################

homozy <- data.matrix(-log10(csq_mod.pval[mars,mars]))
#homozy[lower.tri(homozy, diag = T)] <- NA
#diag(homozy) = 0
#homozy[lower.tri(homozy)] = 0

###########################################################################

incompat <- data.matrix(-log10(csq.pval[mars,mars]))
#incompat[lower.tri(incompat, diag = T)] <- NA
#diag(incompat) = 0
#incompat[lower.tri(incompat)] = 0

###############
map <- map2table(pull.map(rf))
#
#for (i in unique(bin.em.2$map$chr)){
# mars <- which(bin.em.2$map$chr == i)
# mat[mars,mars] <- 0
#}
###########################################################################

#mat.dwn <- which(mat > 0 & mat < 100, arr.ind=T)
#mar_b <- colnames(mat)[mat.dwn[,'col']]
#mar_a <- rownames(mat.dwn)
##lod_ab <- unlist(lod.df)
#lod_ab <- mapply(function(x,y){ lod_ab[x,y, drop = T] }, mar_a, mar_b)
#lod_p <- mapply(function(x,y){ lod_phen[x,y, drop = T] }, mar_a, mar_b)
#rfs <- mapply(function(x,y){ rf.df[x,y, drop = T] }, mar_a, mar_b)
#lod_homz <- mapply(function(x,y){ homozy[x,y, drop = T] }, mar_a, mar_b)
#lod_inco <- mapply(function(x,y){ incompat[x,y, drop = T] }, mar_a, mar_b)




lod_ab <- lod_ab[cbind(mar_a, mar_b), drop = T]
#lod_p <- lod_phen[cbind( mar_a, mar_b), drop = T]
rfs <- rf.df[cbind( mar_a, mar_b), drop = T]
lod_homz <- homozy[cbind( mar_a, mar_b), drop = T]
lod_inco <- incompat[cbind(mar_a, mar_b), drop = T]


chr_a <- map[mar_a, c('chr','pos')]
chr_b <- map[mar_b, c('chr','pos')]

##links <- data.frame(cbind(chr_a,chr_b,lod_ab,lod_p,rfs),stringsAsFactors=F)

links <- data.frame(cbind(chr_a,chr_b,lod_ab,rfs,lod_homz,lod_inco),stringsAsFactors=F)

###########################################################################
##save.image(file.path(mpath,paste0(pop,'_circos_wo_Coef.rsave')))
save.image(file.path(mpath,paste0(pop,'_circos.rsave')))


###########################################################################
###########################################################################
###########################################################################

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


make_lodrf_tables <- function(X,Y,Z){
 df <- Y[which(Y[,1] == X | Y[,3] == X),]
 ch <- unique(c(df[,1],df[,3]))

 xmax <- lapply(ch, function(chx){
   z <- which.max(df[which(df[,1] == chx | df[,3] == chx),][,Z])
   df[which(df[,1] == chx | df[,3] == chx),][z,]
  })
 xmax <- data.frame(do.call(rbind,xmax))
 ##ind <- order(xmin[,Z])
 ##print(head(xmin))
 xmax <- xmax[order(xmax[,Z],decreasing = T),]


 xmin <- lapply(ch, function(chx){
   z <- which.min(df[which(df[,1] == chx | df[,3] == chx),][,Z])
   df[which(df[,1] == chx | df[,3] == chx),][z,]
  })
 xmin <- data.frame(do.call(rbind,xmin))
 ##ind <- order(xmin[,Z])
 ##print(order(xmin[,Z]))
 xmin <- xmin[order(xmin[,Z]),]

 list(xmin=xmin,xmax=xmax)

}
chroms <- unique(c(as.character(links[,1]),as.character(links[,3])))
ab_tables <- lapply(chroms, make_lodrf_tables, Y = links, Z = 'lod_ab')
p_tables <- lapply(chroms, make_lodrf_tables, Y = links, Z = 'lod_p')
rf_tables <- lapply(chroms, make_lodrf_tables, Y = links, Z = 'rfs')
hom_tables <- lapply(chroms, make_lodrf_tables, Y = links, Z = 'lod_homz')
dist_tables <- lapply(chroms, make_lodrf_tables, Y = links, Z = 'lod_inco')


names(ab_tables) <- names(p_tables) <- names(rf_tables) <- chroms


################################################################################
pop <- paste0(pop,'_wCoef')



hz <- quantile(links$lod_homz, 0.9999, na.rm=T)
hi <- quantile(links$lod_inco, 0.9999, na.rm=T)

low_p <- links[which(links$lod_homz > hz ),]
low_p <- links[which(links$lod_inco > hi),]
low_p <- links[which(links$lod_homz > hz | links$lod_inco > hi),]

dim(low_p)
plot_test(paste0(pop,'circ_incompat'))
pc(low_p)
dev.off()


geno.crosstab(cross, '1:291287','19:38646154')


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

plot_test(paste0(pop,'cors_lod'), width= 3000, height=1000)
 par(mfrow=c(1,2))
 plot(links$lod_p, links$lod_ab, pch=19, cex=0.5, col = NA)
 text(links$lod_p, links$lod_ab, links[,1])
 plot(links$lod_p, links$lod_ab, pch=19, cex=0.5, col = NA)
 text(links$lod_p, links$lod_ab, links[,3])
dev.off()

plot_test(paste0(pop,'cors_rf'), width= 3000, height=1000)
 par(mfrow=c(1,2))
 plot(links$lod_p, links$rfs, pch=19, cex=0.5, col = NA)
 text(links$lod_p, links$rfs, links[,1])
 plot(links$lod_p, links$rfs, pch=19, cex=0.5, col = NA)
 text(links$lod_p, links$rfs, links[,3])
dev.off()

#col_fun = colorRamp2(c(0, 4), c("white", "red"), transparency = 0.5)

################################################################################

mhrf <- quantile(links$rfs, 0.99,na.rm=T)
mlrf <- quantile(links$rfs, 0.01,na.rm=T)
hp <- quantile(links$lod_p, 0.99,na.rm=T)

inc_13_18 <- links[which(links$lod_p > 6 & links$rfs > 0.6),]

inc_13_18[order(inc_13_18[,'rfs']),]
