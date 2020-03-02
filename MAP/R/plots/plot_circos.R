#load(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
library('qtl')
library(circlize)
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'


###########################################################################
#load(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))
load(file.path(mpath,paste0(pop,'_csq_scan.rsave')))
cross1 <- est.rf(cross)

rf <- pull.rf(cross1)
lod <- pull.rf(cross1, what='lod')

mf1 <- file.path(mpath,paste0(pop,'_rf.tsv'))
write.table(rf,mf1)

mf2 <- file.path(mpath,paste0(pop,'_lod.tsv'))
write.table(lod,mf2)

test <- read.table(mf1)

###########################################################################

###############
ahr_genes <- get_AHR(cross1)
###############

###############
map <- map2table(pull.map(cross1))
mars <- grep('NW',markernames(cross1), invert=T,value=T)
###############

###############
mat <- rf
mat <- mat[mars,mars]
diag(mat) = 0

#mat[lower.tri(mat)] = 0
n = nrow(mat)
rn = rownames(mat)
###############

mar.names <- matrix(mars, nrow = length(mars), ncol = length(mars))
mat.names <- matrix(mars, nrow = length(mars), ncol = length(mars))
##mar.names <- matrix(mars, nrow = dim(rf)[1], ncol = dim(rf)[2])
##mat.names <- matrix(mars, nrow = dim(rf)[1], ncol = dim(rf)[2])
mat.names <- gsub(":.*","",mat.names)

diag(mat.names) = 0
mat.names[lower.tri(mat.names)] = 0
n = nrow(mat.names)

rownames(mat.names) <- rownames(mat)
colnames(mat.names) <- rownames(mat)
rn = rownames(mat.names)

mar_b <- colnames(mat)
mar_a <- rownames(mat)

###############
lod <- lod[mars,mars]
lod <- data.matrix(lod)
rf <- rf[mars,mars]
rf <- data.matrix(rf)
###############

###############
s1 <- scanone(cross1, pheno.col=5, model="normal", method="mr")
s1 <- s1[mars,'lod']

s1 <- matrix(s1, nrow = length(mars), ncol = length(mars))
rownames(s1) <- colnames(s1) <- mars

#s1l[lower.tri(s1l, diag = T)] <- NA
###############

###############
#lod_phen <- norm.mr.2$lod
#lod_phen[lower.tri(lod_phen, diag = T)] <- NA
#diag(lod_phen) = 0
#lod_phen[lower.tri(lod_phen)] = 0

#load(file.path(mpath,paste0(pop,'_scan2_normal_mr.rsave')))

lod_phen <- data.matrix(norm.mr.2$lod)
rownames(lod_phen) <- markernames(cross1)
colnames(lod_phen) <- markernames(cross1)
lod_phen <- lod_phen[mars,mars]
lod_phen <- data.matrix(lod_phen)

###########################################################################

load(file.path(mpath,paste0(pop,'_csq_scan.rsave')))

lod_hom <- data.matrix(-log10(csq_mod.pval[mars,mars]))
lod_inc <- data.matrix(-log10(csq.pval[mars,mars]))

###########################################################################



ab <- lod
phen <- lod_phen
rfs <- rf.df
homz <- lod_hom
inco <- lod_inc
s1 <- s1



ab <- lod[cbind(mar_a, mar_b), drop = T]
phen <- lod_phen[cbind( mar_a, mar_b), drop = T]
rfs <- rf.df[cbind( mar_a, mar_b), drop = T]
homz <- lod_hom[cbind( mar_a, mar_b), drop = T]
inco <- lod_inc[cbind(mar_a, mar_b), drop = T]
s1 <- s1[cbind(mar_a, mar_b), drop = T]

chr_a <- map[mar_a, c('chr','pos')]
chr_b <- map[mar_b, c('chr','pos')]

##links <- data.frame(cbind(chr_a,chr_b,lod_ab,lod_p,rfs),stringsAsFactors=F)

links <- data.frame(cbind(chr_a,chr_b,ab,rfs,homz,inco,phen,s1),stringsAsFactors=F)

try <- data.frame(cbind(lod,rf))

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


chroms <- unique(c(as.character(links[,1]),as.character(links[,3])))
ab_tables <- lapply(chroms, make_lodrf_tables, Y = links, Z = 'lod_ab')
p_tables <- lapply(chroms, make_lodrf_tables, Y = links, Z = 'lod_p')
rf_tables <- lapply(chroms, make_lodrf_tables, Y = links, Z = 'rfs')
hom_tables <- lapply(chroms, make_lodrf_tables, Y = links, Z = 'lod_homz')
dist_tables <- lapply(chroms, make_lodrf_tables, Y = links, Z = 'lod_inco')
names(ab_tables) <- names(p_tables) <- names(rf_tables) <- chroms

################################################################################

lod_hom[lower.tri(lod_hom, diag = T)] <- NA
lod_inc[lower.tri(lod_inc, diag = T)] <- NA
rfs[lower.tri(rfs, diag = T)] <- NA
lod_phen[lower.tri(lod_phen, diag = T)] <- NA
lod[lower.tri(lod, diag = T)] <- NA
rf.df[lower.tri(rf.df, diag = T)] <- NA
s1[lower.tri(s1, diag = T)] <- NA

hoz <- quantile(lod_hom, 0.9999, na.rm=T)
inc <- quantile(lod_inc, 0.9999, na.rm=T)
rfq <- quantile(rfs, 0.999, na.rm=T)
phe <- quantile(lod_phen, 0.999, na.rm=T)
abq <- quantile(lod_ab, 0.999, na.rm=T)

################################################################################
################################################################################

load(file.path(mpath,paste0(pop,'_circos.rsave')))

################################################################################

lod_gtl <- which(links$s1 > 4)
linked <- which(links$homz > hoz | links$inco > inc)
ind <- intersect(lod_gtl,linked)
toplot <- links[ind,]
dim(toplot)

################################################################################

plot_test(paste0(pop,'circ_incompat'))
pc(toplot)
dev.off()

################################################################################
################################################################################

lod_gtl <- which(links$phen > 10)
linked <- which(links$homz > hoz | links$inco > inc)
ind <- intersect(lod_gtl,linked)
toplot <- links[ind,]
dim(toplot)

################################################################################

plot_test(paste0(pop,'circ_incompat'))
pc(toplot)
dev.off()

################################################################################
################################################################################


which(toplot$chr.1 == 15)

##### INCOMPATABILITY WITH CYP IN NBH TOO?
a <- find.marker(cross,15,6.562248)
b <- '2:21870833'
          2:21870833
15:3296707  - AA AB BB
        -   0  0  0  0
        AA  0 11 18  1
        AB  0 10 14 12
        BB  0 13  2 11

a <- '15:4060195'

### INCOMPAT WITH
a <- find.marker(cross,14,63.452152)
b <- '2:24377353'




a <- links[which.max(links$inco),]
b <- find.marker(cross,a$chr.1,a$pos.1)
a <- find.marker(cross,a$chr,a$pos)
geno.crosstab(cross,a,b)


links[which.max(links$phen),]

a <- links[which.max(links$homz),]
a <- find.marker(cross,a$chr.1,a$pos.1)
geno.crosstab(cross, '1:291287',"24:23237312")
### LOD 4



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

################################################################################
################################################################################
################################################################################

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
