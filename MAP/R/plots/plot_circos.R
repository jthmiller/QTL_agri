#load(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
pop <- 'NBH'
pop <- 'ELR'

library('qtl')
library(circlize)
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'

################################################################################
################################################################################
load(file.path(mpath,paste0(pop,'_circos.rsave')))
################################################################################
################################################################################

counts <- table(mtcars$vs, mtcars$gear)
barplot(counts, main="Car Distribution by Gears and VS",
  xlab="Number of Gears", col=c("darkblue","red"),
  legend = rownames(counts))



################################################################################
################################################################################

plot_test(paste0(pop,'_allele_balance2'), width = 4000, height = 2000)
par(mfrow = c(4,6))
for(i in 1:24){
counts <- geno.table(cross,i)[,c('AA','AB','BB')]
barplot(t(counts), col=c("red","purple","darkblue"),space=c(0,0,0,0), border = NA, xlab=NA)
}
dev.off()

################################################################################
################################################################################
plot_test(paste0(pop,'_allele_balance'), width = 4000, height = 2000)
par(mfrow = c(4,6))

for(i in 1:24){

 ind <- which(gt$chr == i)

 x <- map[names(gt.AB[ind]),'pos']

 #par(mfrow = c(2,1), new=TRUE)
 #plot(x, s1a[rownames(map[names(gt.AB[ind]),]),'lod'], ylim = c(0,15), col = 'cornflowerblue', xlab = NA, ylab = NA, cex.axis=3)

 plot(x, gt.AB[ind], ylim = c(-15,15), type='n', col = 'purple', xlab = NA, ylab = NA, cex.axis=3)

 ## xleft, ybottom, xright, ytop
 rect(0, aa_lm_pred['lwr']  , max(x,na.rm=T), aa_lm_pred['upr'], density = 30, col = 'red',border=F)
 rect(0, bb_lm_pred['lwr']  , max(x,na.rm=T), bb_lm_pred['upr'], density = 30, col = 'blue',border=F)
 rect(0, ab_lm_pred['lwr']  , max(x,na.rm=T), ab_lm_pred['upr'], density = 30, col = 'purple', border=F)

 abline(h = 0, col = 'black')

lapply(

 lines(x, predict(yAB), col = 'purple', lwd = 3)
 lines(x, predict(yAA), col = 'red', lwd = 3)
 lines(x, predict(yBB), col = 'blue', lwd = 3)
 text(0,14,i, cex=4)

 if(any(ahr$chr == i)) {
  sapply(which(ahr$chr == i), function(X) {
   abline(v = ahr[X,'pos'], col = 'black')
  })
 }
}
dev.off()
################################################################################
################################################################################


################################################################################
################################################################################
plot_test(paste0(pop,'_999_circ_homz'))
 circos.par("track.height" = 0.15)
 circos.initialize(factors = map$chr, x = map$pos)

## Track of LOD phenotype
 circos.track(factors = map$chr, y = yP,
    panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"),CELL_META$sector.index)
        circos.axis(labels.cex = 0.6)
     })
 circos.trackLines(facP, x = xP, y = yP, pch = 16, cex = 0.5, type = 'I', area = T, col = 'blue')
 ## Track of LOD phenotype

 ## Track of segregation distortion
 circos.track(factor = map$chr, x = map$pos, ylim=c(0,4))
 circos.trackLines(facp, x = xp, y = yp, col = 'purple', pch = 16, cex = 0.5, type = 'I', area = T)
 ## Track of segregation distortion

 ## Connections
 #pc.plot(top_links(links2 = links_ord,n = 500, col = 'lod_hom'), col = 'red')
 #pc.plot(top_links(links2 = links_ord,n = 500, col = 'lod_inc'), col = 'blue')
 ##pc.plot(top_links(links2 = links_ord,n = 5, col = 'lod_phen'), col = 'green')
 pc.plot(quant_links(links2 = links, n = 0.999975, col = lod_hom), col = 'red')
 pc.plot(quant_links(links2 = links, n = 0.99999, col = lod_inc), col = 'blue')

dev.off()
circos.clear()
################################################################################
################################################################################

################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
plot_test('elr18_17')
effectplot(cross1,pheno.col=1,mname1='18:20273448',mname2='17:11276576')
dev.off()

plot_test('elr_18')
plotPXG(cross1,pheno.col=1,'18:20273448')
dev.off()
################################################################################
################################################################################

### SCANONE PLOTS ##############################################################
cross1$pheno$homa <- pull.geno(cross1)[,'24:23951185']
cross1$pheno$incb1 <- as.numeric(pull.geno(cross1)[,'1:291287'] == 1)
cross1$pheno$incb2 <- as.numeric(pull.geno(cross1)[,'1:291287'] == 2)
cross1$pheno$incb3 <- as.numeric(pull.geno(cross1)[,'1:291287'] == 3)
cross1$pheno$inca <- pull.geno(cross1)[,'22:19409957']
cross1$pheno$incb <- pull.geno(cross1)[,'2:31789344']

### MAKE COVARIATES
mar <- '1:291287'
mar <- '2:31789344'
mar <- '24:23951185'

g.int <- pull.geno(fill.geno(cross1))[,mar]
g.int <- data.frame(cbind(as.numeric(g.int == 1), as.numeric(g.int == 2)))

sone.int1 <- scanone(cross1, pheno.col=1, model="normal", method="em")
sone.int2 <- scanone(cross1, pheno.col=1, model="normal", method="em", intcovar = as.matrix(g.int))
sone.int3 <- scanone(cross1, pheno.col=1, model="normal", method="em", addcovar = as.matrix(g.int))
################################################################################


chrs <- c(1:23)
plot_test('segdist_s1_int',width=1000)
plot(sone.int2, col='red', chr = chrs)
plot(sone.int1, col='green', chr = chrs, add=T)
plot(sone.int3, col='blue', chr = chrs, add=T)
dev.off()


################################################################################

segdist_int_1 <- scan1(pr2, pheno=cross2$pheno[,'bin'] , model = "binary", cores = cores, intcovar = cross2$pheno[,'incb1'])
segdist_int_2 <- scan1(pr2, pheno=cross2$pheno[,'bin'] , model = "binary", cores = cores, intcovar = cross2$pheno[,'incb2'])
segdist_int_3 <- scan1(pr2, pheno=cross2$pheno[,'bin'] , model = "binary", cores = cores, addcovar = cross2$pheno[,'incb3'])

segdist_homa <- scan1(pr2, pheno=cross2$pheno[,'homa'], model = "normal", cores = cores)
segdist_homb <- scan1(pr2, pheno=cross2$pheno[,'homb'], model = "normal", cores = cores)
segdist_inca <- scan1(pr2, pheno=cross2$pheno[,'inca'], model = "normal", cores = cores)
segdist_incb <- scan1(pr2, pheno=cross2$pheno[,'incb'], model = "normal", cores = cores)

plot_test('segdist_int',width=1000)
plot(segdist_int_1, map2, col='red')
plot(segdist_int_2, map2, col='green')
plot(segdist_int_3, map2, col='blue')
dev.off()

################################################################################
################################################################################

hoz <- quantile(lod_hom, 0.999, na.rm=T)
inc <- quantile(lod_inc, 0.999, na.rm=T)
rfq <- quantile(rf, 0.999, na.rm=T)
phe <- quantile(lod_phen, 0.999, na.rm=T)
abq <- quantile(lod, 0.999, na.rm=T)
s1q <- quantile(s1, 0.95, na.rm=T)

lod_gtl <- which(links$s1 > 2)
linked <- which(links$lod_hom > hoz | links$lod_inc > inc)
ind <- intersect(lod_gtl,linked)

###########################################################################
################################################################################

################################################################################
################################################################################
csq.pval.hm <- data.matrix(-log10(csq_mod.pval))

for (i in chrnames(cross)){
 mars2 <- markernames(cross, i)
 ##mars2 <- mars2[which(mars2 %in% mars)]
 csq.pval.hm[mars2, mars2] <- NA
}

csq.pval.hm[lower.tri(csq.pval.hm, diag = T)] <- NA

plot_test('heatmap_homozygotes_elr',height=3000,width=3000)
heatmap(csq.pval.hm, Rowv=NA, Colv=NA, scale="column")
dev.off()
################################################################################
################################################################################

################################################################################
################################################################################
top_links(links2 = links_ord,n = 200, col = 'lod_inc')
top_links(links2 = links_ord,n = 200, col = 'lod')
top_links(links2 = links_ord,n = 200, col = 'lod_hom')
top_links(links2 = links_ord,n = 5, col = 'lod_phen'
top_links(links2 = links_ord,n = 300, col = 's1')
top_links(links2 = links_ord,n = 300, col = 'rf')
################################################################################
################################################################################

################################################################################
################################################################################
geno.crosstab(cross,'18:20273448','17:11276576')
geno.crosstab(cross,'18:19222949','13:20651592')
geno.crosstab(cross,'18:17874376','7:35475725')
geno.crosstab(cross,'24:23237312','1:291287')
geno.crosstab(cross,'24:10521322','1:291287')
geno.crosstab(cross,'17:29388433','2:13880632')
geno.crosstab(cross,'18:20771373','17:14629450')
geno.crosstab(cross,'8:10228388','2:30305156')
geno.crosstab(cross,'18:19222949','2:32379535')
geno.crosstab(cross,'18:20010144','17:12700256')
geno.crosstab(cross,'24:27146556','1:291287')
geno.crosstab(cross,'18:20745539','17:14629450')
geno.crosstab(cross,'17:7480177','8:37635736')
geno.crosstab(cross,'17:33155007','5:1282903')
geno.crosstab(cross,'17:29007925','5:1896342')
geno.crosstab(cross,'24:23951185','1:291287')
geno.crosstab(cross,'19:37446878','1:33855817')
geno.crosstab(cross,'1:857165','24:23363287')
geno.crosstab(cross1,'1:857165','24:24242668')
geno.crosstab(cross1,'1:857165','24:23951185')
geno.crosstab(cross1,'13:22410641','7:946994')
geno.crosstab(cross1,'24:23951185','1:291287')
geno.crosstab(cross1,'5:33290434','1:291287')
geno.crosstab(cross1,'13:25406978','24:23951185')
geno.crosstab(cross,'24:24013348','1:291287')
geno.crosstab(cross,'24:24013348','1:291287')
geno.crosstab(cross,'13:22410641','7:946994')
geno.crosstab(cross,'24:30529145','1:291287')
geno.crosstab(cross,'19:37446878','1:33855817')
geno.crosstab(cross,'22:19409957','2:31789344')
geno.crosstab(cross1,'18:20509237','17:12700256')


HSP and AIP
geno.crosstab(cross1,'22:19528880','2:27373969')
geno.crosstab(cross1,'22:19409957','2:32539238')


geno.crosstab(subset(cross1, ind=cross1$pheno$bin == 0),'14:1738185','24:23951185')
geno.crosstab(subset(cross1, ind=cross1$pheno$bin == 1),'14:1738185','24:23951185')
geno.crosstab(subset(cross1, ind=cross1$pheno$bin == 1),'13:25406978','24:23951185')
geno.crosstab(subset(cross1, ind=cross1$pheno$bin == 0),'13:25406978','24:23951185')
geno.crosstab(subset(cross1, ind=cross1$pheno$bin == 0),'5:33290434','1:291287')
geno.crosstab(subset(cross1, ind=cross1$pheno$bin == 1),'5:33290434','1:291287')
geno.crosstab(subset(cross1, ind=cross1$pheno$bin == 0),'1:857165','24:23951185')
geno.crosstab(subset(cross1, ind=cross1$pheno$bin == 1),'1:857165','24:23951185')

################################################################################
################################################################################

55   13   hspa14  60   1 23552848   81768  13:23471080  HSP
56   13   tmtc2a  34   1 10549074   37468  13:10586542 <NA>
57   13     ARNT  61   2 24528101   71711  13:24456390
28    7    hspb7  62   0 30894102   97788   7:30991891  HSP
29    7    rnd1a  22   1  9736976  198968    7:9538008 <NA>
30    7    rnd1b  48   0 21809146   19110   7:21790036 <NA>
31    7     vdra  66   0 33273622     976   7:33274599 <NA>
32    7    fkbpl  41   0 19289922   56853   7:19346775 <NA>
33    8   atxn1b  56   2 21146220   13720   8:21132501  AHR
1     1     AHR1   0   2   741331  115834     1:857165  AHR

96393   24:23363287  1:5578105 1.956 0.644   4.775   3.079    2.726 1.884
################################################################################
################################################################################

################################################################################
################################################################################
### Incompat in
geno.crosstab(cross,'18:17874376','7:35475725')
intxs.bin('7:35475725','18:17874376',  popchr = "18v7", locbN = 'test', main = 'test2')

geno.crosstab(cross,'8:10228388','2:30305156')
intxs.bin('8:10228388','2:30305156',  popchr = "2v8", locbN = 'test', main = 'test2')

###################################rfgterD######################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

##### INCOMPATABILITY WITH CYP IN NBH TOO?
a <- find.marker(cross,15,6.562248)
b <- '2:21870833'
          2:21870833
15:3296707  - AA AB BB
        -   0  0  0  0
        AA  0 11 18  1
        AB  0 10 14 12
        BB  0 13  2 11
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
