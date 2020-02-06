library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

load(file.path(mpath,paste0(pop,'_scan2_bin_em_noCof.rsave')))
## with coefload(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
library(circlize)

###############
AHR.bed <- read.table(file.path(mpath,"lift_AHR_genes.bed"), stringsAsFactors = F, header = F)
colnames(AHR.bed) <- c("chrom", "str", "stp", "gene")
AHR.bed$chrom <- as.numeric(gsub("chr", "", AHR.bed$chrom))
AHR.bed$str <- as.numeric(AHR.bed$str)
AHR.bed$stp <- as.numeric(AHR.bed$stp)
AHR.notmap <- AHR.bed[is.na(AHR.bed$chrom), ]
AHR.bed <- AHR.bed[!is.na(AHR.bed$chrom), ]
AHR.bed$gene <- gsub(":158640", "", AHR.bed$gene)
AHR.bed <- AHR.bed[!AHR.bed$chr == 5,]
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
nbh_gens <- cnv.ahrs(rf, AHRdf = AHR.bed, EXP = F)
###############

###############
rf <- subset(cross, chr = c(1:4,6:24))
rf <- est.rf(rf, maxit=100000, tol=1e-6)
mars <- find.marker(rf, bin.em.2$map$chr, bin.em.2$map$pos)
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
########################################

colnames(csq.pval) <- rownames(csq.pval) <- colnames(rf.gts)

for (i in unique(bin.em.2$map$chr)){
 mars <- markernames(rf, i)
 csq.pval[mars,mars] <- NA
}

csq.pval[lower.tri(csq.pval)] = NA

##########################################################################################
#########################################################################################################
#########################################################################################################

##save.image(file.path(mpath,paste0(pop,'_csq_pval.rsave')))
load(file.path(mpath,paste0(pop,'_csq_pval.rsave')))

##########################################################################################
##########################################################################################
##########################################################################################


maxdist <- sapply(unique(bin.em.2$map$chr), function(i) {
 mars <- markernames(rf, i)
 a <- which(csq.pval[mars,] == min(csq.pval[mars,], na.rm = T), arr.ind=T)
 b <- markernames(rf)[as.numeric(a[,'col'])]
 cbind(a, b, -log10(csq.pval[cbind(a[,'row'],a[,'col'])]))
})
maxdist <- do.call(rbind,maxdist)
maxdist <- maxdist[order(as.numeric(maxdist[,4])),]


quantile(-log10(csq.pval),0.999, na.rm = T)

mar <- find.marker(rf,1,10)
mar <- '19:42466531'
head(-log10(sort(csq.pval[,"24:6682975"])), 100)



################################################################################
sb1 <- scanone(rf,pheno.col=4,method="imp",model="bin")
sn1 <- scanone(rf,pheno.col=4,method="imp",model="normal")

col <- gsub(":.*","",markernames(rf))

gt.table <- geno.table(rf)

a <- find.marker(rf,2,88)
a <- '24:6682975'

plot_test('pvals', width = 1500, height = 1000)
par(mfrow= c(2,1))

 plot(1:length(csq.pval[a,]), -log10(csq.pval[a,]), pch = 19, col = NA, ylim = c(0,7))

  points(1:length(csq.pval[a,]),  -log10(gt.table[,'P.value']), pch = 19, col = 'grey', cex = 0.5)

  points(1:length(csq.pval[a,]), -log10(csq.pval[a,]), pch = 19, col = factor(col))

  abline(h=quantile(-log10(csq.pval),0.999, na.rm = T), col = 'red')

 plot(1:length(sn1$lod), sn1$lod, pch = 19, col = 'lightgrey', ylim = c(0,18), cex = 0.5)

 points(1:length(sb1$lod), sb1$lod, pch = 19, col = factor(sb1$chr))

dev.off()

plot_test('pval_hist', width = 1500)
hist(-log10(csq.pval), breaks = 50)
dev.off()


which(sn1$lod == Inf)

## greatest epi-distortion
            row  col
24:7289510 1553  161
2:35090193  161 1553


plot_test('2_13_elr')
effectplot(rf, pheno.col=4, mname1 = find.marker(rf, 13, 63) ,mname2 = find.marker(rf, 2, 77.8), ylim=c(0,1))
dev.off()


X <- find.marker(rf,1, 10)
Y <- find.marker(rf,1, 10)

crstb <- function(X,Y) {
print(geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = X, mname2 = Y))
print(geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = X, mname2 = Y))
print(geno.crosstab(rf, mname1 = X, mname2 = Y))
}


crstb(X = '24:7289510', Y = '2:35090193')

crstb(X = "19:42466531" , Y = '1:33916475')

 summary(pull.map(rf))
        n.mar length ave.spacing max.spacing
1          93   99.5         1.1         5.0
2          82  104.9         1.3        10.6
3          76   87.4         1.2         7.3
4          83   93.1         1.1         9.1
6          43   67.5         1.6        10.5
7          74   85.1         1.2         3.6
8          83  130.4         1.6         8.4
9          49   66.7         1.4         5.6
10         39   39.2         1.0         4.5
11         39   56.0         1.5         9.3
12         62   74.3         1.2         5.0
13         60   76.9         1.3        11.0
14         71   80.8         1.2        14.0
15         93  109.8         1.2         7.7
16         70   73.2         1.1         5.0
17         56   62.8         1.1         7.4
18         43   73.5         1.8        11.1
19         92   91.9         1.0         5.7
20         76   80.1         1.1         5.0
21         86  104.7         1.2        11.2
22         70   70.5         1.0         5.6
23         89  100.1         1.1         4.5
24         89   82.6         0.9         5.7
overall  1618 1911.0         1.2        14.0


#####
Deletion on 24?
6394135  6410483

a <- find.marker(rf,1,1)


24:6682975   24       0  7 63 21      0      0 0.0001385034
crstb(X = '24:6682975', Y = '2:35090193')

find.marker(rf,24,54)
"1:3565012" highest lod_ab (only at 1.4475842)
20:36134345 lowest

head(-log10(sort(csq.pval[,"24:6682975"])), 100)


crstb(X = '24:6682975', Y = '20:36134345')

24:6446937  24      17 18 17 46      0      0 7.490393e-11
