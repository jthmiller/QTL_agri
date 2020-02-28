#!/bin/R
pop <- 'ELR'
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

load(file.path(mpath,paste0(pop,'_scan2_bin_em_noCof.rsave')))
################################################################################
#sone.o <- scanone(cross,pheno.col=4, model="binary", method="em")
#sone.a <- scanone(cross,pheno.col=4, model="binary", method="em", addcovar=g[,1])
#sone.i <- scanone(cross,pheno.col=4, model="binary", method="em", addcovar=g[,1],intcovar=g[,1])
#sone.io <- scanone(cross,pheno.col = 4, model="binary", method="em", addcovar=g[,1],intcovar=g[,1])
#cbind(summary(sone.o),summary(sone.a)$lod,summary(sone.i)$lod,summary(sone.io)$lod)

rf <- subset(cross, chr = c(1:4,6:24))
rf <- est.rf(rf, maxit=100000, tol=1e-6)

s1 <- scanone(rf,pheno.col=4, model="binary", method="em")
s1l <- matrix(s1$lod, nrow = dim(rf.df)[1], ncol = dim(rf.df)[1])

mars <- find.marker(rf, bin.em.2$map$chr, bin.em.2$map$pos)
mat.names <- matrix(mars, nrow = dim(rf.df)[1], ncol = dim(rf.df)[2])
mat.names <- gsub(":.*","",mat.names)

rf.df <- pull.rf(rf)
rf.df <- rf.df[mars,mars]

lod.df <- pull.rf(rf, what='lod')
lod.df <- lod.df [mars,mars]

lod_phen <- bin.em.2$lod

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


###########################################################################
chr1 <- gsub(":.*","",colnames(rf.df)) %in% c(1)
chr18 <- gsub(":.*","",colnames(rf.df)) %in% c(18)
chr2 <- gsub(":.*","",colnames(rf.df)) %in% c(2)
#########################


col <- matrix('grey', nrow = dim(rf.df)[1], ncol = dim(rf.df)[1])

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
####################################################################################################





############################################
#plot
############################################
no_qtl_mr <- scanone(cross, pheno.col=4, method="mr", model="binary")
qtl <- summary(no_qtl_mr, 4)

for (i in chrnames(rf)){
 ind <- which(markernames(rf) %in% markernames(rf,i))
 rf$rf[ind,ind] <- NA
}

quantile(pull.rf(rf), 0.999,na.rm=T)

rf.df[lower.tri(rf.df, diag = T)] <- NA

h <- quantile(rf.df, 0.999,na.rm=T)
l <- quantile(rf.df, 0.001,na.rm=T)
h9 <- quantile(rf.df, 0.95,na.rm=T)
l9 <- quantile(rf.df, 0.05,na.rm=T)

plot_test('elr_1_2_8_13_18', width=1250,height=1000)
par(mfrow = c(6,1))

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

dev.off()

################################################################################
################################################################################

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


mat.names <- matrix(mars, nrow = dim(rf.df)[1], ncol = dim(rf.df)[1])
mat.names <- gsub(":.*","",mat.names)


ind1 <- gsub(":.*","",colnames(rf.df)) %in% c(2)
ind2 <- gsub(":.*","",colnames(rf.df)) %in% c(18)

col <- matrix('black', nrow = dim(rf.df)[1], ncol = dim(rf.df)[1])
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




ar_ind <- which(lod_phen > 15 & rfm < 0.45, arr.ind = TRUE)


table(mat.names[which(lod_phen > 20 & rfm < 0.5 & rfm > 0.4, arr.ind = TRUE)])


summary(lod_phen[which(lod_phen > 20 & rfm < 0.5, arr.ind = TRUE)])



sort(table(mat.names[which(lod_phen > 15 & rfm < 0.45 | rfm > 0.65, arr.ind = TRUE)]))

tf.mat <- matrix(lod_phen > 15 & rfm < 0.45, nrow = 1618, ncol = 1618)
tf.mat <- matrix(lod_phen > 10 & rfm < 0.45, nrow = 1618, ncol = 1618)


tf.mat <- matrix(lod_phen < 15 & rfm > 0.65 | rfm < 0.4, nrow = 1618, ncol = 1618)

which(lod_phen == max(lod_phen, na.rm=T), arr.ind = TRUE)

sort(table(gsub(":.*","",markernames(rf)[rowSums(tf.mat, na.rm = T) > 0])))


bin.em.2$lod,


 16   21    3


################################################################################
cross <- subset(cross, chr=c(1:4,6:24))

no_qtl <- scanone(cross, pheno.col=4, method="hk", model="binary", verbose=FALSE, tol=1e-4, maxit=1000)

add.perms <- scanone(cross, pheno.col=4, method="hk", model="binary", n.perm = 1000, n.cluster=6, perm.Xsp=T)

lod <- summary(add.perms)[[1]][1]

qtl <- summary(no_qtl,lod)

cross$pheno <- as.data.frame(cross$pheno)

  add.qtl1 <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")
  add.qtl1 <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="draws")

add.qtl1 <- refineqtl(cross, pheno.col = 4, qtl=add.qtl1, method = "hk", model='binary', incl.markers=F)
     int.em <- addint(cross, pheno.col = 4, qtl = add.qtl1, method='hk', model='binary', covar=data.frame(cross$pheno$sex) ,formula=y~Q1+Q2+Q3, maxit=10000)

add_Q4 <- addqtl(cross, pheno.col=4, qtl = add.qtl1, method="hk", model="binary", formula = y~Q1+Q2+Q3+Q4,
            incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000)

#add_Q5_wInts <- addqtl(cross, pheno.col=4, qtl = add.qtl1, method="hk", model="binary",
#            incl.markers=T, verbose=FALSE, tol=1e-4, maxit=1000,
#            formula = y~Q1*Q3+Q2+Q4)
#
#int.em <- addint(cross, qtl = add.qtl1, formula=y~Q1+Q2+Q3+Q4, method='hk', model='binary')

################################################################################

fit_3 <- fitqtl(cross, pheno.col=4, method="hk", model="binary", qtl = add.qtl1,
  covar=NULL, formula=y~Q1+Q2+Q3+Q1:Q3, dropone=TRUE, get.ests=T,
  run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)

################################################################################



plot_test('nbh_add_Q5_add_Q5_wInts.png',width=1000)
 plot(no_qtl, add_Q5, add_Q5_wInts)
dev.off()

plot_test('nbh_add_Q5_add_Q5_wInts.png',width=1000)
 plot(no_qtl, add_Q5, add_Q5_wInts)
dev.off()

add_Q3Q1_int <- addqtl(cross, pheno.col=4, qtl = add.qtl1, method="hk", model="binary",
            incl.markers=F, verbose=FALSE, tol=1e-4, maxit=1000,
            formula = y~Q1*Q2+Q1*Q3)

add_Q3Q1_int <- addqtl(cross, pheno.col=4, qtl = add.qtl1, method="hk", model="binary",
            incl.markers=F, verbose=FALSE, tol=1e-4, maxit=1000,
            formula = y~Q1*Q2+Q2*Q3)

plot_test('int_Q5.png',width=1000)
 plot(add_Q3, add_Q3Q1_int, add_Q3Q1_int)
dev.off()

plot_test('add_Q5.png',width=1000)
 plot(no_qtl, add_Q3)
dev.off()
################################################################################

################################################################################
no_qtl <- scanone(cross, pheno.col=4, method="hk", model="binary", verbose=FALSE, tol=1e-4, maxit=1000)

AHR.qtl <- makeqtl(cross, chr=1, pos=0, what="prob")

add_Q2_AHR <- addqtl(cross, pheno.col=4, qtl = AHR.qtl, method="hk", model="binary",
            incl.markers=F, verbose=FALSE, tol=1e-4, maxit=1000,
            formula = y~Q1+Q2)

int_Q2_AHR <- addqtl(cross, pheno.col=4, qtl = AHR.qtl, method="hk", model="binary",
            incl.markers=F, verbose=FALSE, tol=1e-4, maxit=1000,
            formula = y~Q1*Q2)

plot_test('only_AHR_qtl.png',width=1000)
 plot(no_qtl, add_Q2_AHR, int_Q2_AHR)
dev.off()
################################################################################


add.perms <- scanone(cross, pheno.col=4, model='binary', method = "hk", n.perm = 1000, n.cluster=6)
lod <- summary(add.perms)[2]
add <- scanone(cross, pheno.col=4, model='binary', method = "hk")
qtl <- summary(add,lod)

add.qtl1 <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")
add.qtl1 <- refineqtl(cross, qtl=add.qtl1, pheno.col=4, model='binary', method = "hk", incl.markers=F)

int.em <- addint(cross, qtl=add.qtl1, formula=y~Q1+Q2+Q3+Q4, method='hk')

bin.add.em.qtls2 <- refineqtl(cross, pheno.col=4, model='binary',
   qtl=add.qtl1, method='hk', incl.markers=F,
   formula=y~Q1+Q2+Q3+Q4+Q1:Q3)





#### EM ##################

no_qtl_em <- scanone(cross, pheno.col=4, method="em", model="binary", maxit=1000)

add.perms <- scanone(cross, pheno.col=4, method="em", model="binary", n.perm = 1000, n.cluster=6)

lod <- summary(add.perms)[1]

qtl <- summary(no_qtl_em,lod)

Q3 <- makeqtl(cross, chr=qtl[['chr']][c(1,3)], pos=qtl[['pos']][c(1,3)], what="draws")

cross <- sim.geno(cross, stepwidth="fixed", step=1,  error.prob=erp, off.end=1, map.function="kosambi", n.draws=100)

fit_3_em <- fitqtl(cross, pheno.col=4, method="imp", model="binary", qtl = Q3,
  covar=NULL, formula=y~Q1+Q2, dropone=TRUE, get.ests=T,
  run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)

fit_3_em_int <- fitqtl(cross, pheno.col=4, method="imp", model="binary", qtl = Q3,
  covar=NULL, formula=y~Q1+Q2+Q3+Q1:Q3, dropone=TRUE, get.ests=T,
  run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)



##### IMP ########
##### To fit imp, must use em to scan for first QTL  (not offered rQTL)
no_qtl_im <- scanone(cross, pheno.col=4, method="imp", model="binary")

imp_perms <- scanone(cross, pheno.col=4, method="imp", model="binary", n.perm = 10000, n.cluster=6)

lod <- summary(imp_perms)[1]

qtl <- summary(no_qtl_im, lod)

Q3 <- makeqtl(cross, chr=qtl[['chr']], pos=qtl[['pos']], what="draws")

Q3 <- refineqtl(cross, pheno.col = 4, qtl=Q3, method = "imp", model='binary',
                incl.markers=F)

int.imp <- addint(cross, pheno.col = 4, qtl = Q3, method='imp', model='binary',
                  covar=data.frame(cross$pheno$sex) ,formula=y~Q1+Q2+Q3, maxit=10000)

fit_3_imp <- fitqtl(cross, pheno.col=4, method="imp", model="binary", qtl = Q3,
                    formula=y~Q1+Q2+Q3, dropone=TRUE, get.ests=T, covar=data.frame(cross$pheno$sex),
                    run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)

scan_3_imp <- scanqtl(cross, pheno.col=4, method="imp", model="binary", qtl = Q3,
                      covar=data.frame(cross$pheno$sex), formula=y~Q1+Q2+Q3,
                      incl.markers=FALSE, verbose=TRUE, tol=1e-4, maxit=1000,
                      forceXcovar=FALSE)



pens <- calc.penalties(perms.2, alpha=0.05)

summary(bin.em.2, perms=perms.2, alphas=0.05, pvalues=F)

summary(bin.em.2, perms=perms.2, alphas=0.3, pvalues=F, what='int')

summary(bin.em.2, thresholds=c(0, Inf, 4, Inf, Inf), what="int")


summary(bin.em.2,perms=perms.2,alphas=0.1, pvalues=T)

summary(bin.em.2,perms=perms,alphas=0.1, pvalues=F)

summary(bin.em.2,perms=perms,alphas=0.2, pvalues=F, what='int')
summary(bin.em.2,perms=perms,alphas=0.05, pvalues=F, what='int')

summary(bin.em.2,perms=perms,alphas=0.1, pvalues=T)

############################################################
##pens
##    main    heavy    light
##3.971724 6.831380 3.617315
## 2nd time
#    main    heavy    light
#3.595578 6.820764 4.083055
############################################################
pop <- 'NBH'
load(file.path(mpath,paste0(pop,'_norm_imp.rsave')))
############################################################
## NBH
summary(full.norm.imp)

      name chr pos n.gen
Q1  2@43.0   2  43     3
Q2  2@87.0   2  87     3
Q3  3@31.0   3  31     3
Q4  3@38.0   3  38     3
Q5 13@31.0  13  31     3
Q6 18@51.0  18  51     3
Q7 19@18.0  19  18     3

Formula: y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q3:Q5 + Q1:Q6

pLOD:  49.014

fitqtl(gg_step2, pheno.col=5, qtl=full.norm.imp, method="imp",model="normal",get.ests=T,covar=data.frame(gg_step2$pheno$sex),
 formula = "y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q3:Q5 + Q1:Q6")

############################################################

$result.full
      df         SS          MS      LOD     %var Pvalue(Chi2) Pvalue(F)
Model 22 469.190321 21.32683275 84.05107 98.57858            0         0
Error 68   6.765348  0.09949041       NA       NA           NA        NA
Total 90 475.955669          NA       NA       NA           NA        NA

############################################################

              df Type III SS       LOD       %var     F value Pvalue(Chi2)
2@43.0          6   32.610873 34.804947  6.8516619   54.629841 0.000000e+00
2@87.0          2  270.988256 73.408315 56.9356084 1361.881232 0.000000e+00
3@31.0          6   50.181181 42.095558 10.5432468   84.063678 0.000000e+00
3@38.0          2    1.939054  4.979921  0.4074023    9.744931 1.047319e-05
13@31.0         6   48.466038 41.491258 10.1828891   81.190465 0.000000e+00
18@51.0         6   78.526263 50.078036 16.4986506  131.547452 0.000000e+00
19@18.0         2    5.337543 11.493253  1.1214370   26.824407 3.211764e-12
2@43.0:18@51.0  4   28.746096 32.763557  6.0396583   72.233331 0.000000e+00
3@31.0:13@31.0  4   49.014383 41.686475 10.2980984  123.163582 0.000000e+00

############################################################

plot_test('nbh_full.norm.imp',width=1500)
plotLodProfile(full.norm.imp,incl.markers=F)
dev.off()

############################################################

summary(full.norm.hk)

 name chr pos n.gen
Q1  2@103.0   2 103     3
Q2   3@17.0   3  17     3
Q3   3@84.0   3  84     3
Q4  11@54.0  11  54     3
Q5  13@28.0  13  28     3
Q6 15@114.0  15 114     3
Q7  18@41.0  18  41     3
Q8  19@49.0  19  49     3

  Formula: y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q7 + Q2:Q5 + Q4:Q8 +
                Q3:Q6

  pLOD:  30.266

fitqtl(gg_step2, pheno.col=5, qtl=full.norm.hk, method="hk",model="normal",get.ests=T,covar=data.frame(gg_step2$pheno$sex),
 formula = "y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q1:Q7 + Q2:Q5 + Q4:Q8 + Q3:Q6")

$result.full
      df        SS         MS      LOD     %var Pvalue(Chi2) Pvalue(F)
Model 32 461.60114 14.4250356 69.18628 96.98406            0         0
Error 58  14.35453  0.2474919       NA       NA           NA        NA
Total 90 475.95567         NA       NA       NA           NA        NA

                df Type III SS       LOD       %var    F value Pvalue(Chi2)
2@103.0          6  185.200867 52.009888 38.9113691 124.718463 0.000000e+00
3@17.0           6   27.285307 21.044661  5.7327413  18.374544 0.000000e+00
3@84.0           6   11.308472 11.480490  2.3759507   7.615381 1.246408e-09
11@54.0          6   13.070904 12.792985  2.7462439   8.802243 7.478684e-11
13@28.0          6   31.630736 23.006152  6.6457316  21.300855 0.000000e+00
15@114.0         6    7.507343  8.312777  1.5773199   5.055615 9.895083e-07
18@41.0          6   43.708858 27.614558  9.1833885  29.434536 0.000000e+00
19@49.0          6   10.667726 10.980855  2.2413277   7.183889 3.615243e-09
2@103.0:18@41.0  4   23.860577 19.348697  5.0131932  24.102381 0.000000e+00
3@17.0:13@28.0   4   20.881997 17.745186  4.3873827  21.093615 1.110223e-16
3@84.0:15@114.0  4    7.185182  8.019416  1.5096327   7.257997 1.861430e-07
11@54.0:19@49.0  4    2.120934  2.723114  0.4456158   2.142427 1.375410e-02



nor_imp_fit <- fitqtl(gg_step2, pheno.col=5, qtl=full.norm.imp, method="imp",model="normal",get.ests=T,covar=data.frame(gg_step2$pheno$sex),
 formula = "y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q3:Q5 + Q1:Q6")

bin_imp_fit <- fitqtl(gg_step2, pheno.col=4, qtl=full.norm.imp, method="imp",model="binary", get.ests=F,covar=data.frame(gg_step2$pheno$sex),
 formula = "y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q3:Q5 + Q1:Q6")


### DROP QTL
############################################################
qtl_drop <- dropfromqtl(full.norm.hk,qtl.name=c('11@54.0','19@49.0','15@114.0','3@84.0'))
qtl_drop_fit <- fitqtl(gg_step2, pheno.col=5, qtl=qtl_drop, method="hk",model="normal",get.ests=T,covar=data.frame(gg_step2$pheno$sex),
 formula = "y ~ Q1 + Q2 + Q3 + Q4 + Q1:Q4 + Q2:Q3")
qtl_drop <- refineqtl(gg_step2,qtl=qtl_drop,keeplodprofile=TRUE)
plot_test('nbh_full.norm.imp',width=1500)
plotLodProfile(qtl_drop,incl.markers=F)
dev.off()
############################################################

## manual add stepwise qtl
bin.add.em.perms <- scanone(gg_step2, pheno.col=4, model='binary', method = "hk", n.perm = 10000, n.cluster=10)
lod <- summary(bin.add.em.perms)[2]
bin.add.em <- scanone(gg_step2, pheno.col=4, model='binary', method = "hk")
qtl <- summary(bin.add.em,lod)
bin.add.em.qtls1 <- makeqtl(gg_step2, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")
bin.add.em.qtls1 <- refineqtl(gg_step2, qtl=bin.add.em.qtls1, pheno.col=4, model='binary', method = "hk", incl.markers=F)

int.em <- addint(gg_step2, qtl=bin.add.em.qtls1, formula=y~Q1+Q2+Q3, method='hk')
bin.add.em.qtls2 <- refineqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls1, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3)
#int.em <- addint(gg_step2, qtl=bin.add.em.qtls, formula=y~Q1+Q2+Q3+Q1:Q3, method='hk')
##scan for an additional additive QTL
### add.em.a <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls2, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3)
### no additional additive
## scan for an additional interactive QTL
int.em <- addint(gg_step2, qtl=bin.add.em.qtls2, formula=y~Q1+Q2+Q3+Q1:Q3, method='hk')
bin.add.em.qtls3 <- refineqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls2, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q2:Q3)
## scan for additional additive
add.em <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls3, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q2:Q3)
### no additional
add.em.1 <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls3, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q1:Q4)
add.em.2 <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls3, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q2:Q4)
add.em.3 <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls3, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q3:Q4)



## Scan for interacting pair to add (long)
add.em.ap <- addpair(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls3, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q2:Q3)
save.image(file.path(mpath,paste0(pop,'_models.rsave')))

##add.em <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3)
##add.Z <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q4)
##qtl <- rbind(qtl,summary(add.em,lod))
##bin.add.em.qtls <- makeqtl(gg_step2, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")
##bin.add.em.qtls <- refineqtl(gg_step2, pheno.col=4, qtl=bin.add.em.qtls,  model='binary', method = "hk", incl.markers=F)
##
###### No more additive
###add.em <- addqtl(gg_step2,pheno.col=4, qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q4+Q5)
###summary(add.em,lod)
##################################################################################
##int.em <- addint(gg_step2, qtl=bin.add.em.qtls, formula=y~Q1+Q2+Q3+Q4, method='hk')
##bin.add.em.qtls <- refineqtl(gg_step2, pheno.col=4, qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q4+Q3:Q4)
##qtl <- summary(bin.add.em.qtls)
##################################################################################

add.em <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q4+Q3:Q4)
qtl <- summary(add.em,lod)
bin.add.em.qtls <- addtoqtl(gg_step2, qtl=bin.add.em.qtls,chr=qtl$chr, pos=qtl$pos)
bin.add.em.qtls <- refineqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q4+Q5+Q3:Q4)

add.em <- addqtl(gg_step2,pheno.col=4, qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q4+Q5+Q3:Q4)
qtl <- summary(add.em,lod)
bin.add.em.qtls <- addtoqtl(gg_step2, qtl=bin.add.em.qtls,chr=qtl$chr, pos=qtl$pos)
bin.add.em.qtls <- refineqtl(gg_step2, pheno.col=4, qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q4+Q5+Q6+Q3:Q4)

### LESS conservative ############################################
add.em <- addqtl(gg_step2, pheno.col=4, qtl=bin.add.em.qtls_0.05, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q4+Q5+Q6+Q3:Q4)
qtl <- summary(add.em, lod)
int.em <- addint(gg_step2, qtl=bin.add.em.qtls, formula=y~Q1+Q2+Q3+Q4+Q5+Q3:Q4, method='hk')
#### Still none
### Get rid of extra chr8 QTL
 ##out.fq <- fitqtl(gg_step2, pheno.col=4,method='hk', qtl=bin.add.em.qtls_0.05,model='binary',formula=y~Q1+Q2+Q3+Q4+Q5+Q6+Q3:Q4)




out.ap <- addpair(hyper, qtl=rqtl, chr=1, formula=y~Q2+Q3*Q4, verbose=FALSE)
