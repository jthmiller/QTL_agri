#!/bin/R

library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)


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

############################################################
qtl_drop <- dropfromqtl(full.norm.hk,qtl.name=c('11@54.0','19@49.0','15@114.0','3@84.0'))
qtl_drop_fit <- fitqtl(gg_step2, pheno.col=5, qtl=qtl_drop, method="hk",model="normal",get.ests=T,covar=data.frame(gg_step2$pheno$sex),
 formula = "y ~ Q1 + Q2 + Q3 + Q4 + Q1:Q4 + Q2:Q3")
qtl_drop <- refineqtl(gg_step2,qtl=qtl_drop,keeplodprofile=TRUE)
plot_test('nbh_full.norm.imp',width=1500)
plotLodProfile(qtl_drop,incl.markers=F)
dev.off()
############################################################

gg_marks <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross)[[X]], 0.75)} ))
gg <- pull.markers(cross,gg_marks)



if(pop == 'ELR.missing') gg_marks <- c(gg_marks,"AHR2a_del")
#

cAZross2 <- sim.geno(cross2, stepwidth="fixed", step=1,off.end=5, error.prob=erp ,map.function="kosambi", n.draws=1000)
na <- stepwiseqtl(cross2, incl.markers=T, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=8)



out.aqi <- addqtl(hyper, qtl=rqtl, formula=y~Q1+Q2+Q3*Q4+Q4*Q5)
out.aq <- addqtl(hyper, qtl=rqtl, formula=y~Q1+Q2+Q3*Q4)

Figure 4: LOD curves for adding one QTL, interacting with the chromosome 15 locus, to the 4-QTL model, with the hyper data.
out.ap <- addpair(hyper, qtl=rqtl, chr=1, formula=y~Q2+Q3*Q4, verbose=FALSE)

addtoqtl, dropfromqtl, and replaceqtl
