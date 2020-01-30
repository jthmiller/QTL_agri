#!/bin/R
pop <- 'NBH'
library('qtl')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

load(file.path(mpath,paste0(pop,1,'_scan_perms_bin_em.rsave')))
perms.2 <- get(paste0('bin.em.perms.2.',1))
#perms.1 <- get(paste0('bin.em.perms.1.',1))

for (i in 2:100){
 arraynum <- i
 load(file.path(mpath,paste0(pop,arraynum,'_scan_perms_bin_em.rsave')))
 perms.2 <- c(perms.2,get(paste0('bin.em.perms.2.',i)))
 #perms_1 <- c(perms,get(paste0('bin.em.perms.2.',i)))
}

pens <- calc.penalties(perms.2, alpha=0.05)

save.image(file.path(mpath,paste0(pop,'_all_perms_bin_em.rsave')))










################################################################################
load(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
################################################################################

summary(bin.em.2,perms=perms,alphas=0.1, pvalues=T)
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
