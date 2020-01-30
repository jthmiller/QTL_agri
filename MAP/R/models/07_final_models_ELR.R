#!/bin/R
pop <- 'ELR'
library('qtl')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)

load(file.path(mpath,paste0(pop,1,'_scan_perms_bin_em.rsave')))
perms <- get(paste0('bin.em.perms.2.',1))

for (i in 2:100){
 arraynum <- i
 load(file.path(mpath,paste0(pop,arraynum,'_scan_perms_bin_em.rsave')))
 perms <- c(perms,get(paste0('bin.em.perms.2.',i)))
 perms_1 <- c(perms,get(paste0('bin.em.perms.2.',i)))
}

pens <- calc.penalties(perms, alpha=0.05)

save.image(file.path(mpath,paste0(pop,'_all_perms_bin_em.rsave')))


################################################################################
load(file.path(mpath,paste0(pop,'_scan2_bin_em.rsave')))
################################################################################


summary(bin.em.2, thresholds=c(0, Inf, 6, Inf, Inf), what="int")

summary(bin.em.2, thresholds=c(0, Inf, 6, Inf, Inf), what="int")
summary(bin.em.2, thresholds=c(9, Inf, 0, Inf, Inf), what="full")


summary(bin.em.2,perms=perms,alphas=0.05, pvalues=T)
summary(bin.em.2,perms=perms,alphas=0.15, pvalues=T, what='int')

summary(bin.em.2,perms=perms,alphas=0.15)













### OLD
############################################################
load(file.path(mpath,paste0(pop,'_norm_imp.rsave')))
############################################################

############################################################
summary(full.norm.imp)

      name chr pos n.gen
Q1  2@92.0   2  92     3
Q2 13@-2.0  13  -2     3
Q3 13@46.0  13  46     3
Q4 18@19.0  18  19     3
Q5 18@24.0  18  24     3
Q6 18@70.0  18  70     3
Q7 23@66.0  23  66     3
Q8 24@69.0  24  69     3
Q9 24@86.0  24  86     3

  Formula: y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9

  pLOD:  3.909

################################################################################
fitqtl(gg_step2, pheno.col=5, qtl=full.norm.imp, method="imp",model="normal",get.ests=F,covar=data.frame(gg_step2$pheno$sex),
 formula = "y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9")
################################################################################

      df        SS         MS      LOD     %var Pvalue(Chi2) Pvalue(F)
Model 18 328.51974 18.2510967 33.17583 85.89615            0         0
Error 59  53.94179  0.9142676       NA       NA           NA        NA
Total 77 382.46153         NA       NA       NA           NA        NA

$result.drop
        df Type III SS       LOD      %var  F value Pvalue(Chi2)    Pvalue(F)
2@92.0   2    36.71084  8.792747  9.598571 20.07664 1.611582e-09 2.233961e-07
13@-2.0  2    75.83076 14.869041 19.827028 41.47077 1.332268e-15 5.661138e-12
13@46.0  2    40.86174  9.551068 10.683883 22.34671 2.811462e-10 5.963132e-08
18@19.0  2    33.76198  8.232624  8.827550 18.46395 5.852968e-09 5.925962e-07
18@24.0  2    22.17884  5.833504  5.798973 12.12929 1.467224e-06 3.867925e-05
18@70.0  2    41.81793  9.721043 10.933892 22.86963 1.900891e-10 4.435101e-08
23@66.0  2    74.56873 14.703520 19.497053 40.78058 1.998401e-15 7.552847e-12
24@69.0  2    65.04305 13.399084 17.006428 35.57112 3.985701e-14 7.324785e-11
24@86.0  2    56.75691 12.176470 14.839901 31.03955 6.661338e-13 6.160177e-10

################################################################################
################################################################################

summary(full.norm.hk)
      name chr pos n.gen
Q1  5@89.0   5  89     3
Q2 12@25.0  12  25     3
Q3 13@10.0  13  10     3
Q4 14@54.0  14  54     3
Q5 15@27.0  15  27     3
Q6 18@53.0  18  53     3
Q7 23@76.0  23  76     3
Q8 24@55.0  24  55     3

  Formula: y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q5:Q8 + Q2:Q4 + Q1:Q7

  pLOD:  9.949
################################################################################
fitqtl(gg_step2, pheno.col=5, qtl=full.norm.hk, method="hk",model="normal",get.ests=T,covar=data.frame(gg_step2$pheno$sex),
 formula = "y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q5:Q8 + Q2:Q4 + Q1:Q7")
###
$result.full
      df        SS         MS      LOD     %var Pvalue(Chi2) Pvalue(F)
Model 28 357.42935 12.7653339 46.17947 93.45498            0         0
Error 49  25.03218  0.5108608       NA       NA           NA        NA
Total 77 382.46153         NA       NA       NA           NA        NA

$result.drop
                df Type III SS       LOD      %var   F value Pvalue(Chi2)
5@89.0           6    27.73328 12.630191  7.251260  9.047891 1.061389e-10
12@25.0          6    46.94578 17.889332 12.274640 15.315905 1.110223e-15
13@10.0          2    85.54899 25.162117 22.367998 83.730230 0.000000e+00
14@54.0          6    42.90348 16.910366 11.217725 13.997120 9.769963e-15
15@27.0          6   115.51479 29.223589 30.202982 37.686320 0.000000e+00
18@53.0          2    75.25607 23.507303 19.676767 73.656133 0.000000e+00
23@76.0          6    19.01314  9.570641  4.971255  6.202973 7.145068e-08
24@55.0          6   111.91893 28.784608 29.262794 36.513183 0.000000e+00
5@89.0:23@76.0   4    13.18779  7.167875  3.448135  6.453710 1.189264e-06
12@25.0:14@54.0  4    40.66522 16.342931 10.632499 19.900342 1.776357e-15
15@27.0:24@55.0  4    99.20406 27.134235 25.938310 48.547497 0.000000e+00

hknorm <- scanone(cross, method = "imp", model = "normal", pheno.col = 5, intcovar= data.frame(gg_step2$pheno$sex))
hknorm <- scanone(gg_step2, method = "hk", model = "normal", pheno.col = 5, intcovar= data.frame(gg_step2$pheno$sex))

plot_test('elr_full.norm.imp.png',width=1500)
plotLodProfile(full.norm.imp,incl.markers=F, main = 'ELR QTL')
dev.off()



############################################################
############################################################
## manual add stepwise qtl
bin.add.em.perms <- scanone(gg_step2, pheno.col=4, model='binary', method = "hk", n.perm = 10000, n.cluster=10)
lod <- summary(bin.add.em.perms)[2]
bin.add.em <- scanone(gg_step2, pheno.col=4, model='binary', method = "hk")
qtl <- summary(bin.add.em,lod)
bin.add.em.qtls1 <- makeqtl(gg_step2, chr=qtl[['chr']], pos=qtl[['pos']], what="prob")
bin.add.em.qtls1 <- refineqtl(gg_step2, qtl=bin.add.em.qtls1, pheno.col=4, model='binary', method = "hk", incl.markers=F)
int.em <- addint(gg_step2, qtl=bin.add.em.qtls1, formula=y~Q1+Q2, method='hk')
bin.add.em.qtls2 <- refineqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls1, method='hk', incl.markers=F, formula=y~Q1+Q2+Q1:Q2)
#int.em <- addint(gg_step2, qtl=bin.add.em.qtls, formula=y~Q1+Q2+Q3+Q1:Q3, method='hk')
##scan for an additional additive QTL
add.em <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q1:Q2)
qtl <- max(add.em)
bin.add.em.qtls3 <- addtoqtl(gg_step2, qtl=bin.add.em.qtls2,chr=qtl$chr, pos=qtl$pos)
bin.add.em.qtls3 <- refineqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls3, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q2)
#####
add.em_a <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls3, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q2)
## no additional additive.
## scan for interaction
add.em_i1 <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls3, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q2+Q1:Q4)
add.em_i2 <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls3, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q2+Q2:Q4)
add.em_i3 <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls3, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q2+Q3:Q4)

qtl <- max(add.em_i1)
bin.add.em.qtls4 <- addtoqtl(gg_step2, qtl=bin.add.em.qtls3,chr=qtl$chr, pos=qtl$pos)
bin.add.em.qtls4 <- refineqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls4, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q4+Q1:Q2+Q1:Q3)
int.em <- addint(gg_step2, qtl=bin.add.em.qtls4, formula=y~Q1+Q2+Q3+Q4+Q1:Q2+Q1:Q3, method='hk')
##no  interactions
add.em_a <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls4, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q4+Q1:Q2+Q1:Q3)


save.image(file.path(mpath,paste0(pop,'_models.rsave')))
### no additional


qtl <- max(add.em_a)
bin.add.em.qtls <- addtoqtl(gg_step2, qtl=bin.add.em.qtls,chr=qtl$chr, pos=qtl$pos)
bin.add.em.qtls <- refineqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q2+Q1:Q3)






out.fq <- fitqtl(gg_step2, pheno.col=4, method='hk', qtl=bin.add.em.qtls, model='binary', formula=y~Q1+Q2+Q3+Q1:Q2+Q1:Q3)



##scan for an additional interactive QTL
add.em.a <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q1:Q4)
add.em.b <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q2:Q4)
add.em.c <- addqtl(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3+Q3:Q4)

## Scan for interacting pair to add (long)
add.em.c <- addpair(gg_step2, pheno.col=4, model='binary', qtl=bin.add.em.qtls, method='hk', incl.markers=F, formula=y~Q1+Q2+Q3+Q1:Q3)
############################################################
############################################################
