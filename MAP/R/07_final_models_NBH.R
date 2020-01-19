#!/bin/R

library('qtl')
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'.mapped.tsp.csv')
fl <- file.path(mpath,fl)


##pens
##    main    heavy    light
##3.971724 6.831380 3.617315

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

############################################################
cross <- jittermap(cross, amount=1e-6)

plot_test('map.png',width=1500)
plot.map(pull.map(cross),est.map(cross))
dev.off()

############################################################
dups <- findDupMarkers(cross, exact.only=F, adjacent.only=F)

cross2 <- pull.markers(cross, names(dups))

gg_marks <- unlist(lapply(1:24,function(X) { pickMarkerSubset(pull.map(cross)[[X]], 0.75)} ))
gg <- pull.markers(cross,gg_marks)



if(pop == 'ELR.missing') gg_marks <- c(gg_marks,"AHR2a_del")
#

cross2 <- sim.geno(cross2, stepwidth="fixed", step=1,off.end=5, error.prob=erp ,map.function="kosambi", n.draws=1000)
na <- stepwiseqtl(cross2, incl.markers=T, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=8)



out.aqi <- addqtl(hyper, qtl=rqtl, formula=y~Q1+Q2+Q3*Q4+Q4*Q5)
out.aq <- addqtl(hyper, qtl=rqtl, formula=y~Q1+Q2+Q3*Q4)

Figure 4: LOD curves for adding one QTL, interacting with the chromosome 15 locus, to the 4-QTL model, with the hyper data.
out.ap <- addpair(hyper, qtl=rqtl, chr=1, formula=y~Q2+Q3*Q4, verbose=FALSE)

addtoqtl, dropfromqtl, and replaceqtl
