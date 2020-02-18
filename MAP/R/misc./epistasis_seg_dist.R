
## 2 x 8 #####
               chr      pos chr      pos   lod_ab    lod_p       rfs
2:28513083.644   2 103.2171   8 32.72228 2.748198 20.03069 0.3470829

geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,8,32), mname2 = find.marker(rf, 2, 103))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,8,32), mname2 = find.marker(rf, 2, 103))
geno.crosstab(rf,mname1 = find.marker(rf,8,32), mname2 = find.marker(rf, 2, 103))

rf_p[1,] 2 103.21706   8 32.72228
ab_p[1,] 2  99.9107   8 36.05828
###############
## 2 x 18 #####
                chr       pos chr      pos   lod_ab    lod_p       rfs
2:37037255.4400   2  89.96982  18 29.90474 1.659472 21.22508 0.3601613





2:19875370.792     2 68.93531  14 44.45958 3.215742  9.9434225 0.3250112

## highest lod linkage
## 7 and 19
7:8827636.880      7 35.96795  19 95.17656 3.665627  1.1024195 0.2919717
geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,7,35), mname2 = find.marker(rf, 19, 95))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,7,35), mname2 = find.marker(rf, 19, 95))
geno.crosstab(rf,mname1 = find.marker(rf,7,35), mname2 = find.marker(rf, 19, 95))

## 14 and 15 and 19
14:32236920.23    14 81.66256  15 21.99959 3.225330  0.2036094 0.3178513
15:7960304.505    15 39.70502  19 63.14575 3.126290  0.2749022 0.6845876

geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,14,81), mname2 = find.marker(rf, 15, 22))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,14,81), mname2 = find.marker(rf, 15, 22))
geno.crosstab(rf,mname1 = find.marker(rf,14,81), mname2 = find.marker(rf, 15, 22))

## 3 and 17
3:33822959.959     3 76.94095  17 26.19020 3.379588  1.6330361 0.6959731
geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,3,81), mname2 = find.marker(rf, 17, 26))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,3,81), mname2 = find.marker(rf, 17, 26))
geno.crosstab(rf,mname1 = find.marker(rf,3,81), mname2 = find.marker(rf, 17, 26))

## 1 and 18
1:1088116.1255     1   2.725096  18 55.639611 2.702734  7.3971236 0.6639803
geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,1,2), mname2 = find.marker(rf, 18, 55))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,1,2), mname2 = find.marker(rf, 18, 55))
geno.crosstab(rf,mname1 = find.marker(rf,1,2), mname2 = find.marker(rf, 18, 55))

##
1:191503.1326      1  2.174072  10 40.335111 1.887185 2.2997859 0.3443818

################################################################################

find.marker(rf,18,48)

#### incompatability ###########################################################
geno.crosstab(rf,mname1 = "18:18316845", mname2 = "10:27910103")
################################################################################

mrks <- c("1:191503",'1:1088116','2:28513083','8:1621500',"8:35172047","10:24350254",'18:18590405','18:17539111')
mrks.gts <- gts[,mrks]

trs <- c('AA','AB','BB')

table(paste0(trs[mrks.gts[,1]],trs[mrks.gts[,7]]))
table(paste0(trs[mrks.gts[,1]],trs[mrks.gts[,8]]))
table(paste0(trs[mrks.gts[,2]],trs[mrks.gts[,7]]))
table(paste0(trs[mrks.gts[,2]],trs[mrks.gts[,8]]))

table(paste0(trs[mrks.gts[,1]],trs[mrks.gts[,6]],trs[mrks.gts[,7]]))
table(paste0(trs[mrks.gts[,1]],trs[mrks.gts[,6]],trs[mrks.gts[,8]]))

table(paste0(trs[mrks.gts[,1]],trs[mrks.gts[,6]],trs[mrks.gts[,8]]))

##cbind(trs[mrks.gts[,1]],trs[mrks.gts[,2]],trs[mrks.gts[,3]],trs[mrks.gts[,4]])

table(paste0(trs[mrks.gts[,1]],trs[mrks.gts[,2]],trs[mrks.gts[,3]],trs[mrks.gts[,4]]))
table(paste0(trs[mrks.gts[,1]],trs[mrks.gts[,4]]))
table(paste0(trs[mrks.gts[,1]],trs[mrks.gts[,3]]))

table(paste0(trs[mrks.gts[,2]],trs[mrks.gts[,3]],trs[mrks.gts[,4]]))

table(paste0(trs[mrks.gts[,2]],trs[mrks.gts[,4]]))
table(paste0(trs[mrks.gts[,3]],trs[mrks.gts[,4]]))
table(paste0(trs[mrks.gts[,1]],trs[mrks.gts[,2]]))
table(paste0(trs[mrks.gts[,4]],trs[mrks.gts[,5]]))
table(paste0(trs[mrks.gts[,5]],trs[mrks.gts[,6]]))

table(paste0(trs[mrks.gts[,4]],trs[mrks.gts[,5]],trs[mrks.gts[,6]]))

table(paste0(trs[mrks.gts[,5]],trs[mrks.gts[,7]]))


tail(links[order(links$lod_ab),],100)



plot_test('18_10')
effectplot(rf, pheno.col=4, mname1 = '18:17539111', mname2 = '10:24350254')
dev.off()

geno.crosstab(rf,mname1 = "18:18316845", mname2 = "10:27910103")

"1:191503",'1:1088116'
plot_test('18_10')
effectplot(rf, pheno.col=4, mname1 = "18:18316845", mname2 = "10:27910103")
dev.off()


"1:191503",'1:1088116'
plot_test('18_10')
effectplot(rf, pheno.col=4, mname1 = "18:18316845", mname2 = "10:27910103")
dev.off()


mrks <- c("1:191503","10:24350254","18:18316845","8:35172047")
mrks <- c("1:191503",'1:1088116','2:28513083','8:1621500',"8:35172047","10:24350254",'18:18590405','18:17539111')
mrks.gts <- gts[,mrks]
trs <- c('AA','AB','BB')

two_loc <- names(table(paste0(trs[mrks.gts[,1]],trs[mrks.gts[,2]])))

table(factor(paste0(trs[mrks.gts[,2]],trs[mrks.gts[,4]]), two_loc))
table(factor(paste0(trs[mrks.gts[,1]],trs[mrks.gts[,4]]), two_loc))
table(factor(paste0(trs[mrks.gts[,1]],trs[mrks.gts[,2]]), two_loc))
table(factor(paste0(trs[mrks.gts[,1]],trs[mrks.gts[,3]]), two_loc))
table(factor(paste0(trs[mrks.gts[,2]],trs[mrks.gts[,3]]), two_loc))
table(factor(paste0(trs[mrks.gts[,3]],trs[mrks.gts[,4]]), two_loc))

thr_loc <- unique(paste0(trs[mrks.gts[,1]],trs[mrks.gts[,2]],trs[mrks.gts[,3]]))
sort(table(factor(paste0(trs[mrks.gts[,1]],trs[mrks.gts[,2]],trs[mrks.gts[,3]]), thr_loc)))

sort(table(factor(paste0(trs[mrks.gts[,1]],trs[mrks.gts[,2]],trs[mrks.gts[,3]]), thr_loc)))




geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,2,11.9), mname2 = find.marker(rf, 18, 50))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,2,11.9), mname2 = find.marker(rf, 18, 50))
geno.crosstab(rf,mname1 = find.marker(rf,2,11.9), mname2 = find.marker(rf, 18, 50))


geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,2,100), mname2 = find.marker(rf, 18, 54))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,2,100), mname2 = find.marker(rf, 18, 54))
geno.crosstab(rf,mname1 = find.marker(rf,2,100), mname2 = find.marker(rf, 18, 54))


geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,2,42.9), mname2 = find.marker(rf, 18, 50.5))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,2,42.9), mname2 = find.marker(rf, 18, 50.5))
geno.crosstab(rf,mname1 = find.marker(rf,2,42.9), mname2 = find.marker(rf, 18, 50.5))


plot_test('18_elr')
effectplot(rf, pheno.col=4, mname1 = find.marker(rf, 18, 50.5) , ylim=c(0,1))
dev.off()

plot_test('13_elr')
effectplot(rf, pheno.col=4, mname1 = find.marker(rf, 13, 25.9) , ylim=c(0,1))
dev.off()


## ELR INCOMPAT
geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,13,74), mname2 = find.marker(rf, 2, 78))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,13,74), mname2 = find.marker(rf, 2, 78))
geno.crosstab(rf,mname1 = find.marker(rf,13,74), mname2 = find.marker(rf, 2, 78))


2:13279585.717    2 54.04262    13  6.695413 1.91791332  5.233144 0.3508211

p_tables[[12]]
ab_tables[[12]]
rf_tables[[12]]

### NBH incompat
geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,13,48), mname2 = find.marker(rf, 2, 89))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,13,48), mname2 = find.marker(rf, 2, 89))
geno.crosstab(rf,mname1 = find.marker(rf,13,48), mname2 = find.marker(rf, 2, 89))

plot_test('2_13_nbh')
effectplot(rf, pheno.col=4, mname1 = find.marker(rf, 13, 48) ,mname2 = find.marker(rf, 2, 89), ylim=c(0,1))
dev.off()



## missing gts in elr
plot_test('2_13_elr')
effectplot(rf, pheno.col=4, mname1 = find.marker(rf, 13, 25) ,mname2 = find.marker(rf, 2, 77.8), ylim=c(0,1))
dev.off()



plot_test('2_18_nbh')
effectplot(rf, pheno.col=4, mname1 = find.marker(rf, 18, 47) ,mname2 = find.marker(rf, 2, 89), ylim=c(0,1))
dev.off()

plot_test('13_18_nbh')
effectplot(rf, pheno.col=4, mname2 = find.marker(rf, 18, 47) ,mname1 = find.marker(rf, 13, 43), ylim=c(0,1))
dev.off()



plot_test('1_13_nbh')
effectplot(rf, pheno.col=4, mname2 = find.marker(rf, 13, 0) ,mname1 = find.marker(rf, 1, 0), ylim=c(0,1))
dev.off()


######################################
## missing gts in elr
plot_test('2_13_elr')
effectplot(rf, pheno.col=4, mname1 = find.marker(rf, 13, 63) ,mname2 = find.marker(rf, 2, 77.8), ylim=c(0,1))
dev.off()

plot_test('2_13_elr')
effectplot(rf, pheno.col=4, mname1 = find.marker(rf, 13, 74) ,mname2 = find.marker(rf, 2, 77), ylim=c(0,1))
dev.off()

plot_test('2_13_elr')
effectplot(rf, pheno.col=4, mname1 = find.marker(rf, 13, 50) ,mname2 = find.marker(rf, 2, 100), ylim=c(0,1))
dev.off()



plot_test('18_2_elr')
effectplot(rf, pheno.col=4, mname2 = "2:35889898", mname1 = find.marker(rf, 18, 0), ylim=c(0,1))
dev.off()

plot_test('18_2_elr')
effectplot(rf, pheno.col=4, mname2 = "2:31789055", mname1 = find.marker(rf, 18, 0), ylim=c(0,1))
dev.off()

2:2018841.5190     2 11.986618    18 50.532743 0.081384390  5.495940 0.4688252

plot_test('18_2_elr')
effectplot(rf, pheno.col=4, mname2 = '2:35889898', mname1 = find.marker(rf, 18, 0), ylim=c(0,1))
dev.off()


plot_test('18_10_elr')
effectplot(rf, pheno.col=4, mname1 = find.marker(rf, 10, 40), mname2 = find.marker(rf, 18, 48), ylim=c(0,1))
dev.off()

plot_test('2_10_elr')
effectplot(rf, pheno.col=4, mname1 = find.marker(rf, 10, 40), mname2 = find.marker(rf, 2, 103), ylim=c(0,1))
dev.off()

plot_test('1_18_elr')
effectplot(rf, pheno.col=4, mname1 = find.marker(rf, 1, 0), mname2 = find.marker(rf, 18, 100), ylim=c(0,1))
dev.off()

######################################

geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,1, 10), mname2 = find.marker(rf,  24, 57))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,1, 10), mname2 = find.marker(rf,  24, 57))
geno.crosstab(rf,mname1 = find.marker(rf,1, 10), mname2 = find.marker(rf,  24, 57))


p_tables[[1]]
ab_tables[[1]]
rf_tables[[1]]


plot_test('2_19_nbh')
effectplot(rf, pheno.col=4, mname1 = find.marker(rf, 19, 18) ,mname2 = find.marker(rf, 2, 103), ylim=c(0,1))
dev.off()

summary(bin.em.2, thresholds=c(0, Inf, 8, Inf, Inf), what="int")

2:15573468.548     2 56.795857    10 28.2079147 0.006681401 4.5348628 0.5076915
10:24350254.519   10  40.335111    18  48.191587 2.8393153 7.87401639 0.6932931


### there are no chr18 AA x AA chr10
plot_test('18_10_elr')
effectplot(rf, pheno.col=4, mname1 = find.marker(rf, 10, 40), mname2 = find.marker(rf, 18, 48), ylim=c(0,1))
dev.off()

geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = find.marker(rf,18, 60), mname2 = find.marker(rf,  2, 80))
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = find.marker(rf,18, 60), mname2 = find.marker(rf,  2,80))
geno.crosstab(rf,mname1 = find.marker(rf,18, 60), mname2 = find.marker(rf,  2,80))



np <- scanone(rf,pheno.col=1,method="em",model="normal",addcovar=g)


np <- scanone(rf,pheno.col=4,method="em",model="bin",intcovar=g)


intcovar=
np <- scanone(rf,pheno.col=1,method="em",model="normal",addcovar=g)


s1n <- scanone(rf, pheno.col=1, model="normal", method="imp", addcovar=g, intcovar=int)


s1ig <- scanone(rf, pheno.col=1, model="normal", method="imp", addcovar=int, intcovar=g)


s1ig <- scanone(rf, pheno.col=4, model="binary", method="imp", addcovar=g, intcovar=int)
s1ig <- scanone(rf, pheno.col=4, model="binary", method="imp", addcovar=int, intcovar=int)


so <- summary(np)[cov,]
top_2 <- order(so$lod,decreasing =T)[1]

mar <- find.pseudomarker(rf, 24, 35.17745)
int <- lapply(mar,function(X){ pull.argmaxgeno(cross)[,X] } )
names(int) <- mar
int <- lapply(int, function(X,Y){ cbind(as.numeric(X==1), as.numeric(X==2))} )
int <- data.frame(do.call(cbind,int))

mar <- find.pseudomarker(rf, 2, 101)
g <- lapply(mar,function(X){ pull.argmaxgeno(cross)[,X] } )
names(g) <- mar
g <- lapply(g, function(X,Y){ cbind(as.numeric(X==1), as.numeric(X==2))} )
g <- data.frame(do.call(cbind,g))

### there are no chr18 AA x AA chr10
plot_test('1_nbh')
effectplot(rf, pheno.col=4, mname1 = find.marker(rf, 1, 0), ylim=c(0,1))
dev.off()



sb1 <- scanone(cross,pheno.col=4,method="em",model="bin")
sn1 <- scanone(rf,pheno.col=4,method="em",model="normal")

sb1 <- scanone(rf,pheno.col=4,method="em",model="bin", intcovar=g, addcovar=g)




'1:3269009','2:37728488','3:29581996','8:6576898','18:17816908','24:33090064'


qtl <- c('2:28513083','18:17816908','13:22712588','8:6576898','24:33090064','10:24350254','1:3269009')

qtl <- c('1:3269009','2:13279585','2:37728488','3:29581996','8:6576898','13:20756769','18:17816908','24:33090064')

cbind(cross$pheno,pull.geno(cross)[,qtl])


qtl <- c('1:3269009','2:37728488','2:13279585','3:29581996','8:6576898','13:12740805','18:17816908','24:33090064')


z <- pull.markers(cross,qtl)
z$pheno$aip <- pull.geno(z)[,'2:13279585']
z$pheno$ahr <- pull.geno(z)[,'18:17816908']
z$pheno$arnt <- pull.geno(z)[,'13:12740805']
za <- subset(z,ind = z$pheno$Pheno > 3 )
zb <- subset(z,ind = z$pheno$Pheno < 3)

#zb <- subset(zb,ind = zb$pheno$ID[order(zb$pheno$aip, zb$pheno$ahr)])
#za <- subset(za,ind = za$pheno$ID[order(za$pheno$aip, za$pheno$ahr)])


plot_test('qtl_genos')
par(mfrow = c(2,1))
 geno.image(za, chr=c(2,8,13,18,24), reorder=c(6), main="Genotype data", alternate.chrid=FALSE, col = c(NA,'orange','black','blue'))
 geno.image(zb, chr=c(2,8,13,18,24), reorder=c(6), main="Genotype data", alternate.chrid=FALSE, col = c(NA,'orange','black','blue'))
dev.off()


, col=NULL, ...)
