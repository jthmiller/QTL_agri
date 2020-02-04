high <- c(10,2,1,18)
med <- c(3,8,24,19)
low <- c(13,11,12,17,15,23)

no_qtl_mr <- scanone(rf, pheno.col=4, method="mr", model="binary")

nbh_gens

ahra <- find.marker(rf,1,0)
aip <- find.marker(rf,2,82)
arnt1 <- find.marker(rf,10,39)
arnt2 <- find.marker(rf,8,62)
ahrb <- find.marker(rf,18,67)

chr10 <- find.marker(rf,10,125)


ahra <- rownames(summary(no_qtl_mr))[1]
aip <- rownames(summary(no_qtl_mr))[2]
arnt1 <- rownames(summary(no_qtl_mr))[9]
arnt2 <- rownames(summary(no_qtl_mr))[7]
ahrb <- rownames(summary(no_qtl_mr))[17]

h <- quantile(rf.df, 0.999,na.rm=T)
l <- quantile(rf.df, 0.001,na.rm=T)
h9 <- quantile(rf.df, 0.95,na.rm=T)
l9 <- quantile(rf.df, 0.05,na.rm=T)

plot_test('nbh_1_2_10_18', width=1250,height=750)
par(mfrow = c(6,1))

plot(pull.rf(rf), ahra, ylim=c(0.3,0.75))
abline(h=h, col='red')
abline(h=h9, col='grey')
abline(h=0.5, col='black')
abline(h=l9, col='grey')
abline(h=l, col='red')

plot(pull.rf(rf), aip , ylim=c(0.3,0.75))
abline(h=h, col='red')
abline(h=h9, col='grey')
abline(h=0.5, col='black')
abline(h=l9, col='grey')
abline(h=l, col='red')

plot(pull.rf(rf), arnt1 , ylim=c(0.3,0.75))
abline(h=h, col='red')
abline(h=h9, col='grey')
abline(h=0.5, col='black')
abline(h=l9, col='grey')
abline(h=l, col='red')

plot(pull.rf(rf),arnt2 , ylim=c(0.3,0.75))
abline(h=h, col='red')
abline(h=h9, col='grey')
abline(h=0.5, col='black')
abline(h=l9, col='grey')
abline(h=l, col='red')


plot(pull.rf(rf), ahrb, ylim=c(0.3,0.75))
abline(h=h, col='red')
abline(h=h9, col='grey')
abline(h=0.5, col='black')
abline(h=l9, col='grey')
abline(h=l, col='red')

plot(pull.rf(rf), "7:8827636", ylim=c(0.3,0.75))
abline(h=h, col='red')
abline(h=h9, col='grey')
abline(h=0.5, col='black')
abline(h=l9, col='grey')
abline(h=l, col='red')

dev.off()

########################################################################################
summary(no_qtl_mr)
########################################################################################

getrf <- function(X, chr, mat){
 mat.order <- order(mat[X,])
 rfs <- mat[X,][mat.order]
 rfs[!is.na(rfs)]
}

chr1_rf <- getrf(X = ahra, chr = 1, mat = pull.rf(rf))

chr18_rf <- getrf(X = ahrb, chr = 18, mat = lod.df)

chr10_rf <- getrf(X = arnt1, chr = 10, mat = rf.df)

tail(sort(lod.df[,'18:17539111']))


geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = ahrb, mname2 = names(chr18_rf)[1363])
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = ahrb, mname2 = names(chr18_rf)[1363])
geno.crosstab(rf,mname1 = ahrb, mname2 = names(chr18_rf)[1363])


geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = arnt1, mname2 = ahrb)
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = arnt1, mname2 = ahrb)
geno.crosstab(rf,mname1 = arnt1, mname2 = ahrb)


geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = ahra, mname2 = ahrb)
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = ahra, mname2 = ahrb)
geno.crosstab(rf,mname1 = ahra, mname2 = ahrb)

geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = aip, mname2 = ahrb)
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = aip, mname2 = ahrb)
geno.crosstab(rf,mname1 = aip, mname2 = ahrb)

geno.table(rf,chr=1)

geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = '18:15890517', mname2 = '1:4188554')
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = '18:15890517', mname2 = '1:4188554')
geno.crosstab(rf,mname1 = '18:15890517', mname2 = '1:4188554')

geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = '7:8827636', mname2 = "19:42466531")
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = '7:8827636', mname2 = "19:42466531")
geno.crosstab(rf,mname1 = '7:8827636', mname2 = "19:42466531")

geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = '18:17539111', mname2 = '10:24350254')
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = '18:17539111', mname2 = '10:24350254')
geno.crosstab(rf,mname1 = '18:17539111', mname2 = '10:24350254')

geno.crosstab(subset(rf,ind=rf$pheno$bin == 0),mname1 = '18:18590405', mname2 = '1:4188554')
geno.crosstab(subset(rf,ind=rf$pheno$bin == 1),mname1 = '18:18590405', mname2 = '1:4188554')
geno.crosstab(rf,mname1 = '18:18590405', mname2 = '1:4188554')
