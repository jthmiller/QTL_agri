

#chr10 23926932 23951650 XM_012876030.1  LOC105935564    aryl hydrocarbon receptor nuclear translocator-like protein 1
#NW_012234311.1 1350282 1390498 XM_012867564.1   arntl   aryl hydrocarbon receptor nuclear translocator-like%2C transcript variant X1
#NW_012225110.1 22934 115575 XM_012856879.1      exception=annotated by transcript or proteomic data     aryl-hydrocarbon receptor nuclear translocator 2

## ON CHR 10
## NW_012234311.1:1152974 (arntl)
## NW_012234311.1:706431 (arnt2)


### CHR2 most often has distorted 2 locus interactions

## MANUAL MODEL
## 2:27373969, 55.63432
## 18:20723840, 53.1464
plot_test('dsf')
effectplot(cross,pheno.col=5,mname2='2:27373969',mname1='18:20723840')
dev.off()

plot_test('dsf')
effectplot(cross,pheno.col=4,mname2='2:36080762',mname1='18:20723840')
dev.off()

plot_test('dsf')
effectplot(cross,pheno.col=4,mname2='2:36080762',mname1='18:17874376')
dev.off()

plot_test('dsf')
effectplot(cross,pheno.col=4,mname2='2:27373969',mname1='18:17874376')
dev.off()

18:17874376

plot_test('dsf')
plotPXG(cross,pheno.col=4,c('2:34937180','24:6767161'))
dev.off()

           chr missing AA AB BB not.BB not.AA     P.value
2:34937180   2       0 38 37 16      0      0 0.001001163
24:6767161  24       1  7 63 20      0      0 0.000114175




disto <- lapply(c(1:4,6:24), function(i) {
 gt <- geno.table(cross,i)
 gt[which.max(-log10(gt$P.value)),]
})
disto <- do.call(rbind,disto)
disto$P.value <- -log10(disto$P.value)
disto <- disto[order(disto$P.value),]




plot_test('sdf')
plotPXG(cross,pheno.col=4,c('2:34936971','24:809918'))
dev.off()

AIP x HSP
geno.crosstab(cross,'24:10291747','2:27757494')
geno.crosstab(cross,'22:19528880','2:27757494')
geno.crosstab(cross,'7:31714010','2:27757494')

AHRa x AHRb (AHRa seems to have higher lod with AHRb)
geno.crosstab(cross,'18:20367780','1:811175')

AIP x AHRb
geno.crosstab(cross,'18:20367780','2:27757494')

AIP x ARNT
geno.crosstab(cross,'8:16768182','2:27757494')
geno.crosstab(cross,'8:9792760','2:27757494')

AIP x AHRa
geno.crosstab(cross,'1:811175','2:27757494')

ARNT x ARNT (high segdist)
geno.crosstab(cross,'13:24355608','8:16768182')

ARNT and 20 (maxdist 8:9131715 20:19045496)
geno.crosstab(cross,'20:19045496','8:16768182')
geno.crosstab(cross,'20:19045496','8:9131715')
-log10(csq.pval['20:19045496','8:9131715'])

HSP and AIP
geno.crosstab(cross,'2:33395185','22:19528880')

which.max(-log10(csq.pval['2:27757494',])) == 8:9792760
which.max(-log10(csq.pval['2:33882931',])) == 22:17595037
################################################################################
## ARE THE AHRS LINKED?
pull.rf(cross,what='lod')['1:857165','18:20565637']
geno.crosstab(cross,'13:22726743','2:34484101')
geno.crosstab(cross,'24:30295052','3:29581625')
geno.crosstab(cross,'22:17595037','2:33395185')
geno.crosstab(cross,'22:19528880','2:33395185')
geno.crosstab(cross,'13:23470876','2:27757494')
geno.crosstab(cross,'24:10652460','2:27956730')
geno.crosstab(cross,'6:15097723', '18:19515505')
geno.crosstab(cross,'24:809918','2:34936971')


plot_test('sdf')
effectplot(cross,pheno.col=1,mname2='2:27757494',mname1='7:31714010')
dev.off()

geno.crosstab(cross,'7:31714010','2:27757494')

geno.crosstab(cross,'17:29007925','5:3014891')
geno.crosstab(cross,'8:37635736','17:7480177')

14:28023254 19:16701029
22:10368582  14:1719508


get_genes(chr,pos,ngens=2)
get_genes(17,29007925,ngens=10)
get_genes(17,7480177,ngens=10)
get_genes(8,37635736,ngens=10)
get_genes(5,3014891,ngens=10)

### MOST DISTORTED MARKER ON 17
17:29007925 17:29388433 17:14629450

geno.crosstab(cross,'17:14629450','24:2123083')

## MANUAL MODEL
## 2:27373969, 55.63432
## 18:20723840, 53.1464
geno.crosstab(cross,'2:36080762','2:27373969')

### ELR 2 and 13 appear linked
