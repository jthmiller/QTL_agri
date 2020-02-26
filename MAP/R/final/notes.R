


26326053 26382053

get_genes(6,26326053,5)
get_genes(5,mean(c(16000000,21000000)),30)

get_genes(5,16000000,50)
get_genes(5,mean(c(14817605,16605791)),20)


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

### NBH 2 and 13 have highest lod effect
c2 :c13 52.88 30.54    26.62    6.04    8.89  17.724 -2.8542
cX:c24  7.84 13.9     9.79    5.82    6.88   2.910  -1.065

b <- find.marker(cross,13,30.5)

a <- find.marker(cross,2,52.88)

geno.crosstab(cross,a,b)


plot_test('sdf')
effectplot(cross,pheno.col=4,mname2=a,mname1=b, ylim=c(0,1))
dev.off()

plot_test('sdf')
effectplot(cross,pheno.col=1,mname2=a,mname1='13:24456390')
dev.off()


################################################################################
### NOTES FROM FIGURES (SEGDIST)
### ELR
CHR1 - D/F/NSF peaks at both ends (at aHR2a). Lod tops out at 2 here.
CHR2 - D/F/NSF at AIP
CHR3 - Narrow div and fst outlier at 19 MB
CHR4 -
CHR5 - Wide d/f doutlier at 16 and 18 MB
CHR6 - Narrow d/f outlier centered over max lod (only lod 2 or so)
CHR7 - none
CHR8 - Wide d/f doutlier at 18 MB, but lod is minimized there. Lod peak ~38 MB
CHR9 - A N/S fst peak, but not much else
CHR10 - NSF peak at 4 MB
CHR11 - NONE
CHR12 - NONE
CHR13 - LOD maximized at sensetive background 1 MB
CHR14 - Wide D/F/NSF at 22-24 MB. No phen LOD
CHR15 - NSF indicatess possible deletion at 1-2 MB
CHR16 - Wide D/F at 13-15 MB. No lod
CHR17 - None
CHR18 - DIVERSITY DOESNT CHANGE AT AHRb position
CHR19 - NONE
CHR20 - FNS outlier (no div change) at 34 MB
CHR21 -
CHR22 - D/F/FNS outlier at 19 MB
CHR23 - FNS outlier at 3 MB, Deletion? at 24 MB
CHR24 -

CHR1
CHR2
CHR3
CHR4
CHR5
CHR6
CHR7
CHR8
CHR9
CHR10
CHR11
CHR12
CHR13
CHR14
CHR15
CHR16
CHR17
CHR18
CHR19
CHR20
CHR21
CHR22
CHR23
CHR24
