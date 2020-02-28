

## AHR x CYP1b1 + (incompatable AIP)
## CYP1b1 uncoupled more important?
### ELR (18 (ahr) is interacting with 15 (cyp1b1) (effect on phenotype)
### ELR (2 (AIP) is incompatable with 13 (arnt) (interaction has no affect on phenotype, but may alter efffects of AIP)

### NBH AIP x ARNT X AHR
### NBH (2 (AIP) interacts with 13 (arnt) (effect on phenotype, but is this due to seg dist?)
### NBH (18 does NOT appear to interact with 15
###

### THE OUTLIERS











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


geno.crosstab(cross,'2:36080762','2:27373969')

geno.crosstab(cross,'2:29113987','13:22410721')


cross_s <- switchAlleles(cross_s, markernames(cross_s,2))
geno.crosstab(cross,'2:29113987','13:22410721')

cross_s <- switchAlleles(cross_s,markernames(cross_s,2))
plot_test('sdf')
effectplot(cross_s,pheno.col=5,mname2=qtl2,mname1='13:24113176', ylim=c(0,5))
dev.off()

sone <- scanone(cross_s, pheno.col=4, model="normal", method="hk")
summary(sone)

qtl2 <- '2:29113987'
rn13 <- '13:22410721'
rn2 <- '2:33578802'

aip <- '2:27502173'
arnt <- '13:24472329'
hsp <- '13:23536332'
tmtc2a <- '13:10592117'

13:21347307  13       2 20 43 15      0      0 0.4815384
13:22410721  13       0 20 45 15      0      0 0.3916056
13:21996131  13       3 20 42 15

markernames(cross,2)

which.max(-log10(csq.pval['13:22410721',markernames(cross,2)]))

geno.crosstab(cross,'2:35401205','13:22410721')

geno.crosstab(cross,'2:23727548','13:22410721')

sort(-log10(csq.pval['13:22410721',markernames(cross,2)]))

sort(-log10(csq.pval['2:28101458',markernames(cross,13)]))

geno.crosstab(cross,'2:28101458','13:1369425')



geno.crosstab(cross,'18:20492770','15:4202438')


pop <- 'NBH'
library('qtl')
library('snow')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
cores <- 6
load(file.path(mpath,paste0(pop,'_scan2_normal_em.rsave')))


geno.crosstab(cross,'2:33640721','13:14492248')

#### NBH INCOMPATABILITY AT 13 v 2
#### WHAT ELSE IS 13 incompatable with?
c2 :c13 80.59 39.42    27.32    4.68    9.38  17.937 -4.7046
################################################################################


##### HEATMAP ######################################################################
csq.pval.hm <- data.matrix(-log10(csq.pval[markernames(cross,2),markernames(cross,13)]))
plot_test('heatmap_dist_elr',height=1000,width=1000)
heatmap(csq.pval.hm)
dev.off()
######################################################################

ELR (popgen outliers, low/no LOD): 1,2
ELR (popgen outliers, med LOD): 5,8,
ELR (popgen outliers, hi LOD):

WIDE: 1,2,5,8,14,16,22

## AHR x CYP1b1 + (incompatable AIP)
## CYP1b1 uncoupled more important?
### ELR (18 (ahr) is interacting with 15 (cyp1b1) (effect on phenotype)
### ELR (2 (AIP) is incompatable with 13 (arnt) (interaction has no affect on phenotype, but may alter efffects of AIP)

### NBH AIP x ARNT X AHR
### NBH (2 (AIP) interacts with 13 (arnt) (effect on phenotype, but is this due to seg dist?)
### NBH (18 does NOT appear to interact with 15
###

################################################################################
### NOTES FROM FIGURES (SEGDIST)
### ELR
**CHR1 - D/F/NSF peaks at both ends (at aHR2a). Lod tops out at 2 here.
**CHR2 - D/F/NSF at AIP
CHR3 - Narrow div and fst outlier at 19 MB
CHR4 - NONE
**CHR5 - Wide d/f doutlier at 16 and 18 MB
CHR6 - Narrow d/f outlier centered over max lod (only lod 2 or so)
CHR7 - NONE
**CHR8 - Wide d/f doutlier at 18 MB, but lod is minimized there. Lod peak ~38 MB
CHR9 - A N/S fst peak, but not much else
CHR10 - NSF peak at 4 MB
CHR11 - NONE
CHR12 - NONE
**CHR13 - LOD maximized at sensetive background 1 MB
**CHR14 - Wide D/F/NSF at 22-24 MB. No phen LOD
CHR15 - NSF indicatess possible deletion at 1-2 MB
CHR16 - Wide D/F at 13-15 MB. No lod
CHR17 - None
**CHR18 - DIVERSITY DOESNT CHANGE AT AHRb position
CHR19 - NONE
CHR20 - FNS outlier (no div change) at 34 MB
CHR21 - NONE - Possible FNS
**CHR22 - Wide D/F/FNS outlier at 19 MB
CHR23 - FNS outlier at 3 MB, Deletion? at 24 MB
CHR24 - NONE

## NBH
**CHR1 - LOW LOD MAXIMIZED NEAR AHRa
**CHR2 - WIDE F/D/FNS
CHR3 - NONE
CHR4 - NONE
**CHR5 - THE WIDE OUTLIER PRESENT IN ELR IS NOT IN NBH
CHR6 -
CHR7 - SUSPICIOUS FST AND DELTA PI
CHR8 - LESS PRONOUNCED INTERVAL THAN ELR
CHR9
CHR10 - SUSPICOUS FST and DELTA pi
CHR11
CHR12
CHR13 - SOMEWHAT MAXIMIZED OVER ARNT
CHR14 - POSSIBLE FST/D OUTLIERS
CHR15
CHR16
CHR17
CHR18 - FS, But diversity does not tank as it does with other outliers
CHR19
CHR20
CHR21
CHR22 - some possible outliers
CHR23 - some good outliers
CHR24 - nothing compelling
