



### CHR2 most often has distorted 2 locus interactions


### 2 and 24 have large 2 locus distortion
### Fall over HSP and AIP
geno.crosstab(cross,'24:6637513','2:27956730')
geno.crosstab(cross,'22:18254439','2:33399485')

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


geno.crosstab(cross,'24:10652460','2:27956730')
 geno.crosstab(cross,'6:15097723', '18:19515505')

geno.crosstab(cross,'24:809918','2:34936971')


plot_test('sdf')
plotPXG(cross,pheno.col=4,c('2:34936971','24:809918'))
dev.off()


AIP x HSP
geno.crosstab(cross,'22:19409957','2:28570451')


AHRa seems to have higher lod with AHRb
geno.crosstab(cross,'1:2698959','18:20565637')

geno.crosstab(cross,'1:857165','18:20565637')

## ARE THE AHRS LINKED?
pull.rf(cross,what='lod')['1:857165','18:20565637']


geno.crosstab(cross,'13:22726743','2:34484101')

geno.crosstab(cross,'24:30295052','3:29581625')
