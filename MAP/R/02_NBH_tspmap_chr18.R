#!/bin/R

i <- 18
pop <- 'NBH'
source("/home/jmiller1/QTL_agri/MAP/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- file.path(paste0(pop,'_unmapped_filtered.csv'))
##fl <- file.path(paste0(pop,'_unmapped_reassigned_markers')
##mapfile <- paste0(pop,'_all_mark_',i,'_tsp')

filename <- file.path(mpath,mapfile)

libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap")
suppressMessages(sapply(libs2load, require, character.only = TRUE))

################################################################################
cross <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
sex <- read.table(file.path(mpath,'sex.txt'),stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross$pheno$ID), 'sex']
cross$pheno$sex <- sex.vec
################################################################################
toss.missing <- c("NBH_5525","NBH_6177","NBH_5528")
gt <- geno.table(subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F')))
bfixA <- rownames(gt[which(gt$P.value > 0.0001 & gt$missing < 1),])
################################################################################
## 30k markers before hand.
###### FILTER #######################################################
cross <- pull.markers(cross,bfixA)
cross <- subset(cross,ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))
################################################################################

fl <- file.path(mpath,'NBH_pre_reform_markers')
write.cross(cross,filestem=fl,format="csv")
################################################################################
fl <- file.path(paste0(pop,'_pre_reform_markers.csv'))
cross <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)



reorg <- subset(cross2,chr=c(18,i))
reorg <- formLinkageGroups(cross, max.rf = 0.1, min.lod = 10, reorgMarkers = TRUE)



################################################################################
reorg <- formLinkageGroups(cross, max.rf = 0.1, min.lod = 10, reorgMarkers = TRUE)
fl <- file.path(mpath,'NBH_unmapped_reassigned_markers')
write.cross(reorg,filestem=fl,format="csv")
################################################################################








pull.map(cross2,

#cross <- subset(cross2,chr=i)

sex <- read.table(file.path(mpath,'sex.txt'),stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross$pheno$ID), 'sex']
cross$pheno$sex <- sex.vec

################################################################################

ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross)[[as.character(i)]]))))
cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
 maxit = 1, tol = 0.1, sex.sp = F)

#### CHR8 #######



cross <- subset(cross,ind=!cross$pheno$ID %in% c('NBH_NBH1M','NBH_NBH1F'))


png(paste0('~/public_html/NBH_gts_preclean',i,'.png'),height=2500,width=4500)
 plotGeno(cross, chr=i, cex=2)
dev.off()

################################################################################
cross <- cleanGeno_jm(cross, chr=i, maxdist=100, maxmark=8, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.05)
cross <- removeDoubleXO(cross, chr=i)
cross <- calc.errorlod(cross, err=0.05)
cross <- cleanGeno_jm_2(cross, chr=i, maxdist=50, maxmark=4, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.05)
################################################################################

################################################################################
gt <- geno.table(cross)
bfixA <- rownames(gt[which(gt$P.value > 0.0001 & gt$missing < 6),])
cross <- pull.markers(cross,bfixA)
################################################################################

################################################################################
cross <-tspOrder(cross = cross,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

cross <- shiftmap(cross, offset=0)

cross_map <-  est.map(cross, error.prob=0.04,map.function="kosambi",maxit=100000,tol=1e-7, sex.sp=FALSE, verbose=FALSE)

cross <- qtl:::replace.map(cross,cross_map)

write.cross(cross,chr=i,filestem=filename,format="csv")

png(paste0('~/public_html/NBH_RF_concord',i,'_tsp.png'))
  plotRF(cross)
dev.off()
