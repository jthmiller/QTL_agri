
###### ADD UNMAPPED TO TEST FOR LINKAGE
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- 'NBH_imputed_2unmapped_tsp.csv'
fl <- file.path(mpath,fl)

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
 estimate.map = FALSE
)

 toss.missing <- c("NBH_5525","NBH_6177","NBH_5528","NBH_6137","NBH_6125")
 cross <- subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))

 cross <- formLinkageGroups(cross, max.rf = 0.5, min.lod = 0.01, reorgMarkers = TRUE)
 names(cross$geno) <- 'NW'




 mpath <- '/home/jmiller1/QTL_agri/data'
 mapfile <- paste0(pop,'_imputed_2unmapped_tsp')
 filename <- file.path(mpath,mapfile)
 write.cross(cross,filestem=filename,format="csv")

rf <- est.rf(cross)


 pval <- 1.0e-5
 mis <- misg(rf,0.15)
 mis <- 10
 bfixA <- rownames(gt[which(gt$P.value > pval & gt$missing < mis),])
 cross2 <- pull.markers(rf,bfixA)

 cross1 <- formLinkageGroups(cross2, max.rf = 0.1, min.lod = 10, reorgMarkers = TRUE)



unmapped NW_012234311.1 are on chr10 and


gt <- geno.table(cross1)

chr <- unique(gt[grep('NW',rownames(gt)),'chr'])

markernames(cross1, chr=chr)

c2 <- subset(cross1, chr=1:29)
