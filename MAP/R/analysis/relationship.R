



########## Remove the samples related by more than 80% of genotypes ############
wh <- which(cpgt > 0.8, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
mats <- cbind(rownames(wh),colnames(cpgt)[as.numeric(wh[,2])])
toss.relat <- unique(apply(mats,1,function(X){
 X[which.max(c(nmissing(cross.1)[X[1]],nmissing(cross.1)[X[2]]))]
}))
#################################################################################


################################################################################
#### TEST SAMPLE GT SIMILARITY ##################################################
#gt_nopar <- geno.table(subset(cross,ind=!cross$pheno$ID %in% c('NBH_NBH1M','NBH_NBH1F')))
#parABxAB <- intersect(rownames(gt_nopar[which(gt_nopar$P.value > 0.0001),]) ,parABxAB)
#cross.1 <- pull.markers(cross,parABxAB)
##cross.1 <- subset(cross.1,ind=!cross.1$pheno$ID%in%c('NBH_NBH1M','NBH_NBH1F'))
#cpgt <- comparegeno(cross.1)
#colnames(cpgt) <- cross.1$pheno$ID
#rownames(cpgt) <- cross.1$pheno$ID
#cpgt[cpgt==NaN] <- NA
#diag(cpgt) <- NA
## remove completely empty
#cpgt <- cpgt[rowSums(is.na(cpgt)) < nind(cross.1),colSums(is.na(cpgt)) < nind(cross.1)  ]
#################################################################################
