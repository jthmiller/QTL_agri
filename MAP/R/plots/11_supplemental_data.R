#!/bin/R

####################################################################################





################################################################################
## FUNCTIONS ###
################################################################################
################################################################################

plot_ef <- function(crs,map,pr,ahr,popgen,chs,main,model=c("bin","pheno_norm"),...){

 for (chr in chs){

  c2eff <- scan1coef(pr[,as.character(chr)], crs$pheno[,model])

  plot(c2eff, map[as.character(chr)], columns=1:3, col=col, ylim=c(0,1), cex.axis = 2,main=main,...)

    if(any( chr %in% ahr$chr )) {
      indx <- which(ahr$chr %in% chr)
      abline(v=as.numeric(ahr[indx,'pos1']), col='red',lwd=0.5)
      #xleft, ybottom, xright, ytop,

    }
    #if(any( chr %in% popgen$chr )) {
    #  indx <- which(popgen$chr %in% chr)
    #  abline(v=as.numeric(popgen[indx,'pos1']), col='red')
    #}


  last_coef <- unclass(c2eff)[nrow(c2eff),] # pull out last coefficients

  for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
  }

}


################################################################################
################################################################################

plot_pgen <- function(crs,chrs,stat, map, ahr, ahr_clm, colnm, popgen, ylimo,rank_clm,stat_name,...){

 for (chr in chrs){

  xl <- summary(pull.map(crs))[chr,'length']
  ind <- which(stat$chr == chr)

  Y <- stat[ind,colnm]
  X <- stat[ind,map]/1000000
##  plot(X, Y, col='blue', cex.axis = 2, ylim = ylimo, xlim = c(0,xl), main=paste('CHR',chr), cex.main=2)

  plot(X, Y, col='black',type="n",xlim=c(0,max(X)), ylim = ylimo, main=NULL,
   xlab='physical position', ylab=stat_name, xaxs="i",yaxs="i", mgp = c(1, 0.5, 0),...)

    if(any( chr %in% ahr$chr )) {
      indx <- which(ahr$chr %in% chr)
      #rect(ahr[indx,ahr_clm]/1000000,ylimo[1],ahr[indx,'stp']/1000000,ylimo[2],lwd=0.5,col=alpha('lightgrey',.5))
      abline(v=as.numeric(ahr[indx,ahr_clm])/1000000,
       col='red',lwd=0.5)
    }

    if(any( chr %in% popgen$chr )) {
      indx <- which(popgen$chr %in% chr)
      rect(popgen[indx,'start']/1000000,ylimo[1],popgen[indx,'end']/1000000,ylimo[2],
       border = NA,lwd=0,col=alpha('lightgrey',.5))

      #abline(v=as.numeric(popgen[indx,rank_clm])/1000000, col='grey',lwd=2)
    }
  points(X, Y, col='black',...)
 }
}
################################################################################
################################################################################
################################################################################
################################################################################
## END FUNCTIONS
################################################################################
################################################################################
################################################################################

################################################################################
##pop <- 'NBH'

source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
library("ggridges")
library("plyr")
library("scales")
library("ggrepel")
library('qtl')
library('RColorBrewer')
library('parallel')
mpath <- '/home/jmiller1/QTL_agri/data'
setwd(mpath)

#load(file.path(mpath,paste0(pop,'_scan1_imputed.rsave')))


plogen_data <- function(...){

 ################################################################################

 names(cross$geno) <- ifelse(names(cross$geno) == "X","5",names(cross$geno))
 attr(cross$geno[["5"]], 'class') <- 'A'

 ################################################################################

 scan <- scanone(cross, method = "em", model = "binary", pheno.col = 4)

 ################################################################################

 ################################################################################
 ###### plot_ef rQTL2

 col <- c("slateblue", "violetred", "green3")

 cross2 <- convert2cross2(cross)
 cross2_map <- insert_pseudomarkers(cross2$gmap, step=1)
 cross2_pr <- calc_genoprob(cross2, cross2_map, error_prob=0.025, cores=4)

 ################################################################################
 ################################################################################
 #### AHRs #####

 AHR.bed <- read.table(file.path(mpath,"lift_AHR_genes.bed"), stringsAsFactors = F, header = F)
 colnames(AHR.bed) <- c("chrom", "str", "stp", "gene")
 AHR.bed$chrom <- as.numeric(gsub("chr", "", AHR.bed$chrom))
 AHR.bed$str <- as.numeric(AHR.bed$str)
 AHR.bed$stp <- as.numeric(AHR.bed$stp)
 AHR.notmap <- AHR.bed[is.na(AHR.bed$chrom), ]
 AHR.bed <- AHR.bed[!is.na(AHR.bed$chrom), ]


 AHR.bed <- AHR.bed[which(!AHR.bed$gene == 'EXP'),]
 AHR.bed <- AHR.bed[which(!AHR.bed$gene %in% grep('rna',AHR.bed$gene, value = T)),]


 # add arnts (forgot to scan for them)
 # ahr_genes <- get_AHR(cross)
 cands <- c("AHR1","aip","ARNT","ARNT2","ahrr","ahr1b","AHR2b")
 ################################################

 cross_gens <- cnv.ahrs(cross, AHRdf = AHR.bed, EXP = F)
 ahr <- cross_gens[which(cross_gens$gene %in% cands),]

 ################################################################################
 ### ggplot popgen locations
 if(pop == 'ELR') pops <- 'ER'
 if(pop == 'NBH') pops <- 'NBH'
 cross.popgen <- read.table(file.path(mpath,paste0("outliers",pops,".txt.ncbi.lifted")), sep = "\t", header = T)

 ################################################################################
 ### Use nbh coords but elr and new popgen

 if(pop == 'ELR') top = 200
 if(pop == 'NBH') top = 200

 pop.rank <- cnv.popgen(cross, cross.popgen, top = top)

 ################################################################################

 ## ALL GENES
 genes.bed <- read.table(file.path(mpath,"lifted_genes.bed"), stringsAsFactors = F, header = T)
 genes.bed$chr <- gsub('chr','',genes.bed$chr)
 genes.bed <- genes.bed[genes.bed$chr %in% c(1:24),]
 genes.bed$mid <- round(apply(genes.bed[c('start','end')],1,mean))

 pop_gene_models <- conv_popstat(cross, popgen=genes.bed, whichcol='start',newname='cm_start')
 pop_gene_models$cm_end <- conv_popstat(cross, popgen=genes.bed, whichcol='end',newname='cm_end')[,'cm_end']
 pop_gene_models$cm_mid <- conv_popstat(cross, popgen=genes.bed, whichcol='end',newname='cm_mid')[,'cm_mid']
 ################################################################################

 ################################################################################
 ################################################################################
 pbs <- file.path(mpath, 'pbstat.txt.ncbi.lifted')
 pbs <- read.table(pbs, sep = "\t", header = T)
 pbs$mid <- pbs$V2 + (abs(pbs$V3 - pbs$V2) * .5)
 pbs$V1 <- gsub('chr',"",pbs$V1)
 colnames(pbs)[1:3] <- c('chr','start','end')
 pbs <- pbs[!is.na(as.numeric(pbs$chr)),]
 pbs$pop_cm <- conv_popstat(cross, popgen=pbs, whichcol='mid',newname='pop_cm')$pop_cm

 ################################################################################

 pfst <- file.path(mpath, 'pfst.txt.ncbi.lifted')
 pfst <- read.table(pfst, sep = "\t", header = T)
 pfst$mid <- pfst$start + (abs(pfst$end - pfst$start) * .5)
 pfst$Scaffold <- gsub('chr',"",pfst$Scaffold)
 colnames(pfst)[1] <- 'chr'
 pfst <- pfst[!is.na(as.numeric(pfst$chr)),]
 pfst$pop_cm <- conv_popstat(cross, popgen=pfst, whichcol='mid',newname='pop_cm')$pop_cm

 ################################################################################

 taj <- file.path(mpath, 'tajstat.txt.ncbi.lifted')
 taj <- read.table(taj, sep = "\t", header = T)
 taj$mid <- taj$start + (abs(taj$end - taj$start) * .5)
 taj$Scaffold <- gsub('chr',"",taj$Scaffold)
 colnames(taj)[1] <- 'chr'
 taj <- taj[!is.na(as.numeric(taj$chr)),]
 taj$pop_cm <- conv_popstat(cross, popgen=taj, whichcol='mid',newname='pop_cm')$pop_cm

 ################################################################################

 pi <- file.path(mpath, 'piper.txt.ncbi.lifted')
 pi <- read.table(pi, sep = "\t", header = T)
 pi$mid <- pi$start + (abs(pi$end - pi$start) * .5)
 pi$Scaffold <- gsub('chr',"",pi$Scaffold)
 colnames(pi)[1] <- 'chr'
 pi <- pi[!is.na(as.numeric(pi$chr)),]

 pi$F.NBH <- pi$NBH - pi$F
 pi$BI.NBH <- pi$NBH - pi$BI
 pi$ER.KC <- pi$ER - pi$KC
 pi$ER.SH <- pi$ER - pi$SH

 pi$pop_cm <- conv_popstat(cross, popgen=pi, whichcol='mid',newname='pop_mp')$pop_cm

 ################################################################################
 ################################################################################

 ################################################################################
 ################################################################################
 ################################################################################
 ## Correlate lod and segregation distortion

 c2eff <- lapply(c(1:24),function(X) { scan1coef(cross2_pr[ ,as.character(X)], cross2$pheno[,"bin"]) })
 c2eff <- do.call(rbind, c2eff)
 seg <- geno.table(cross)[rownames(c2eff),'P.value']

 ################################################################################
 ##save.image(file.path(mpath,paste0(pop,'_supplemental_plot_env.rsave')))
 ################################################################################

 sapply(1:24, plogen, stat = 'pfst')
}

plogen_data(pop = 'NBH', cross = cross, mpath = mpath)
