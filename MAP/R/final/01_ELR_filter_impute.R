#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','BRP','NEW','ELR','ELR.missing')]

source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")

mpath <- '/home/jmiller1/QTL_agri/data'

libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
library(scales)

################################################################################
## read in the QTL cross
cross <- read.cross(file = file.path(mpath, paste0(pop, ".unphased.f2.csvr")),
format = "csvr", geno = c(1:3), estimate.map = FALSE)
################################################################################

################################################################################
### Pull names from plinkfile
path <- file.path(mpath, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross$pheno$ID <- paste(popname, indname, sep = "_")
################################################################################

#### PHENO #####################################################################
cross$pheno$bin <- ifelse(cross$pheno$Pheno > 2, 1 , 0)
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno))
################################################################################

#### SEX #######################################################################
sex <- read.table(file.path(mpath,'sex.txt'),stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross$pheno$ID), 'sex']
cross$pheno$sex <- sex.vec
################################################################################

toss.badata <- c("ELR_10869","ELR_10987","ELR_11580")
cross <- subset(cross,ind=!cross$pheno$ID %in% c(toss.badata,'BLI_BI1124M','ELR_ER1124F'))

### TOSS MARKERS WITH HIGH PERCENTAGE OF MISSING DATA ##########################
misg <- function(X,perc) { nind(cross) * perc }
mis <- misg(cross,0.10)
drop <- names(which(colSums(is.na(pull.geno(cross))) > mis))
cross <- drop.markers(cross,drop)
################################################################################

### READ THE UNMAPPED MARKER ASSIGNMENT TABLE
movefl <- file.path(mpath,'NBH_NW_scaffold_assignments.tsv')
move <- read.table(movefl, stringsAsFactors = F, header=T, sep = " ")
move <- move[which(move$nw_marks_assign %in% markernames(cross4)),]

### ASSIGN UNMAPPED MARKERS
cross5 <- cross4
for (i in 1:length(move[,1])){
 cross5 <<- movemarker(cross5, marker = move[i,'nw_marks_assign'], newchr = move[i,'nw_ch'], newpos = as.numeric(move[i,'nw_pos']))
}
cross5 <- subset(cross5,chr=1:24)
################################################################################




### PLOTS ######################################################################
sm <- scanone(cross, pheno.col=4, model="binary",method="mr")
Y <- c(0, as.numeric(gsub(".*:","",markernames(cross))))/1000000
X <- 1:length(Y)
gt <- geno.table(cross)
plot_test('elr_mar_regression', width = 5500, height = 750)
par(mfrow=c(3,1))
 plot(1:length(sm$lod), sm$lod, pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
 plot(1:length(gt[,1]), -log10(gt[,'P.value']), pch = 19, col = factor(sm$chr), ylim = c(0,18), cex = 0.25)
 abline(h=6)
 plot(c(1,length(X)),c(0,max(Y)),type="n", ylab='physical position')
  points(X,Y)
dev.off()
################################################################################

mapit_noimpute <- function(i){

 erprob <- 0.05
 Z <- i

 cross1 <- subset(cross,chr=Z)
 cross1 <- use_phys_map(cross1)
 cross1 <- est.rf(cross1)

 ### MAKE LOW/NO RECOMB GROUPS FOR CLEANUP #####################################
 ## high LOD initial check phase ###############################################
 RF <- 3/nind(cross1)
 LOD <- 15
 cross2 <- formLinkageGroups(cross1, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)

 ### REMOVE NON AB AB ###################################################
 dist <- sapply(chrnames(cross2), function(X) { mean(-log10(geno.table(cross2, chr=X)$P.value)) })
 keep <- names(which(dist < 4))
 cross2 <- subset(cross2,chr=keep)
 ####################################################################

 ### PVAL filt #################################################################
 gt <- geno.table(cross2)
 pval <- 1.0e-5
 mis <- misg(cross2,0.15)
 bfixA <- rownames(gt[which(gt$P.value > pval & gt$missing < mis),])
 cross2 <- pull.markers(cross2,bfixA)
 ###############################################################################

 ## Switch phase of groups that are out of phase with LG1
 rf <- pull.rf(est.rf(cross2))
 chr <- chrnames(cross2)[-1]
 rf.mean <- sapply(chr, function(X) { mean(rf[markernames(cross2,chr=1),markernames(cross2,chr=X)],na.rm=T) })
 flips <- names(which(rf.mean > 0.5))
 cross2 <- switchAlleles(cross2, markers = markernames(cross2,chr=flips))

 ## RETURN ALL MARKERS TO THE SAME LG ##########################################
 RF <- 8/nind(cross2)
 LOD <- 8
 cross2 <- formLinkageGroups(cross2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
 cross2 <- subset(cross2, chr=chrnames(cross2)[1])
 ###############################################################################

 cross2 <- use_phys_map(cross2)
 ord <- order(as.numeric(unlist(pull.map(cross2))))
 cross2 <- switch.order(cross2, chr = 1, ord, error.prob = 0.01, map.function = "kosambi", maxit = 1, tol = 0.1, sex.sp = F)

 cross2 <- removeDoubleXO(cross2)
 cross2 <- fill.geno(cross2, method="no_dbl_XO", error.prob = 0.05, min.prob=0.995)

 ## Take the best marker every 1kb #############################################
 cross2 <- thin_by_distortion(cross2,10)
 ###############################################################################

 #############################
 ## REMOVE MARKERS WITH HIGH MISSING DATA
 mis <- misg(cross2,0.20)
 bfixA <- names(which(colSums(is.na(pull.geno(cross2))) > mis))
 print(paste('dropped',length(bfixA),'markers due to missing data'))
 cross2 <- drop.markers(cross2, bfixA)
 #############################
 pos <- as.numeric(gsub(".*:","",markernames(cross2)))
 map <- as.numeric(pull.map(cross2)[[1]])

 if(cor(pos,map, use="complete.obs") < 0) cross2 <- flip.order(cross2, 1)
 cross2 <- shiftmap(cross2, offset=0)

 mapfile <- paste0(pop,'_unmapped_noimput_',i,'_tsp')
 filename <- file.path(mpath,mapfile)
 write.cross(cross2,filestem=filename,format="csv")

 cross3 <- fill.geno(cross2, method="maxmarginal", error.prob = 0.05, min.prob=0.98)
 cross3 <- fill.geno(cross3, method="no_dbl_XO", error.prob = 0.05, min.prob=0.98)
 cross3 <- tspOrder(cross = cross3, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
 pos <- as.numeric(gsub(".*:","",markernames(cross3)))
 map <- as.numeric(pull.map(cross3)[[1]])
 if(cor(pos,map, use="complete.obs") < 0) cross3 <- flip.order(cross3, 1)
 cross3 <- shiftmap(cross3, offset=0)
 mapfile <- paste0(pop,'_imputed_',i,'_tsp')
 filename <- file.path(mpath,mapfile)
 write.cross(cross3,filestem=filename,format="csv")

 cross4 <- tspOrder(cross = cross2, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
 mapfile <- paste0(pop,'_reordered_noImp_',i,'_tsp')
 filename <- file.path(mpath,mapfile)
 write.cross(cross4,filestem=filename,format="csv")

 plotit(cross2,nme='no_imp')
 plotit(cross3,nme='imputed')
 plotit(cross4,nme='reorder_noimp')
}


library(qtl)
library(doParallel)
cl <- makeCluster(5)
registerDoParallel(cl)
foreach(i = 1:24, .inorder = F, .packages = libs2load) %dopar% mapit_noimpute(i)

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
 swits <- markernames(cross2, chr=chrnames(cross2)[1])
 cross2 <- switchAlleles(cross2, markers = swits)
 cross2 <- formLinkageGroups(cross2, max.rf = RF, min.lod = LOD, reorgMarkers = TRUE)
