#!/bin/R
################################################################################
pop <- 'NBH'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales")
suppressMessages(sapply(libs2load, require, character.only = TRUE))

library(doParallel)
################################################################################
##fl <- "NBH_4822_imputed_NW_tsp.csv"
#fl <- "NBH_5755_imputed_NW_tsp.csv"
#cross <- read.cross(file=fl , format = "csv", dir=mpath,
# genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

################################################################################
### IMPUTED
fl <- 'NBH_8498_imputed_high_confidence_tsp_mapped.csv'
cross <- read.cross(file=fl , format = "csv", dir=mpath,
 genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
#################################################################################


## High missing
drop <- 'NBH_6137'
cross <- subset(cross, ind = !cross$pheno$ID == drop)
#################################################################################

drop <- c('NW_012234461.1:562796','NW_012234461.1:1059423','NW_012234494.1:872524','8:12257474','NW_012234558.1:739273','NW_012234311.1:628953','NW_012234311.1:166994','NW_012234311.1:628953','NW_012234311.1:166994','NW_012234520.1:169476','NW_012225526.1:51452','NW_012234326.1:397833','NW_012234326.1:380172')
cross <- drop.markers(cross,drop)

drop2 <- c('10:9960308','NW_012234311.1:3468744', 'NW_012234311.1:3445534','NW_012234311.1:3383007','NW_012234326.1:2162519','NW_012224981.1:180658','NW_012224824.1:3425')
cross <- drop.markers(cross,drop2)





###
NWs <- grep('NW',markernames(cross), value=T)
crossX <- drop.markers(cross, NWs)

crossX <- fill.geno(crossX, method="maxmarginal", error.prob = 0.05, min.prob=0.95)
gt <- geno.table(crossX)
toss <- rownames(gt[which(gt[,'P.value'] < 3.16e-4),])
crossX <- drop.markers(crossX,toss)

## CANNOT FILTER ON 2 AND 13 ###################################################
gt <- geno.table(crossX)
gt <- gt[! gt$chr %in% c(2,13),]
toss <- rownames(gt[which(gt[,'P.value'] < 1.0e-2),])
crossX <- drop.markers(crossX,toss)
################################################################################

### ANALYYTICAL PLOT ###########################################################
gt <- geno.table(crossX)
plot_test('missing.nbh')
hist(gt$missing, breaks=80)
dev.off()
################################################################################

################################################################################
drop <- names(which(colSums(is.na(pull.geno(crossX))) > 10))
crossX <- drop.markers(crossX,drop)
pull <- markernames(crossX)
cross <- pull.markers(cross, pull)
################################################################################

### READ THE UNMAPPED MARKER ASSIGNMENT TABLE ##################################
movefl <- file.path(mpath,'NBH_NW_scaffold_assignments.tsv')
move <- read.table(movefl, stringsAsFactors = F, header=T, sep = " ")
rownames(move) <- move$nw_marks_assign
move <- move[which(move$nw_marks_assign %in% markernames(cross)),]
################################################################################

### ASSIGN MARKERS TO GROUPS 1-24 ##############################################
for(i in move$nw_marks_assign){
 gi <- pull.geno(cross)[,i]
# add marker to cross
 crossX <<- addmarker(crossX, gi, i, move[i,'nw_ch'], move[i,'nw_pos'])
}
cross <- subset(crossX,chr=1:24)
################################################################################

## CANNOT FILTER ON 2 AND 13 ###################################################
gt <- geno.table(cross)
gt <- gt[! gt$chr %in% c(2,13),]
toss <- rownames(gt[which(gt[,'P.value'] < 1.0e-2),])
cross <- drop.markers(cross,toss)
################################################################################

################################################################################
### TOSS MARKERS WITH HIGH PERCENTAGE OF MISSING DATA ##########################
misg <- function(X,perc) { nind(cross) * perc }
mis <- misg(cross,0.10)
drop <- names(which(colSums(is.na(pull.geno(cross))) > mis))
cross <- drop.markers(cross,drop)
################################################################################

################################################################################
### REORDER MARKERS
cross <- est.rf(cross,maxit=1000)
cross <- tspOrder(cross = cross, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
save.image(file.path(mpath,paste0(pop,'_estrf_no_impute.rsave')))
################################################################################


################################################################################
### REORDER MARKERS
cross <- est.rf(cross,maxit=1000)
cross <- tspOrder(cross = cross, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
save.image(file.path(mpath,paste0(pop,'_estrf_no_impute.rsave')))
################################################################################

### ESTIMATE ERROR RATE AND MAP ################################################
loglik <- err <- c(0.0001, 0.0005, 0.001, 0.01, 0.05)

cross.sub <- subset(cross,chr=c(1,3,4,7,9,12,15,17,20,23))

update.lik <- function(z){
  cat(z, "of", length(err), "\n")
  tempmap <- est.map(cross.sub,maxit=1000, error.prob=err[z])
  loglik[z] <- sum(sapply(tempmap, attr, "loglik"))
}

cl <- makeCluster(5)
registerDoParallel(cl)

loglik <- foreach(z = seq(along=err), .inorder = T, .export = c("loglik"), .packages = c("qtl")) %dopar% update.lik(z)

loglik <- unlist(loglik)

lod <- (loglik - max(loglik))/log(10)

erprob <- err[which.max(lod)]

print(paste('error lod =',erprob))

##erprob <- 1e-04
par.est.map <- function(X) { est.map(cross, chr = X, error.prob=erprob, map.function="kosambi",maxit=1000, tol=1e-7, sex.sp=FALSE, verbose=FALSE) }

### MAP ########################################################################
cl <- makeCluster(20)
registerDoParallel(cl)
maps <- foreach(X = 1:24, .inorder = T, .packages = libs2load) %dopar% par.est.map(X)
################################################################################

maps <- lapply(maps,"[[",1)
names(maps) <- c(1:24)
attr(maps, 'class') <- 'map'

cross <- qtl:::replace.map(cross,maps)

direc <- sapply(1:24,function(i) {
 pos <- as.numeric(gsub(".*:","",markernames(cross,i)))
 map <- as.numeric(pull.map(cross)[[i]])
 cor(pos,map, use="complete.obs")
})

if(any(direc < 0)) cross <- flip.order(cross,which(direc < 0))

mapfile <- paste0(pop,'_',sum(nmar(cross)),'_imputed_tsp_mapped')
filename <- file.path(mpath,mapfile)
write.cross(cross,filestem=filename,format="csv")

print(paste(filename, 'cross written'))

################################################################################
save.image(file.path(mpath,paste0(pop,'_estrf_no_impute.rsave')))
################################################################################

################################################################################
png(paste0('~/public_html/',pop,'high_conf_no_imput_rf.png'),height=2500,width=2500)
par(mfrow=c(6,4))
for (B in 1:24){
  plotRF(cross,chr=B)
}
dev.off()
################################################################################

################################################################################
plot_test('nbh_high_confid_rf')
 plotRF(cross,zmax=10,col.scheme="redblue")
dev.off()
################################################################################
################################################################################
