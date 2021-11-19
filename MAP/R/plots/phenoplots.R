pop <- 'NBH'
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap","scales","doParallel")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
mpath <- '/home/jmiller1/QTL_agri/data'

cross <- elr_cross
qtlscan <- elr_sone_bin
chr <- 18

plot_pheno_effect <- function(cross, qtlscan, chr, pltnme){

 qtlpos <- summary(qtlscan)[chr,]

 mark <- find.marker(cross, chr,  qtlpos$pos)

 genos <- pull.geno(cross, chr)[,mark]

 pheno <- pull.pheno(cross,4)

 gp <- cor(genos, pheno)

 if (gp < 0 ) genos <- factor(genos, labels = c('SS','RS','RR'))

 if (gp > 0 ) genos <- factor(genos, labels = c('RR','RS','SS'))

 pheno <- ifelse(pheno == 1,0,1)

 gnt <- c('RR','RS','SS')

 propr <- sapply(gnt,function(x){
    tab <- table(factor(pheno[which(genos == x)], levels = c(0,1)))
    answer <- tab['1']/(tab['0']+tab['1'])
    return(ifelse(answer == 'Inf', 1, answer))
 })

 names(propr) <- gnt
 genos <- as.numeric(factor(genos, levels = gnt))


 plot_test(pltnme)
 plot(x = c(0.75,3.25), y = c(0,1), ylim=c(-0.1,1.1), type = 'n', lwd = 5, yaxt="n", xaxt="n", ylab = '', xlab = '')
 ## title( ylab = 'Proportion Resistant', xlab = 'Embryo Genotype', cex = 2 )
 points(jitter(genos, amount = 0.1), jitter(pheno, amount = 0.1), pch = 19)
 lines(x = c(1,2,3) + 0.1, y = as.numeric(propr), lwd = 5)
 points(1:3 + 0.1, as.numeric(propr) , pch=22, lwd = 3, bg = 'white', cex = 4)
 box(lwd=3)
 axis(1, at=c(1,2,3),labels=c('RR','RS','SS'), las=0, cex.axis=2, lwd = 2, line = 1)
 axis(2, at=c(0,0.5,1),labels=c(0,0.5,1), las=2, cex.axis=2, lwd = 2)
 axis(4, at=c(0,1),labels=c('S','R'), las=2, cex.axis=2, lwd = 2)
 dev.off()
}

brp_mapfile <- 'BRP_reorder_imp_mapped'
new_mapfile <- 'NEW_reorder_imp_mapped'
nbh_mapfile <- 'NBH_reorder_imp_nopar'
elr_mapfile <- 'ELR_reorder_imp_mapped'

nbh_cross <- read.cross(file = paste0(nbh_mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
new_cross <- read.cross(file = paste0(new_mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
elr_cross <- read.cross(file = paste0(elr_mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
brp_cross <- read.cross(file = paste0(brp_mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

## scan with marker regression
nbh_sone_bin <- scanone(nbh_cross, pheno.col=4, method="mr", model="bin", n.cluster = cores)
elr_sone_bin <- scanone(elr_cross, pheno.col=4, method="mr", model="bin", n.cluster = cores)
brp_sone_bin <- scanone(brp_cross, pheno.col=4, method="mr", model="bin", n.cluster = cores)
new_sone_bin <- scanone(new_cross, pheno.col=4, method="mr", model="bin", n.cluster = cores)

## plot chr2
plot_pheno_effect(cross = nbh_cross, qtlscan = nbh_sone_bin, chr = 2, pltnme = 'nbh_chr2_phenoplot')
plot_pheno_effect(cross = new_cross, qtlscan = new_sone_bin, chr = 2, pltnme = 'new_chr2_phenoplot')
plot_pheno_effect(cross = brp_cross, qtlscan = brp_sone_bin, chr = 2, pltnme = 'brp_chr2_phenoplot')
plot_pheno_effect(cross = elr_cross, qtlscan = elr_sone_bin, chr = 2, pltnme = 'elr_chr2_phenoplot')

## plot chr18
plot_pheno_effect(cross = nbh_cross, qtlscan = nbh_sone_bin, chr = 18, pltnme = 'nbh_chr18_phenoplot')
plot_pheno_effect(cross = new_cross, qtlscan = new_sone_bin, chr = 18, pltnme = 'new_chr18_phenoplot')
plot_pheno_effect(cross = brp_cross, qtlscan = brp_sone_bin, chr = 18, pltnme = 'brp_chr18_phenoplot')
plot_pheno_effect(cross = elr_cross, qtlscan = elr_sone_bin, chr = 18, pltnme = 'elr_chr18_phenoplot')

## plot chr13
plot_pheno_effect(cross = nbh_cross, qtlscan = nbh_sone_bin, chr = 13, pltnme = 'nbh_chr13_phenoplot')
plot_pheno_effect(cross = new_cross, qtlscan = new_sone_bin, chr = 13, pltnme = 'new_chr13_phenoplot')
plot_pheno_effect(cross = brp_cross, qtlscan = brp_sone_bin, chr = 13, pltnme = 'brp_chr13_phenoplot')
plot_pheno_effect(cross = elr_cross, qtlscan = elr_sone_bin, chr = 13, pltnme = 'elr_chr13_phenoplot')

## plot chr11 (newark may have minor there)
plot_pheno_effect(cross = nbh_cross, qtlscan = nbh_sone_bin, chr = 11, pltnme = 'nbh_chr11_phenoplot')
plot_pheno_effect(cross = new_cross, qtlscan = new_sone_bin, chr = 11, pltnme = 'new_chr11_phenoplot')
plot_pheno_effect(cross = brp_cross, qtlscan = brp_sone_bin, chr = 11, pltnme = 'brp_chr11_phenoplot')
