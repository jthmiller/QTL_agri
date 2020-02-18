library('qtl')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
mpath <- '/home/jmiller1/QTL_agri/data'
fl <- paste0(pop,'temp.imp.mapped.tsp.csv')
fl <- file.path(mpath,fl)

source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")


cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)


data <- pull.geno(cross)



gt <- geno.table(cross)

gts <- nind(cross) - gt$missing

HET <- gts*0.5
HOM <- gts*0.25


AA <- HOM - gt$AA
AB <- HET - gt$AB
BB <- HOM - gt$BB


plot_test('ab')
plot(1:length(BB), AB, col = as.factor(gt$chr))
dev.off()
