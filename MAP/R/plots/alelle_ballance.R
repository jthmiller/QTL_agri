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



gt <- geno.table(cross,1:24)
gts <- nind(cross) - gt$missing
HET <- gts*0.5
HOM <- gts*0.25
AA <- gt$AA - HOM
AB <- gt$AB - HET
BB <- gt$BB - HOM
all <- c(AA,AB,BB)
pv <- -log10(gt$P.value)

plot_test('ab', width=5000, height = 5000)
par(mfrow=c(5,5))
for (i in 1:24){
ind <- which(gt$chr == i)
map <- as.numeric(unlist(pull.map(cross,i)))
y <- c(min(all),max(all))
x <- c(min(map),max(map))
plot(x, y, type='n', col = NA, pch=19)
abline(h=0, lwd=5)
points(map, AB[ind], col = 'black', pch=19, cex = 5)
points(map, AA[ind], col = 'red', pch=19, cex = 5)
points(map, BB[ind], col = 'green', pch=19, cex = 5)
}
dev.off()

plot_test('dist', width=1500, height = 1000)
par(mfrow=c(5,5))
for (i in 1:24){
ind <- which(gt$chr == i)
map <- as.numeric(unlist(pull.map(cross,i)))
x <- c(min(map),max(map))
y <- c(min(pv),max(pv)+1)
plot(x, y, type='n', col = NA, pch=19)
points(map, pv[ind], col = 'black', pch=19,cex=1)
text(10,1.5,i, cex=2)
}
dev.off()



### PLOTS ######################################################################
sm <- scanone(cross, pheno.col=4, model="binary",method="mr")
gt <- geno.table(cross)

y <- as.numeric(gsub(".*:","",markernames(cross)))/1000000
x <- 1:length(y)

Y <- c(0, max(y))
X <- c(0, length(y))

ord <- order(as.numeric(gt$chr),as.numeric(gsub(".*:","",rownames(gt))))

lodo <- order(as.numeric(sm$chr))
lod <- sm$lod[lodo]

plot_test('nbh_mar_regression_after_map', width = 5500, height = 750)
par(mfrow=c(3,1))
 plot(1:length(sm$lod), lod, pch = 19, col = factor(sm$chr[lodo]), ylim = c(0,18), cex = 0.25)
 plot(1:length(gt[,1]), -log10(gt[ord,'P.value']), pch = 19, col = factor(sm$chr[ord]), ylim = c(0,2), cex = 0.25)
 abline(h=6)
 plot(X,Y,type="n", xlab=paste('chr',i), ylab='physical position')
  points(x,y[ord], col = as.factor(gt$chr[ord]))
dev.off()
################################################################################
