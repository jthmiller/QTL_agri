#!/bin/R
### first run combine pops for multi-pop cross objects

load(file.path(mpath,'supplemental_plot_env.rsave'))

################################################################################
#bottom, left, top and right
#location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks


plogen <- function(ch, stat = 'pfst', crossA, erprob = 0.001, qtl = hk.qtl.4_perms){

 ##load(file.path(mpath,paste0(pop,'_supplemental_plot_env.rsave')))

 ################################################################################
 ################################################################################

 cross <- sim.geno(crossA, n.draws = 150, error.prob = erprob, map.function="kosambi", stepwidth="fixed")
 cross <- calc.genoprob(cross, error.prob = erprob, map.function="kosambi", stepwidth="fixed")

 erp <- 0.001
 norm <- scanone(cross, method = "imp", model = "normal", pheno.col = 5)
 ## bin <- scanone(cross, method = "hk", model = "binary", pheno.col = 4)
 bin <- scanone(cross, method = "em", model = "binary", pheno.col = 4)
 gt <- geno.table(cross)
 map <- pull.map(cross)
 map_sum <- summary(pull.map(cross))

 ################################################################################

 col <- c("slateblue", "violetred", "green3")
 cross2 <- convert2cross2(cross)
 map2 <- insert_pseudomarkers(cross2$gmap, step=0.5)
 pr2 <- calc_genoprob(cross2, map2, error_prob=0.0025, cores=cores)
 bin2 <- scan1(pr2, pheno=cross2$pheno[,'bin'] , model = "binary", cores = cores)
 peaks <- find_peaks(bin2, map2)

 cross2$pmap <- cross2$gmap
 cross2$pmap <- lapply(cross2$pmap, function(X) {
   return(as.numeric(gsub("[0-9]+:", "", names(X))))
 })
 for (i in chrnames(cross2)) {
   names(cross2$pmap[[i]]) <- names(cross2$gmap[[i]])
 }

 ################################################################################

 ##rank <- pop.rank

 if(pop == 'NBH') { pbsname <- 'NBH'; pfstNSname <- 'F.NBH' ; piname <- pfstname <- 'BI.NBH' ; rank <- nbh.rank ; ahr <- ahr_nbh }
 if(pop == 'ELR') { pbsname <- 'ER'; pfstNSname <- 'ER.SH' ; piname <- pfstname <- 'ER.KC' ; rank <- elr.rank ; ahr <- ahr_elr }

 ## AHR GENES IN MAP2 SPACE
 pos <- mapply(mean, ahr$str, ahr$stp)
 AHR.list <- split(pos, as.numeric(ahr$chr))
 AHR.pos <- lapply(AHR.list, as.numeric)
 AHR.pos <- interp_map(AHR.pos,  cross2$pmap, cross2$gmap)
 AHR.name <- split(ahr$gene, as.numeric(ahr$chr))

 print(paste('plotting pfst for',pfstNSname))

################################################################################
pdf(paste0("/home/jmiller1/public_html/",pop,"_all",ch,"segdist.pdf"), width=4.5,height=6)
 mat<-matrix(c(1:8),8,1, byrow=T)
 layout(mat, widths=1, heights= c(0.1, 0.1, 0.1, 0.025, 0.1, 0.1, 0.1, 0.05))
 par(mar=c(0.5,2.25,0,1)+0.1, oma = c(1, 1, 0, 1))

 len <- map_sum[ch,'length']
 ind <- which(gt$chr == ch)
 segX <- map[[ch]][rownames(gt)[ind]]
 segY <- -log10(gt[ind,'P.value'])
 out <- which(!gt$chr == ch)
 segB <- -log10(gt[out,'P.value'])
 segA <- rescale(unlist(map, use.names=F),to = c(0,len))[out]
 ab <- ahr[which(ahr$chr==ch),'pos1']

 ################################################################################
 plot_pgen(
  crs = cross, chrs = ch, stat = pfst, map = 'mid', mgp = c(1.25, 0.5, 0),
  ahr = ahr, ahr_clm = 'stp',  colnm = pfstNSname, popgen = rank,
  rank_clm='end', ylimo=c(-0.025, 1), stat_name= paste0('pfst NxS (',pfstNSname,')'), pch=16,
  cex=0.15, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75, xaxt="n"
 )

 if(stat == 'pfst'){

  print(paste('plotting pfst for',pfstNSname))

  plot_pgen(
   crs = cross, chrs=ch, stat = pfst, map = 'mid', mgp = c(1.25, 0.5, 0),
   ahr = ahr, ahr_clm = 'stp',  colnm = pfstname, popgen = rank,
   rank_clm = 'end', ylimo = c(-0.025,1), stat_name = 'pfst', pch=16,
   cex=0.15, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75, xaxt="n"
  )
 } else {
  plot_pgen(
   crs = cross, chrs=ch, stat = pbs, map = 'mid', mgp = c(1.25, 0.5, 0),
   ahr = ahr, ahr_clm= 'stp',  colnm = pbsname, popgen = rank,
   rank_clm='end', ylimo=c(-0.08,0.5), stat_name='pbs', pch=16,
   cex=0.15, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75, xaxt="n"
  )
 }

 plot_pgen(
  crs = cross, chrs=ch, stat = pi, map = 'mid', mgp = c(1.25, 0.5, 0),
  ahr = ahr, ahr_clm= 'stp', colnm = piname, popgen = rank,
  rank_clm = 'end', ylimo=c(-0.03,0.03), stat_name='delta pi', pch=16,
  cex=0.15, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75
 )

 ## spacer empty plot ## ## ## ## ##
 plot(0,type='n',axes=FALSE,ann=FALSE)
 ## ## ## ## ## ## ## ## ## ## ## ##

 ## CM
 plot(bin2, map,chr=ch, lodcolumn=1, ylim=c(0,8),xaxt="n", mgp = c(1.25, 0.5, 0))
 ## add qtls
 if (any(qtl$chr==ch)) {
  #bi <- bayesint(bin, chr = ch, prob=0.95, lodcolumn=1, expandtomarkers=TRUE)
  peaks <- peaks[which(peaks$chr == ch),]
  ah.pos <- unname(unlist(AHR.pos[as.character(ch)]))
  ah.gen <- unname(unlist(AHR.name[as.character(ch)]))
  ind <- which(qtl$chr==ch)
  ## 1.5-LOD support interval
  bi <- bayes_int(bin2,  map2, chr = ch, expand2markers=TRUE)
  rect(bi[1], 0, bi[3], 8, col = alpha('lightgrey',.5), border = 'black', lwd=0.5)
  ## sapply(ah.pos, function(x) { abline(v = x, col='red',lwd=0.5) })
  abline(v = ah.pos, col='red',lwd=0.5)
  abline(v = bi[2], col='green',lwd=0.5)
  ## sapply(ah.pos, function(x) { text(x = x, 2, labels = ah.gen) })

 }

 if (any(ahr$chr==ch)) {
  ah.pos <- unname(unlist(AHR.pos[as.character(ch)]))
  ah.gen <- unname(unlist(AHR.name[as.character(ch)]))
  abline(v = ah.pos, col='red',lwd=0.5)
  text(x = ah.pos, 2, labels = ah.gen)
 }

 plot_ef(
  crs = cross2, map = map, pr = pr2 , ahr = ahr,
  popgen = rank, chs=ch, main=NULL, model="bin", xaxt="n",
  cex.axis = 0.5, cex.lab= 0.75, cex.main = 0.75
 )

 effectscan(
  cross, pheno.col=4, chr=ch, get.se=T, draw=TRUE,
  gap=25, mtick="line",add.legend=F, alternate.chrid=T, ylim=c(-0.5,0.5),
  xlim=c(0,len), main=NULL,ylab='Effect Est. (+/- 1 SE)', xaxs="i",, xaxt="n",
  cex.axis = 0.75, cex.lab= 0.75, cex.main = 0.75, mgp = c(1.5, 0.5, 0)
 )
 abline(v=ab,col='red',lwd=0.5)

 plot(segA,segB, ylim=c(0,4), xlim=c(0,max(segX)), pch=16, col='lightgrey',
  xlab="",ylab='-log10 pval',xaxs="i", mgp = c(1.25, 0.5, 0),
  cex=0.25, cex.axis = 0.75, cex.lab=0.75, cex.main=0.75)

 points(segX,segY,pch=16,cex=0.25)

 abline(v=ab,col='red',lwd=0.5)
 dev.off()
}
################################################################################

cores <- 12

plogen(18, crossA = cross_NBH, stat = 'pfst')
plogen(2, crossA = cross_NBH, stat = 'pfst')
plogen(8, crossA = cross_NBH, stat = 'pfst')

sapply(1:24, plogen, crossA = cross_NBH, stat = 'pfst')
sapply(1:24, plogen, crossA = cross_ELR, stat = 'pfst')
################################################################################


for i in {1..24} ; do
scp jmiller1@farm.cse.ucdavis.edu:/home/jmiller1/public_html/NBH_all2segdist.pdf .
done
