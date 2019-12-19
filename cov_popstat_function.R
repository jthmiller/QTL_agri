



conv_popstat <- function(cross2, popgen) {

  nbhmap <- convert2cross2(cross2)
  nbhmap$pmap <- nbhmap$gmap
  nbhmap$pmap <- lapply(nbhmap$pmap, function(X) {
    return(as.numeric(gsub("[0-9]+:", "", names(X))))
  })
  for (i in 1:24) {
    names(nbhmap$pmap[[i]]) <- names(nbhmap$gmap[[i]])
  }


  #popgen$chrom <- gsub("chr", "", popgen$chrom)
  #colnames(popgen)[1] <- "chr"

  pop.list <- split(popgen, popgen$chr)
  pop.list <- pop.list[as.character(c(1:24))]
  popgen.pos1 <- lapply(pop.list, "[[", 2)
  popgen.pos2 <- lapply(pop.list, "[[", 3)
  popgen.pos1 <- lapply(popgen.pos1, as.numeric)
  popgen.pos2 <- lapply(popgen.pos2, as.numeric)

  pos1 <- interp_map(popgen.pos1, nbhmap$pmap, nbhmap$gmap)
  pos2 <- interp_map(popgen.pos2, nbhmap$pmap, nbhmap$gmap)

  pop.list <- mapply(cbind, pop.list, pos1 = pos1, SIMPLIFY = FALSE)
  pop.list <- mapply(cbind, pop.list, pos2 = pos2, SIMPLIFY = FALSE)
  fraps <- ldply(pop.list, data.frame, .id = NULL)
  return(fraps)
}
