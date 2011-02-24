snp.compare <- function(obj1, obj2) {
  if((mode(obj1) != "raw") || (mode(obj2) != "raw"))
    stop("non-raw input(s)\n")
  if(!all.equal(rownames(obj1), rownames(obj2)) || !all.equal(colnames(obj1), colnames(obj2)))
    stop("row/col names not the same\n")
  
  snps.names <- colnames(obj1)
  snp.count <- length(snps.names)
  sample.count <- dim(obj1)[1] #rownames(obj1)
  result <- .C('count_gt',
              as.raw(obj1@.Data), as.raw(obj2@.Data),
              as.integer(snp.count),
              as.integer(sample.count),
              count = integer(snp.count),
              count.signed = integer(snp.count),
              PACKAGE="chopsticks"
              )[c('count', 'count.signed')]
  names(result$count) <- snps.names
  names(result$count.signed) <- snps.names
  result
}

.label.clusters.by.gtype <- function(cluster, gtype) {
  aa.select <- gtype == 0x01
  ab.select <- gtype == 0x02
  bb.select <- gtype == 0x03
  nc.select <- gtype == 0x00
  list(aa=cluster[aa.select,,drop=FALSE],
       ab=cluster[ab.select,,drop=FALSE],
       bb=cluster[bb.select,,drop=FALSE],
       nc=cluster[nc.select,,drop=FALSE])
}

.plot.clust <- function(clus.list, title="test") {
  grid.quantile <- quantile(c(clus.list$aa, clus.list$ab, clus.list$bb, clus.list$nc),
                            probs = c(0.5, 0.97, 1), na.rm=TRUE)
  max.xy <- min(grid.quantile[['100%']], (grid.quantile[['97%']] * 2 - grid.quantile[['50%']]))
  plot(clus.list$aa, col='red', xlim=c(0, max.xy), ylim=c(0, max.xy), main=title)
  points(clus.list$ab, col='blue')
  points(clus.list$bb, col='yellow4')
  points(clus.list$nc, col='gray15')
}

snp.clust.plot <- function(cluster, gtype, title="test") {
  coloured.cluster <- .label.clusters.by.gtype(cluster, gtype)
  .plot.clust(coloured.cluster, title=title)
}
