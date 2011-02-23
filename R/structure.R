# Routines for structure analyses

xxt <- function(snps, correct.for.missing=FALSE, lower.only=FALSE) {
  .Call("xxt", snps, correct.for.missing, lower.only,
        PACKAGE="snpMatrix")
}

ibsCount <- function(snps) {
  .Call("ibs_count", snps, PACKAGE="snpMatrix")
}

ibsDist <- function(counts) {
  .Call("ibs_dist", counts, PACKAGE="snpMatrix")
}

snp.pre <- function(snps,  mat, frequency=NULL) {
  .Call("snp_pre", snps, mat, frequency, PACKAGE="snpMatrix")
}

snp.post <- function(snps,  mat, frequency=NULL) {
  .Call("snp_post", snps, mat, frequency, PACKAGE="snpMatrix")
}

snp.cor <- function(x, y) {
  .Call("corsm", x, as.matrix(y), PACKAGE="snpMatrix")
}
