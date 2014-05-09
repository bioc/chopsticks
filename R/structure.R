# Routines for structure analyses

xxt <- function(snps, correct.for.missing=FALSE, lower.only=FALSE) {
  .Call("xxt", snps, correct.for.missing, lower.only,
        PACKAGE="chopsticks")
}

ibsCount <- function(snps) {
  .Call("ibs_count", snps, PACKAGE="chopsticks")
}

ibsDist <- function(counts) {
  .Call("ibs_dist", counts, PACKAGE="chopsticks")
}

snp.pre <- function(snps,  mat, frequency=NULL) {
  .Call("snp_pre", snps, mat, frequency, PACKAGE="chopsticks")
}

snp.post <- function(snps,  mat, frequency=NULL) {
  .Call("snp_post", snps, mat, frequency, PACKAGE="chopsticks")
}

snp.cor <- function(x, y) {
  .Call("corsm", x, as.matrix(y), PACKAGE="chopsticks")
}
