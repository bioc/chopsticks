# Routines for structure analyses

xxt <- function(snps, strata=NULL,
                correct.for.missing=FALSE, lower.only=FALSE) {
  if (!is.null(strata) && !is.factor(strata))
    strata <- as.factor(strata)
  .Call("xxt", snps, strata, correct.for.missing, lower.only,
        PACKAGE="snpMatrix")
}

ibsCount <- function(snps) {
  .Call("ibs_count", snps, PACKAGE="snpMatrix")
}

ibsDist <- function(counts) {
  .Call("ibs_dist", counts, PACKAGE="snpMatrix")
}

snp.pre.multiply <- function(snps,  mat, frequency=NULL) {
  .Call("snp_pre", snps, mat, frequency, PACKAGE="snpMatrix")
}

snp.post.multiply <- function(snps,  mat, frequency=NULL) {
  .Call("snp_post", snps, mat, frequency, PACKAGE="snpMatrix")
}

snp.cor <- function(x, y) {
  .Call("corsm", x, as.matrix(y), PACKAGE="snpMatrix")
}
