snp.imputation<- function(X, Y, pos.X, pos.Y, phase=FALSE, try=50,
                          r2.stop = 0.9, max.X=4) {
  if (phase)
    stop("phase=TRUE option not yet implemented")
  if (!(is(X, "snp.matrix") && is(Y, "snp.matrix")))
    stop("X and Y arguments must be of class snp.matrix or X.snp.matrix")
  if (class(X)[1]!=class(Y)[1])
    stop("X and Y arguments must be of the same class")
  if (nrow(X)!=nrow(Y) || length(pos.X)!=ncol(X) || length(pos.Y)!=ncol(Y))
    stop("Inconsistent argument dimensions")
  order.X <- order(pos.X)
  order.Y <- order(pos.Y)
  if (try>ncol(X))
    try <- ncol(X)
  if (max.X>try)
    max.X <- try
  .Call("snp_impute", X, Y, order.X, order.Y,
        as.double(pos.X[order.X]), as.double(pos.Y[order.Y]),
        as.logical(phase), as.integer(try), as.double(r2.stop),
        as.integer(max.X), PACKAGE="snpMatrix")
}

impute.snps <- function(rules, snps, subset=NULL) {
  if (class(rules)!= "snp.reg.imputation")
    stop("incorrect class for `rules' argument")
  if (class(snps)!= "snp.matrix" && class(snps)!="X.snp.matrix")
    stop("incorrect class for `snps' argument")
  if (!is.null(subset)) {
    nr <- nrow(snps)
    if (is.logical(subset)) {
      if (length(subset)!=nr)
        stop("incompatible dimensions for snps and subset")
      subset <- (1:nr)[subset]
    }
    subset <- as.integer(subset);
    if (any((subset<1) | (subset>nr)))
      stop("row number out of range in subset argument")
  }
  .Call("impute_snps", rules, snps, subset, PACKAGE="snpMatrix")
}

# Returns names of rules (imputed snps) which feature excluded SNPs

filter.rules <- function(rules, snps.excluded, exclusions=TRUE) {
  if (class(rules)!= "snp.reg.imputation")
    stop("incorrect class for `rules' argument")
  if (!is.character(snps.excluded))
    stop("excluded SNPs must be referenced by name")
  if.excl <- sapply(rules, FUN=function(rule, exc){
                             any(rule$snps %in% exc)
                           }, snps.excluded)
  if (exclusions)
    names(rules)[if.excl]
  else
    names(rules)[!if.excl]
}

# List extraction functions

imputation.maf <- function(rules) sapply(rules, function(x) x$maf)

imputation.r2 <- function(rules) sapply(rules, function(x) x$r.squared)
