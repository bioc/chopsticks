test.allele.switch <- function(snps, snps2=NULL, split=NULL, prior.df=1) {
  if ((class(snps)!="snp.matrix") && (class(snps)!="X.snp.matrix"))
    stop("argument `snps' is not a snp.matrix or X.snp.matrix")
  if (is.null(snps2) && is.null(split))
    stop("either second snp.matrix or split vector must be specified")
  if (is.null(snps2)) {
    if (length(split) != nrow(snps)) 
      stop("argument `split' not equal to number or rows in snp.matrix")
    if (length(unique(split[!is.na(split)]))!=2)
      stop("argument `split' cannot be coerced into a 2-level factor")
    split <- as.factor(split)
  }
  else {
    if ((class(snps2)!="snp.matrix") && (class(snps2)!="X.snp.matrix"))
      stop("argument `snps2' is not a snp.matrix or X.snp.matrix")
    if (ncol(snps)!=ncol(snps2))
      stop("two snp.matrix objects have different numbers of columns")
  }
  .Call("test_switch", snps, snps2, split, prior.df, PACKAGE="snpMatrix")
}

