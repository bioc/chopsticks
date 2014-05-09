single.snp.tests <- function(phenotype, stratum, data=sys.parent(),
                             snp.data, subset, snp.subset) {
  m <- match.call()
  smiss <- missing(stratum)
  ssmiss <- missing(subset)
  if (smiss)
    stratum <- NULL
  if (!is(snp.data, "snp.matrix"))
    stop("snp.data argument must be of class snp.matrix")
  nr.snps = nrow(snp.data)
  nc.snps = ncol(snp.data)
  if (missing(data)) { # phenotype and stratum are in calling environment
    phenotype <- as.numeric(phenotype)
    complt <- !is.na(phenotype)
    if (length(phenotype)!=nr.snps)
      stop("incompatible length for phenotype vector")
    if (!smiss) {
      stratum <- as.integer(as.factor(stratum))
      if (length(stratum)!=nr.snps)
        stop("incompatible length for stratum vector")
      complt <- complt & !is.na(stratum)
    }
    if (!ssmiss) {
      if (is.logical(subset)) {
        if (length(subset)!=nr.snps)
          stop("incompatible length for subset vector")
        else
          inss <- subset
      }
      else if (is.integer(subset)) {
        if (length(subset)>nr.snps ||
            any(subset>nr.snps | duplicated(subset)))
          stop("illegal subset vector")
        inss <- rep(FALSE, nr.snps)
        inss[subset] <- TRUE
      }
      else {
        stop("illegal type for subset vector")
      }
      complt <- complt & inss
    }
  }
  else { # phenotype, stratum, and subset  are in data dataframe
    nm.snps <- rownames(snp.data)
    nm.data <- rownames(data)
    nr.data <- nrow(data)
    # which points to rows in the data matrix
    which <- match(nm.snps, nm.data, nomatch=NA)
    phenotype <- as.numeric(eval(m$phenotype, envir=data))[which]
    complt <- !is.na(phenotype) 
    if (!smiss) {
      stratum <- as.integer(as.factor(eval(m$stratum, envir=data)))[which]
      complt <- complt & !is.na(stratum)
    }
    if (!ssmiss) {
      subset <- eval(m$subset, envir=data)
      if (is.logical(subset)) {
        if (length(subset)!=nr.data)
          stop("incompatible length for subset vector")
        inss <- subset
      }
      else {
        if (is.character(subset)) {
          subset <- match(subset, nm.data, nomatch=0)
          if (any(subset==0))
            stop("unmatched name in subset vector")
        }
        else if (is.integer(subset)) {
          if (length(subset)>nr.data ||
              any(subset>nr.data | duplicated(subset)))
            stop("illegal subset vector")
        }
        else {
          stop("illegal type for subset vector")
        }
        inss <- rep(FALSE, nr.data)
        inss[subset] <- TRUE
      }
      inss <- inss[which]
      inss[is.na(inss)] <- FALSE
      complt <- complt & inss
    }
  }
  subset <- (1:nr.snps)[complt]
  if (missing(snp.subset))
    snp.subset <- NULL
  else if (is.character(snp.subset)) {
    snp.subset <- match(snp.subset, colnames(snp.data))
    if (any(is.na(snp.subset)))
      stop("unmatched snps in snp.subset argument")
  }
  else if (is.logical(snp.subset)){
    if (length(snp.subset)!=nc.snps)
      stop("snp.subset element(s) out of range")
    snp.subset <- (1:nc.snps)[snp.subset]
  }
  else if (is.integer(snp.subset)) {
    if (any(snp.subset<1 | snp.subset>nc.snps))
      stop("snp.subset element(s) out of range")
  }
  else 
    stop("illegal type for snp.subset")
  .Call("one_at_a_time", phenotype, stratum, snp.data, subset, snp.subset,
        PACKAGE="chopsticks")
}
