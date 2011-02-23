single.snp.tests <- function(phenotype, stratum, data=sys.parent(), snp.data,
                             rules=NULL, subset, snp.subset, score=FALSE) {
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
  if (class(snp.data)=="X.snp.matrix" && any(is.na(snp.data@Female))) {
    warning("There are ", sum(is.na(snp.data@Female)),
            " subjects with NA for sex. These were ignored")
    complt <- complt & !is.na(snp.data@Female)
  }
  subset <- (1:nr.snps)[complt]


  if (missing(snp.subset))
    snp.subset <- NULL
  else {
    if (is.null(rules))
      n.targ <- nc.snps
    else
      n.targ <- length(rules)
    if (is.character(snp.subset)) {
      if (is.null(rules))
        snp.subset <- match(snp.subset, colnames(snp.data))
      else
        snp.subset <- match(snp.subset, names(rules))
      if (any(is.na(snp.subset)))
        stop("unmatched snps in snp.subset argument")
    }
    else if (is.logical(snp.subset)){
      if (length(snp.subset)!=n.targ)
        stop("snp.subset array has incorrect length")
      snp.subset <- (1:n.targ)[snp.subset]
    }
    else if (is.integer(snp.subset)) {
      if (any(snp.subset<1 | snp.subset>n.targ))
        stop("snp.subset element(s) out of range")
    }
    else 
      stop("illegal type for snp.subset")
  }
  scores <- .Call("score_single", phenotype, stratum, snp.data, rules, subset,
                  snp.subset, PACKAGE="snpMatrix")
  chisq <- .Call("chisq_single", scores, PACKAGE="snpMatrix")
  if (is.null(rules)) {
    if (is.null(snp.subset))
      tested <- colnames(snp.data)
    else
      tested <- colnames(snp.data)[snp.subset]
  } else {
    if (is.null(snp.subset))
      tested <- names(rules)
    else
      tested <- names(rules)[snp.subset]
  }
  if (score)
    res <- new("snp.tests.single.score", snp.names=tested,
               chisq=chisq, N=scores$N, N.r2=scores$N.r2,
               U=scores$U, V=scores$V)
  else
    res <- new("snp.tests.single", snp.names=tested, chisq=chisq,
               N=scores$N,  N.r2=scores$N.r2)
  res
}

pool <- function(..., score=FALSE) {
  argl <- list(...)
  na <- length(argl)
  while (na > 0 && is.null(argl[[na]])) {
    argl <- argl[-na]
    na <- na - 1
  }
  if (na<2)
    stop("need at least two test sets to pool")
  if (na>2) {
    p2 <- pool2(..1, ..2, score=TRUE)
    r <- do.call(pool, c(p2, argl[3:na], score=score))
  }
  else
    r <- pool2(..1, ..2, score=score)
  r
}

switch.alleles.single <- function(tests, switch) {
  if (class(tests)!="snp.tests.single.score") 
    stop("Object should be of class snp.tests.single.score")
  if (is.character(switch)) {
    switch <- match(switch, tests@snp.names)
    if (any(is.na(switch)))
      warning(sum(is.na(switch)), " SNP names were not found in tests object")
    switch <- switch[!is.na(switch)]
  } 
  ntest <- length(tests@snp.names)
  if (is.logical(switch)) {
    if (length(switch)!=ntest)
      stop("incompatible arguments")
    if (sum(switch)==0)
      return(tests)
  } else if (is.numeric(switch)) {
    if (length(switch)==0)
      return(tests)
    if (max(switch)>ntest || min(switch)<1)
      stop("incompatible arguments")
  } else {
    stop("SNPs to be switched must be indicated by name, position, or by a logical vector")
  }
  res <- tests
  res@U[switch,1] <- - tests@U[switch,1]
  if (ncol(tests@U)==3) {
    # X chromosome tests
    res@U[switch,2] <- - tests@U[switch,2]
    res@U[switch,3] <- tests@U[switch,3] - tests@U[switch,2]
    res@V[switch,4] <- tests@V[switch,4] - 2*tests@V[switch,3] +
                       tests@V[switch,2]
    res@V[switch,3] <- tests@V[switch,3] - tests@V[switch,2]
  } else {
    # Autosome tests
    res@U[switch,2] <- tests@U[switch,2] - tests@U[switch,1]
    res@V[switch,3] <- tests@V[switch,3] - 2*tests@V[switch,2] +
                       tests@V[switch,1]
    res@V[switch,2] <- tests@V[switch,2] - tests@V[switch,1]
  }
  res
}
