# Set control parameters 

glm.test.control <-
function(maxit, epsilon, R2Max) 
  list(maxit=as.integer(maxit),
       epsilon=as.double(epsilon),
       R2Max=as.double(R2Max))

# Logit tests with snps on LHS */

snp.lhs.tests <-
  
function(snp.data, base.formula, add.formula, subset, snp.subset,
         data=sys.parent(), robust=FALSE,
         control=glm.test.control(maxit=20, epsilon=1.e-4, R2Max=0.98))
{   
  if (is(snp.data,"snp.matrix")) 
    snames <- rownames(snp.data)
  else
    stop("snp.data must be stored as a snp.matrix object")
  
  call <- match.call()
  temp <- c("", "base.formula", "subset", "data")
  m <- call[match(temp, names(call), nomatch=0)]
  data.missing <- is.null(m$data)
  special <- c("strata", "cluster")
  Terms <- terms(base.formula, special, data=data)

  # Separate intercept term from rest of model formula
  # Note as of now I ignore this and fit an intercept anyway!

  INTERCEPT <- attr(Terms, "intercept")
  if (!INTERCEPT) # put one back temporarily
    attr(Terms, "intercept") <- 1

  # Model frame for base model

  names(m)[2] <- "formula"
  m[[1]] <- as.name("model.frame")
  m[[2]] <- Terms
  m <- as.call(c(as.list(m), list(na.action=as.symbol("na.omit"))))
  m <- eval(m)
  bnames <- rownames(m)

  # Base model matrix -- minus any cluster or strata variables

  strats <- attr(Terms, "specials")$strata
  clust <- attr(Terms, "specials")$cluster
  dropx <- NULL
  if (!is.null(clust)) {
    # Some check on valid clusters?
    if (missing(robust))
      robust <- TRUE
    tempc <- untangle.specials(Terms, "cluster", 1)
    clust <- as.integer(cluster(m[, tempc$vars]))  # unused 'shortlabel=TRUE' was here - possibly a mistake
    dropx <- tempc$terms
  }
  if (!is.null(strats)) {
    temps <- untangle.specials(Terms, "strata", 1)
    strats <- as.integer(strata(m[, temps$vars], shortlabel=TRUE))
    dropx <- c(dropx, temps$terms)
  }

  if (length(dropx)) 
    newTerms <- Terms[-dropx]
  else
    newTerms <- Terms
  X <- model.matrix(newTerms, m)

  # Remove unit column vector
  
  ones <- match("(Intercept)", colnames(X), nomatch=0)
  if (ones)
    X <- X[,-ones, drop=FALSE]
  
  # Additional model matrix
  
  temp <- c("", "add.formula", "subset", "data")
  m <- call[match(temp, names(call), nomatch=0)]
  Terms <- terms(add.formula, data=data)
     
  names(m)[2] <- "formula"
  m[[1]] <- as.name("model.frame")
  m <- as.call(c(as.list(m), list(na.action=as.symbol("na.omit"))))
  m <- eval(m)
  anames <- rownames(m) 
  Z <- model.matrix(Terms, m)

  # If there is an intercept term in the model already, remove this from Z

  ones <- match("(Intercept)", colnames(Z), nomatch=0)
  if (ones)
    Z <- Z[,-ones, drop=FALSE]
 
  
  # Find common rows and mapping to line up  snp.data

  common <- intersect(anames, bnames)
  if (!data.missing) {
    common <- intersect(common, snames)
    ix <- match(common, bnames)
    iz <- match(common, anames)
    iy <- match(common, rownames(snp.data))
  }
  else {
    ix <- match(common, bnames)
    iz <- match(common, anames)
    iy <- as.numeric(common) # common will refer to row numbers in snp.data  
  }
  N <- as.integer(length(common))
  if (N==0)
    stop("no observations")
  else {
    snp.data <- snp.data[iy,,drop=FALSE]
    X <- X[ix,, drop=FALSE]
    if (!is.null(strats))
      strats <- strats[ix]
    if (!is.null(clust))
      clust <- clust[ix]
    Z <- Z[iz,, drop=FALSE]
  }

  if (missing(snp.subset))
    snp.subset <- NULL
  else if (is.character(snp.subset)) {
    snp.subset <- match(snp.subset, colnames(snp.data))
    if (any(is.na(snp.subset)))
      stop("unmatched snps in snp.subset argument")
  }
  else if (is.logical(snp.subset)){
    if (length(snp.subset)!=ncol(snp.data))
      stop("snp.subset element(s) out of range")
    snp.subset <- (1:ncol(snp.data))[snp.subset]
  }
  else if (is.integer(snp.subset)) {
    if (any(snp.subset<1 | snp.subset>ncol(snp.data)))
      stop("snp.subset element(s) out of range")
  }
  else 
    stop("illegal type for snp.subset")

  .Call("snp_lhs_score",
        snp.data, X, strats, Z, snp.subset, robust, clust,
        control, PACKAGE="chopsticks")
}

# GLM tests with SNPs on RHS

snp.rhs.tests <-
function (formula, family="binomial", link, weights, subset,
          data=parent.frame(), snp.data,
          tests=NULL, robust=FALSE,
          control=glm.test.control(maxit=20, epsilon=1.e-4, R2Max=0.98),
          allow.missing=0.01) {
  
  call <- match.call()

  # Family and link function
  
  fam <- pmatch(tolower(family),
                c("binomial", "poisson", "gaussian", "gamma"))
  if (is.na(fam))
    stop("Unrecognized family argument")
  if (fam==0)
    stop("Ambiguous family argument")

  if (missing(link))
    lnk <- fam # Canonical link is default
  else {
    lnk <- pmatch(tolower(lnk), c("logit", "log", "identity", "inverse"))
    if (is.na(fam))
      stop("Unrecognized link argument")
    if (lnk==0)
      stop("Ambiguous link argument")
  }
  
  if (missing(data))
    data <- environment(formula)
  temp <- c("", "formula", "data", "weights", "subset")
  m <- call[match(temp, names(call), nomatch=0)]
  special <- c("strata", "cluster")
    
  Terms <- if (missing(data))
    terms(formula, special)
  else terms(formula, special, data=data)

  # Separate intercept term from rest of model formula

  INTERCEPT <- attr(Terms, "intercept")
  if (!INTERCEPT) # put one back temporarily
    attr(Terms, "intercept") <- 1

  # Model frame for null model

  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())

  # Find common rows and mapping to line up  snp.data

  if (is(snp.data, "snp.matrix")) 
    snames <- rownames(snp.data)
  else
    stop("snp.data must be stored as a snp.matrix")
  
  mnames <- rownames(m)
  map <- match(mnames, snames, nomatch=0)
  use <- map>0
  N <- sum(use)
  if (N==0)
    stop("no observations")
  

  # Response vector

  if(is.factor(m[[1]])){
    # This is just to avoid the warning that conversion
    # to numeric(=double) from factor is ignored.
    Y <- model.response(m, "any")
  } else {
    Y <- model.response(m, "numeric")
  }
  if (fam==1) { # Binomial family
    Y <- as.factor(Y)
    if (nlevels(Y)==2)
      Y <- as.numeric(Y)-1
    else
      stop("Response is not on two levels")
  }

  # Model matrix -- minus any cluster or strata variables

  strats <- attr(Terms, "specials")$strata
  clust <- attr(Terms, "specials")$cluster
  dropx <- NULL
  if (length(clust)) {
    # Some check on valid clusters?
    warning("cluster(...) not yet implemented")
    if (missing(robust))
      robust <- TRUE
    tempc <- untangle.specials(Terms, "cluster", 1)
    clust <- strata(m[, tempc$vars], shortlabel=TRUE)
    dropx <- tempc$terms
  }
  else
    clust <- NULL
  if (length(strats)) {
    temps <- untangle.specials(Terms, "strata", 1)
    strats <- strata(m[, temps$vars], shortlabel=TRUE)
    dropx <- c(dropx, temps$terms)
  }
  else
    strats <- NULL

  if (length(dropx)) 
    newTerms <- Terms[-dropx]
  else
    newTerms <- Terms
  X <- model.matrix(newTerms, m)

  # Remove unit column vector
  
  ones <- match("(Intercept)", colnames(X), nomatch=0)
  if (ones)
    X <- X[, -ones, drop=FALSE]
  if (ncol(X) == 0)
    X <- NULL
  
  # Prior weights

  if (!is.null(model.weights(m))) {
    cat("One\n")
    weights <- model.weights(m)
    if (any(weights<0))
      stop("negative prior weights not allowed")
  }
  else
    weights <- NULL

  # Offsets

  offset <- model.offset(m)
  if(!is.null(offset)) {
    offset <- offset
    stop("Can't deal with offsets at the moment")
  }

  # Sort and reshape arrays if necessary

  if (N != length(use)) {
    Y <- Y[use]
    if (!is.null(X))
      X <- X[use,,drop=FALSE]
    if (!is.null(weights))
      weights <- weights[use]
    if (!is.null(strats))
      strats <- strats[use]
    if (!is.null(clust))
      clust <- clust[use]
    snp.data <- snp.data[map,,drop=FALSE]
  }
  else if (any(map != 1:N))
    snp.data <- snp.data[map,,drop=FALSE]
    
    

  # Coerce tests argument to correct form #
  if(is.null(tests))
    tests <- as.integer(1:ncol(snp.data))

  # Private utility function to 
  # match character array  to list of names
  .col.numbers <- function(inc, names) {
    res <- match(inc, names, nomatch=NA)
    if (any(is.na(res))) {
      warning("Unmatched SNP names in tests argument")
      res <- tests[!is.na(res)]
    }
    res
  }
  
  if(is.character(tests)) {
    tests <- .col.numbers(tests, colnames(snp.data))
  } else {
    if(is.list(tests)) {
      first <- tests[[1]]
      if (is.character(first)) {
        tests <- lapply(tests, .col.numbers, colnames(snp.data))
      } else {
        if(!is.numeric(first))
          stop("illegal first element of tests argument")
      }
      tests<-lapply(tests,as.integer)
    } else {
      if(!(is.integer(tests) || is.null(tests)))
        stop("illegal tests argument")
    }
  }
  
  .Call("snp_rhs_score", Y, fam, lnk, X, strats, snp.data, weights,
        tests, robust, clust, control, allow.missing, PACKAGE="chopsticks")
}

