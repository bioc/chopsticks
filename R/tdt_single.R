tdt.snp <- function(ped, id, father, mother, affected,
                         data=sys.parent(), snp.data, rules=NULL,
                         snp.subset=NULL,
                         check.inheritance=TRUE, robust=FALSE,
                         score=FALSE) {
  if (!is.null(rules) || !check.inheritance) 
    robust <- TRUE
  mcall <- match.call()
  if (!is(snp.data, "snp.matrix"))
    stop("snp.data argument must be of class snp.matrix")
  nr.snps = nrow(snp.data)
  nc.snps = ncol(snp.data)
  if (missing(data)) { # ped data are in calling environment
    ped <- as.character(ped)
    nped <- length(id)
    id <- as.character(id)
    if (length(id)!=nped)
      stop("incompatible length for `id' and `ped' arguments")
    father <- as.character(father)
    if (length(father)!=nped)
      stop("incompatible length for `father' and `ped' arguments")
    mother <- as.character(father)
    if (length(mother)!=nped)
      stop("incompatible length for `mother' and `ped' arguments")
    affected <- as.logical(affected)
    if (length(affected)!=nped)
      stop("incompatible length for `affected' and `ped' arguments")
    subject.names <- id
  }
  else { # ped data are in data dataframe
    data <- as.data.frame(data)
    nped <- nrow(data)
    subject.names <- rownames(data)
    if (missing(ped))
      ped <- as.character(data[,1])
    else
      ped <- as.character(eval(mcall$ped, envir=data))
    
    if (missing(id))
      id <- as.character(data[,2])
    else
      id <- as.character(eval(mcall$id, envir=data))
    
    if (missing(father))
      father <- as.character(data[,3])
    else
      father <- as.character(eval(mcall$father, envir=data))
    
    if (missing(mother))
      mother <- as.character(data[,4])
    else
      mother <- as.character(eval(mcall$mother, envir=data))
    
    if (missing(affected))
      affected <- (data[,6]==2)
    else
      affected <- as.logical(eval(mcall$affected, envir=data))
  }

  # Treat subjects with "affected" missing as not affected
  
  affected[is.na(affected)] <- FALSE

  # Correspondence between ped data and snp data

  in.snp <- match(subject.names, rownames(snp.data))
  have.snps <- !is.na(in.snp)

  # Father and mother locations in ped file

  unique <- paste(ped, id, sep=":")
  f.unique <-  paste(ped, father, sep=":")
  fpos <- match(f.unique, unique)
  m.unique <-  paste(ped, mother, sep=":")
  mpos <- match(m.unique, unique)

  # Potentially complete trios
  
  trio <- have.snps & affected & (!is.na(fpos)) & (have.snps[fpos]) &
                                 (!is.na(mpos)) & (have.snps[mpos])
  ntrio <- sum(trio, na.rm=TRUE)
  if (ntrio==0) {
    cat("No potentially complete trios to analyse\n")
    return(NULL)
  }
  pd.snps <- in.snp[trio] # Proband rows in snp.matrix
  fr.snps <- in.snp[fpos[trio]] # Fathers' rows in snp.matrix
  mr.snps <- in.snp[mpos[trio]] # Mothers' rows in snp.matrix
  
  clust <-   as.integer(factor(ped[trio]))
  cord <- order(clust)
  cat("Analysing", ntrio, "potentially complete trios in", max(clust),
      "different pedigrees\n")

  # Calculate scores and score variances

  scores <- .Call("score_tdt", pd.snps[cord], fr.snps[cord], mr.snps[cord],
                  clust[cord], snp.data, rules, snp.subset,
                  check.inheritance, robust, PACKAGE="snpMatrix")
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

