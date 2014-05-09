read.snps.long <-function(files, sample.id=NULL, snp.id=NULL, female=NULL,
                         fields=c(sample=1, snp=2, genotype=3, confidence=4),
                         codes=c("0", "1", "2"), threshold=0.9, lower=TRUE,
                         sep=" ", comment="#",  skip=0,
                         simplify=c(FALSE, FALSE),
                         verbose=FALSE, every=1000) {
  if (!is.null(female) && any(is.na(female)))
    stop("female argument contains one or more NA's")
  if (!is.integer(fields))
    mode(fields) <- "integer"
  .Call("insnp_new", files, sample.id, snp.id, female, fields,
        codes, threshold, lower, sep, comment, as.integer(skip),
        simplify, verbose, as.integer(every), PACKAGE="chopsticks")
}

read.snps.long.old <- function(file, chip.id, snp.id, codes, female,
                           conf=1, threshold=0.9, drop=FALSE,
                           sorted=FALSE, progress=interactive()) {
  nch <- length(chip.id)
  if (length(unique(chip.id))!=nch)
    stop("\"chip.id\" error - identifiers must be unique")
  nsnp <- length(snp.id)
  if (length(unique(snp.id))!=nsnp)
    stop("\"snp.id\" error - identifiers must be unique")
  if (!missing(female)) {
    if (!is.logical(female))
      stop("indicator array for female sex must be of type \"logical\"")
    if (length(female)!=nch)
      stop("incorrect length of indicator array for female sex")
    ifX <- TRUE
  }
  else {
    ifX <- FALSE
    female <- logical(1)
  }
  if (missing(codes)) {
    codes <- c("0", "1", "2")
    if (ifX)
      codes <- c(codes, "0", "2")
  }
  else {
    if (ifX) {
      if (!is.character(codes) || length(codes)!=5)
        stop("\"codes\" argument should be a character vector of length 5")
    }
    else {
      if (!is.character(codes) || length(codes)!=3)
        stop("\"codes\" argument should be a character vector of length 3")
    }
  }
  size <- nch*nsnp
  if (sorted)
    res <- .C("insnp_long", file,
            as.integer(nch), as.character(sort(chip.id)),
            as.integer(nsnp), as.character(sort(snp.id)),
            as.integer(ifX), as.integer(female),
            as.character(codes),
            as.integer(conf), as.double(threshold),
            as.integer(progress),
            matrix=raw(size), integer(9), integer(1), PACKAGE="chopsticks")
  else
    res <- .C("insnp_long_unsorted", file,
            as.integer(nch), as.character(chip.id),
            as.integer(nsnp), as.character(snp.id),
            as.integer(ifX), as.integer(female),
            as.character(codes),
            as.integer(conf), as.double(threshold),
            as.integer(progress),
            matrix=raw(size), integer(9), integer(1), PACKAGE="chopsticks")
  error <- res[[14]]
  if (error == (-1)) 
    stop("error sorting input file")
  else if (error == (-2))
    stop("error opening input file")
  else if (error == (-3))
    warning("end-of-file reached before matrix filled")
  else if (error > 0)
    stop("read error on input file after ", error, " records")
  counts <- res[[13]]
  cat("\n", counts[4], " genotype calls were read from input\n", sep="")
  if (counts[7]>0) {
    skip.percent <- sprintf("%.2f", counts[7]/counts[4] * 100)
    cat(counts[7], " calls (", skip.percent, "%) of input read were skipped due to non-matching snp.id/chip.id\n", sep="")
  }
  if (counts[5]>0)
    cat(counts[5], "calls were repeats\n")
  if (counts[8]>0)
    cat(counts[8], "calls were \"no calls\"\n");
  if (counts[9]>0) {
    cat(counts[9], "calls were rejected due to ")
    if (conf>0)
      cat("low")
    else
      cat("high")
    cat(" confidence score\n");
  }
  
  if (counts[2]>0)
    cat(counts[2], "genotypes were not called\n")
  if (counts[6]>0)
    cat(counts[6], "genotypes had conflicting calls\n")
  if (counts[3]>0)
    cat(counts[3], "genotypes could not be found on input file\n")
  called.percent <- sprintf("%.2f",100*counts[1]/size)
  cat(counts[1], " genotypes were successfully called (", called.percent,
      "%)\n", sep="")
  dim(res$matrix) <- c(nch, nsnp)
  rownames(res$matrix) <- chip.id
  colnames(res$matrix) <- snp.id

  if (drop) {
    entries <- .C("empty", as.integer(nch), as.integer(nsnp),
                  as.raw(res$matrix), rows=logical(nch), cols=logical(nsnp),
                  PACKAGE="chopsticks")
    if (all(entries$rows, entries$cols)) {
      if (ifX) {
        new("X.snp.matrix", res$matrix, Female=female)
      }
      else {
        new("snp.matrix", res$matrix)
      }
    }
    else {
      nch1 <- sum(entries$rows)
      if (nch1==0)
        stop("No relevant genotype was in the input file")
      else if (nch1 < nch)
        warning("No data for ", nch - nch1, " chip.id entries") 
      nsnp1 <- sum(entries$cols)
      if (nsnp1 < nsnp)
        warning("No data for ", nsnp - nsnp1 , " snp.id entries")
      if (ifX) {
        new("X.snp.matrix", res$matrix[entries$rows, entries$cols],
            Female=female)
      }
      else {
        new("snp.matrix", res$matrix[entries$rows, entries$cols])
      }
    }
  }
  else {
    if (ifX) {
      new("X.snp.matrix", res$matrix, Female=female)
    }
    else {
      new("snp.matrix", res$matrix)
    }
  }
}

