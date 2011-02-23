read.pedfile.info <- function(file) {
  # do exactly 3 columns, even if there are less or more
  result <- read.table(file, col.names=c("snp.names", "position", "chromosome"),
                       as.is=TRUE, fill=TRUE, flush=TRUE)
  rownames(result) <- result$snp.names
  result
}

read.snps.pedfile <- function(file, snp.names=NULL, assign=NULL, missing=NULL, X=FALSE, sep="." , low.mem=FALSE) {
  # If there is no input names, try to see if we can find the accompanying info file,
  # and if possible, load it for the snp names
  join.info <- FALSE
  if(is.null(snp.names)) {
    info.file <- sub('.ped$', '.info', file)  
    if (!(file.access(info.file,mode=4))) {
      # file.access() return 0 for success
      cat("Found accompanying info file, reading it first\n")
      snp.info <- read.pedfile.info(info.file)
      snp.names <- rownames(snp.info)
      join.info <- TRUE
    }
  }
  if (!is.null(snp.names)) {
    snp.names <- as.character(snp.names)
  }
  if (!is.null(missing)) {
    missing <- as.character(missing)
  }
  
  # if file starts with http:// or ftp://, download it and replace the input
  # with the downloaded file
  if (length(grep("^(ftp|http|file)://", file)) > 0) { 
  # mode = b is needed for windows
    saved.file <- tempfile()
    status <- download.file(file, destfile=saved.file, mode="wb")    
  # download.file() supposedly should throw error already, or return 0/1
    if ((status != 0) && (status != 1)) {
      stop("Download has gone wrong\n");
    }
    file <- saved.file
  }
  
  if (!low.mem) {
    result <- .Call("read_pedfile", file, snp.names, missing,
                    as.logical(X), as.character(sep), PACKAGE="snpMatrix")
    if (join.info) {
      snp.info[['assignment']] <- as.factor(result$snp.support)
      result$snp.support <- snp.info
    } else {
      result$snp.support <- as.factor(result$snp.support)
    }
  } else {
    if (is.null(missing)) {
      missing <- "0"
    }
    result <-  .Call("readped", file, as.character(snp.names), as.character(missing),
                     as.logical(X), as.character(sep), PACKAGE="snpMatrix")
    if (join.info) {
      result$snp.support <- snp.info
    }
  }

  result
}

