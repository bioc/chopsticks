read.plink <- function(bed, bim, fam) {
  lb <- nchar(bed);
  ext <- substr(bed, lb-3, lb);
  if (ext==".bed") {
    stub <- substr(bed, 1, lb-4)
  } else {
    stub <- bed
    bed <- paste(bed, ".bed", sep="")
  }
  if (missing(bim))
    bim <- paste(stub, ".bim", sep="")
  if (missing(fam))
    fam <- paste(stub, ".fam", sep="")
  df.bim <- read.table(bim)
  snps <- as.character(df.bim[,2])
  cat("Number of SNPs: ", length(snps), "\n")
  if (any(duplicated(snps)))
    stop("Duplicated SNP name(s)")
  cat(head(snps), " ...\n")
  df.fam <- read.table(fam)
  ped <- as.character(df.fam[,1])
  id <- as.character(df.fam[,2])
  cat("Number of samples: ", length(id), "\n")
  if (any(duplicated(id))) {
    cat("Subject IDs not unique - concatenating with pedigree IDs\n")
    id <- paste(ped, id, sep=":")
    if (any(duplicated(id)))
      stop("Sample IDs are still not unique")
  }
  cat(head(id), "...\n")
  .Call("readbed", bed, id, snps, PACKAGE="snpMatrix") 
}
