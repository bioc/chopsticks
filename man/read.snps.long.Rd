\name{read.snps.long}
\alias{read.snps.long}
\title{Read SNP data in long format}
\description{
  Reads SNP data when organized in free format
  as one call per line. Other than the one
  call per line requirement, there is considerable flexibility. Multiple
  input files can be read, the input fields can be in any order on the
  line, and irrelevant fields can be skipped. The samples and SNPs
  to be read  must be pre-specified, and define rows and columns of an
  output object of class \code{"snp.matrix"}.  
}
\usage{
read.snps.long(files, sample.id = NULL, snp.id = NULL, female = NULL,
              fields = c(sample = 1, snp = 2, genotype = 3, confidence = 4),
              codes = c("0", "1", "2"), threshold = 0.9, lower = TRUE,
              sep = " ", comment = "#", skip = 0, simplify = c(FALSE,FALSE),
              verbose = FALSE, every = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{A character vector giving the names of the input files}
  \item{sample.id}{A character vector giving the identifiers of the
    samples to be read}
  \item{snp.id}{A character vector giving the names of the SNPs to be read}
  \item{female}{If the SNPs are on the X chromosome and the data are to
    be read as such, this logical vector (of the same length as
    \code{sample.id} should specify whether each sample was from a
    female subject}
  \item{fields}{A integer vector with named elements specifying the
    positions of the required fields in the input record. The fields are
    identified by the names \code{sample} and \code{snp} for the sample
    and SNP identifier fields, \code{confidence} for a call confidence
    score (if present) and either \code{genotype} if genotype calls
    occur as a single field, or \code{allele1} and \code{allele2} if the
    two alleles are coded in different fields}
  \item{codes}{Either the single string \code{"nucleotide"} denoting
    that coding in terms of nucleotides
    (\code{A}, \code{C}, \code{G} or \code{T}, case insensitive),
    or a character vector
    giving genotype or allele codes (see below)}
  \item{threshold}{A numerical value for the calling threshold on the
    confidence score}
  \item{lower}{If \code{TRUE}, then \code{threshold} represents a lower
    bound. Otherwise it is an upper bound}
  \item{sep}{The delimiting character separating fields in the input record}
  \item{comment}{A character denoting that any remaining input on a line
    is to be ignored}
  \item{skip}{An integer value specifying how many lines are to be
    skipped at the beginning of each data file}
  \item{simplify}{If \code{TRUE}, sample and SNP identifying strings
    will be shortened by removal of any common leading or trailing
    sequences when they are used as row and column names of the output
    \code{snp.matrix} }
  \item{verbose}{If \code{TRUE}, a progress report is generated as
    every \code{every} lines of data are read}
  \item{every}{See \code{verbose}}
}
\details{
  If nucleotide coding is not used, the \code{codes} argument
  should be a character array giving the valid codes. 
  For genotype coding of autosomal SNPs, this should be
  an array of length 3 giving the codes
  for the three genotypes, in the order homozygous(AA), heterozygous(AB),
  homozygous(BB). All other codes will be  treated
  as "no call". The default codes are \code{"0"},  \code{"1"},
  \code{"2"}.  For X SNPs, males are assumed to be coded as homozygous,
  unless an additional two codes are supplied (representing the 
  AY and BY genotypes). For allele coding, the
  \code{codes} array should be of length 2 and should specify the codes
  for the two alleles. Again, any other code is treated as
  "missing" and, for X SNPs, males should be coded either as
  homozygous or by omission of  the second allele.

  Although the function allows for reading of data for the X chromosome
  directly into an object of class \code{"X.snp.matrix"},
  it will often be preferable to read such data as a \code{"snp.matrix"}
  (i.e. as autosomal) and to coerce it to an object of type
  \code{"X.snp.matrix"} later using \code{as(..., "X.snp.matrix")} or
  \code{new("X.snp.matrix", ..., female=...)}.

  The vectors \code{sample.id} and \code{snp.id} must be in the same
  order as they  vary on the input file(s) and this ordering must be
  consistent. However, there is
  no requirement that either SNP or sample should vary fastest; this is
  detected from the input. 
  Each file may represent a separate sample or SNP, in which case the
  appropriate \code{.id} argument can be omitted and row or column names
  taken from the file names.  
}
\value{
  An object of class \code{"snp.matrix"} or \code{"X.snp.matrix"}.
}
\note{
  The function will read gzipped files.
  
  This function has replaced and earlier version which was much less
  flexible. Because all features have not been fully tested, the older
  version has been retained as \code{\link{read.snps.long.old}}.
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{read.HapMap.data}}\code{\link{read.snps.pedfile}},
  \code{\link{read.snps.chiamo}},\code{\link{read.snps.long}},
  \code{\link{snp.matrix-class}}, \code{\link{X.snp.matrix-class}}}
\keyword{manip}
\keyword{IO}
\keyword{file}
\keyword{utilities}
