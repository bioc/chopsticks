\name{read.snps.long.old}
\alias{read.snps.long.old}
\title{Read SNP input data in "long" format (old version) }
\description{
This function reads SNP genotype data and creates an object of class
\code{"snp.matrix"} or \code{"X.snp.matrix"}.
Input data are assumed to be arranged as one line per
SNP-call (without any headers). This function can read gzipped files.
}
\usage{
read.snps.long.old(file, chip.id, snp.id, codes, female,
                          conf = 1, threshold = 0.9, drop=FALSE,
                          sorted=FALSE, progress=interactive())
}
\arguments{
  \item{file}{Name of file containing the input data. Input files
    which have been compressed by the \code{gzip} utility are recognized}
  \item{chip.id}{Array of type \code{"character"} containing (unique)
    identifiers for the chips, samples, or subjects for which calls are
    to be read. Other samples in the input data will be ignored}
  \item{snp.id}{Array of type \code{"character"} containing (unique)
    identifiers of the SNPs for which data will be read. Again, further
    SNPs in the input data will be ignored}
  \item{codes}{For autosomal SNPs, an array of length 3 giving the codes
    for the three genotypes, in the order homozygous(AA), heterozygous(AB),
    homozygous(BB). For X SNPs, an additional two codes for the male
    genotypes (AY and BY) must be supplied. All other codes will be  treated
    as "no call". The default codes are \code{"0"},  \code{"1"},
    \code{"2"} [,\code{"0"},  \code{"2"}]}
  \item{female}{If the data to be read refer to SNPs on the X chromosome, this
    argument must be supplied and should indicate whether each row of
    data refers to a female (\code{TRUE}) or to a male
    (\code{FALSE}). The output object will then be of class
    \code{"X.snp.matrix"}.}
  \item{conf}{Confidence score. See details }
  \item{drop}{If \code{TRUE}, any rows or columns without genotype calls
    will be dropped from the output matrix. Otherwise the full matrix,
    with rows and columns
    defined by the \code{chip.id} and \code{snp.id} arguments, will
    be returned}
  \item{threshold}{Acceptance threshold for confidence score}
  \item{sorted}{Is input file already sorted into the correct order (see
    details)? }
  \item{progress}{If \code{TRUE}, progress will be reported to the
    standard output stream}
}
\details{
  Data are assumed to be input with one line per call, in free
  format:\cr\cr
  \emph{<chip-id>} \emph{<snp-id>} \emph{<code for genotype call>}
  [\emph{<confidence>}] ... \cr\cr
  Currently, any fields following the first three (or four) are
  ignored. If the argument \code{sorted} is \code{TRUE}, the file is
  assumed to be sorted
  with \emph{snp-id} as primary key and
  \emph{chip-id} as secondary key using the current locale. The rows and
  columns of the returned matrix will also be ordered in this manner. If 
  \code{sorted} is set to \code{FALSE}, then an algorithm which avoids
  this assumption is used. The rows and columns of the returned matrix
  will then be in the same order as the input \code{chip_id} and
  \code{snp_id} vectors. Calls in which both id fields match elements in the
  \code{chip.id} and
  \code{snp.id} arguments are read in, after (optionally) checking that
  the level of confidence achieves a given threshold.
  Confidence level checking is
  controlled by the \code{conf} argument. \code{conf=0} indicates that
  no confidence score is present and no checking is done. \code{conf>0}
  indicates that calls with scores  \emph{above} \code{threshold} are accepted,
  while \code{conf<0} indicates that  only calls with scores \emph{below}
  \code{threshold} should be accepted. 

  The routine is case-sensitive and it is important that the
  \emph{<chip-id>} and \emph{<snp-id>} match the cases of
  \code{chip.id} and \code{snp.id} exactly. 
}
\note{If more than one instance of any
  combination of \code{chip_id} element and \code{snp_id} element
  passes the confidence threshold, the called to be used is decided by
  the following rules:
  \enumerate{
    \item{1}{Any call trumps "no-call"}
    \item{2}{In the event of call conflict, "no-call" is returned}
  }

  Use of \code{sorted=TRUE} is usually discouraged since the alternative
  algorithm is safer and, usually, not appreciably slower. However, if
  the input file is to be read multiple times and there is a reasonably
  close correspondence between cells of the matrix to be returned and
  lines of the input file, the sorted option can be faster.

  This function has been replaced by the more flexible function
  \code{\link{read.snps.long}}.
}
\value{
  An object of class \code{snp.matrix}.
}
\references{\url{http://www-gene.cimr.cam.ac.uk/clayton}}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk} and Hin-Tak Leung}
\seealso{\code{\link{snp.matrix-class}}, \code{\link{X.snp.matrix-class}}}
\keyword{manip}
\keyword{IO}
\keyword{file}
\keyword{utilities}
