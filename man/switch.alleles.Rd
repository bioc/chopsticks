\name{switch.alleles}
\alias{switch.alleles}
\title{Switch alleles in columns of a snp.matrix}
\description{
  Switch the allele coding for
  specified columns of an object of class \code{"snp.matrix"} or
  \code{"X.snp.matrix"}. Thus \emph{AA} genotypes become \emph{BB} and
  vice-versa while \emph{AB} remain unchanged.
}
\usage{
switch.alleles(snps, cols)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snps}{The input \code{snp.matrix} or \code{X.snp.matrix}}
  \item{cols}{A vector of type integer, character or logical specifying
  the column selection}
}
\value{
  A \code{snp.matrix} in which the specified columns have had their
  alleles switched.
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{test.allele.switch}},
  \code{\link{snp.matrix-class}},
  \code{\link{X.snp.matrix-class}}}
\examples{
library(snpMatrix)
data(testdata)
which <-  c("173774", "173811")
Asw <- switch.alleles(Autosomes, which)
summary(Autosomes)[which,]
summary(Asw)[which,]
}
\keyword{utilities}

