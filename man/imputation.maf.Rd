\name{imputation.maf}
\alias{imputation.maf}
\alias{imputation.r2}
\title{Extract MAF and R-squared values from imputation rules}
\description{
These functions extract key characteristics of regression-based
imputation rules stored as an object of class
\code{"snp.reg.imputation"}. \code{imputation.maf} extracts the minor
allele frequencies of the imputed SNPs and \code{imputation.r2} extracts
the prediction \eqn{R^2}.
}
\usage{
imputation.maf(rules)
imputation.r2(rules)
}

\arguments{
  \item{rules}{An object of class \code{"snp.reg.imputation"}}
}
\value{
  A numeric vector containing the extracted values
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{snp.reg.imputation-class}},
  \code{\link{snp.imputation}}}
\examples{
# These functions are currently defined as
function (rules) sapply(rules, function(x) x$maf)
function (rules) sapply(rules, function(x) x$r2)
}
\keyword{manip}

