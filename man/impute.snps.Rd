\name{impute.snps}
\alias{impute.snps}
\title{Impute snps}
\description{
Given SNPs stored in an object of class \code{"chopsticks"} or
\code{"X.snp.matrix"} and a set of imputation equations in as object of
class \code{"snp.reg.imputation"}, this function calculates imputed values.
}
\usage{
impute.snps(rules, snps, subset=NULL)
}
\arguments{
  \item{rules}{The imputation equations; an object of class
    \code{"snp.reg.imputation"}}
  \item{snps}{The object of class \code{"chopsticks"} or
    \code{"X.snp.matrix"} containing the observed SNPs}
  \item{subset}{A vector describing the subset of subjects
    to be used. If \code{NULL} (default), then use all subjects}
}
\value{
  A matrix with imputed SNPs as columns. The imputed values are the
  estimated expected values of each SNP when coded 0, 1 or 2.
}
\references{
  Chapman J.M., Cooper J.D., Todd J.A. and Clayton D.G. (2003)
  \emph{Human Heredity}, \bold{56}:18-31.
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\note{
  Because the imputation is based on a linear model, the imputed value
  may lie outside the range 0 to 2.
}
\seealso{\code{\link{snp.imputation}}}
\examples{
# Remove 5 SNPs from a datset and derive imputation rules for them
library(chopsticks)
data(for.exercise)
sel <- c(20, 1000, 2000, 3000, 5000)
to.impute <- snps.10[,sel]
impute.from <- snps.10[,-sel]
pos.to <- snp.support$position[sel]
pos.fr <- snp.support$position[-sel]
imp <- snp.imputation(impute.from, to.impute, pos.fr, pos.to)
# Now calculate the imputed values
imputed <- impute.snps(imp, impute.from)
}
\keyword{models}
\keyword{regression}

