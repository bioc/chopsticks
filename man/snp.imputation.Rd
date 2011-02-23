\name{snp.imputation}
\alias{snp.imputation}
\title{Calculate regression-based imputation rules}
\description{
  Given two set of SNPs typed in the same subjects, this function
  calculates regression equations which can be used to impute one set
  from the other in a subsequent sample.
}
\usage{
snp.imputation(X, Y, pos.X, pos.Y, phase=FALSE, try = 50, r2.stop = 0.9, max.X = 4)
}
\arguments{
  \item{X}{An object of class \code{"snpMatrix"} or
    \code{"X.snp.matrix"} containing observations
    of the SNPs to be used for imputation ("regressor SNPs")}
  \item{Y}{An object of same class as \code{X} containing observations
  of the SNPs to be imputed in a future sample ("target SNPs")}
  \item{pos.X}{The positions of the regressor SNPs}
  \item{pos.Y}{The positions of the target SNPs}
  \item{phase}{See "Details" below}
  \item{try}{The number of potential regressor SNPs to be
    considered in the stepwise regression procedure around each target
    SNP . The nearest \code{try} regressor SNPs to each target SNP
    will be considered}
  \item{r2.stop}{The value of \eqn{R^2} at which to stop entering new
    variables in the step-wise regression. If set to 1.0 or greater, the
    AIC will be used}
  \item{max.X}{The maximum number of regressor SNPs to be used for each
    imputation}
}
\details{
  The routine carries out a series of step-wise regression analyses in
  which each Y SNP is regressed on the nearest \code{try} regressor SNPs. If
  \code{phase} is \code{TRUE}, the regressions are calculated at the
  chromosome (haplotype) level, variances being simply  \eqn{p(1-p)} and
  covariances estimated 
  using the same algorithm used in \code{\link{ld.snp}}. Otherwise, the
  analysis is carried out at the diplotype level based on  
  conventional variance and covariance estimates using the
  \code{"all.obs"} missing value treatment (see \code{\link{cov}}).
}
\value{
  An object of class
  \code{"snp.reg.imputation"}.
}
\references{
  Chapman J.M., Cooper J.D., Todd J.A. and Clayton D.G. (2003)
  \emph{Human Heredity}, \bold{56}:18-31.
}
\note{The \code{phase=TRUE} option is not yet implemented} 
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{snp.reg.imputation-class}}, \code{\link{ld.snp}}}
\examples{
# Remove 5 SNPs from a datset and derive imputation rules for them
library(snpMatrix)
data(for.exercise)
sel <- c(20, 1000, 2000, 3000, 5000)
to.impute <- snps.10[,sel]
impute.from <- snps.10[,-sel]
pos.to <- snp.support$position[sel]
pos.fr <- snp.support$position[-sel]
imp <- snp.imputation(impute.from, to.impute, pos.fr, pos.to)
}
\keyword{models}
\keyword{regression}
