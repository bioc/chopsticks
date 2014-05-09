\name{snp.pre}
\alias{snp.pre}
\alias{snp.post}
\title{Pre- or post-multiply a snp.matrix object by a general matrix}
\description{
  These functions first standardize the input \code{snp.matrix} in the
  same way as does the function \code{\link{xxt}}. The standardized
  matrix is then either pre-multiplied (\code{snp.pre}) or
  post-multiplied (\code{snp.post}) by a general matrix. Allele
  frequencies for standardizing the input snp.matrix may be supplied
  but, otherwise, are calculated from the input snp.matrix
}
\usage{
snp.pre(snps, mat, frequency=NULL)
snp.post(snps, mat, frequency=NULL)
}
\arguments{
  \item{snps}{An object of class \code{"snp.matrix"} or \code{"X.snp.matrix"}}
  \item{mat}{A general (numeric) matrix}
  \item{frequency}{A numeric vector giving the allele (relative)
    frequencies to be used for standardizing the columns of \code{snps}.
    If \code{NULL}, allele frequencies will be calculated
    internally. Frequencies should refer to the second (\code{B}) allele
  }
}
\details{
The two matrices must be conformant, as with standard matrix
multiplication. The main use envisaged for these functions is the
calculation of factor loadings in principal component analyses of large
scale SNP data, and the application of these loadings to other
datasets. The use of externally supplied allele frequencies for
standardizing the input snp.matrix is required when applying loadings
calculated from one dataset to a different dataset
}
\value{
  The resulting matrix product
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{xxt}}}
\examples{
##--
##-- Calculate first two principal components and their loading, and verify
##--	
# Make a snp.matrix with a small number of rows
data(testdata)
small <- Autosomes[1:20,]
# Calculate the X.X-transpose matrix
xx <- xxt(small, correct.for.missing=FALSE)
# Calculate the first two principal components and corresponding eigenvalues
eigvv <- eigen(xx, symmetric=TRUE)
pc <- eigvv$vectors[,1:2]
ev <- eigvv$values[1:2]
# Calculate loadings for first two principal components
Dinv <- diag(1/sqrt(ev))
loadings <- snp.pre(small,  Dinv \%*\% t(pc))
# Now apply loadings back to recalculate the principal components
pc.again <- snp.post(small, t(loadings) \%*\% Dinv)
print(cbind(pc, pc.again))
}
\keyword{array}
\keyword{multivariate}

