\name{snp.cor}
\alias{snp.cor}
\title{Correlations with columns of a snp.matrix}
\description{
This function calculates Pearson correlation coefficients  between columns of a
\code{snp.matrix} and columns of an ordinary matrix. The two matrices
must have the same number of rows. All valid pairs are used in the
computation of each correlation coefficient.
} 
\usage{
snp.cor(x, y)
}
\arguments{
  \item{x}{An \var{N} by \var{M} \code{snp.matrix}}
  \item{y}{An \var{N} by \var{P} general matrix}
}
\details{
  This can be used together with \code{\link{xxt}} and
  \code{\link[base]{eigen}} to calculate standardized loadings in the principal
  components
}
\value{
An  \var{M} by \var{P} matrix of correlation coefficients
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\note{
  This version cannot handle X chromosomes
}
\seealso{\code{\link{xxt}}}
\examples{
# make a snp.matrix with a small number of rows
data(testdata)
small <- Autosomes[1:100,]
# Calculate the X.X-transpose matrix
xx <- xxt(small, correct.for.missing=TRUE)
# Calculate the principal components
pc <- eigen(xx, symmetric=TRUE)$vectors
# Calculate the loadings in first 10 components,
# for example to plot against chromosome position
loadings <- snp.cor(small, pc[,1:10])
}
\keyword{array}
\keyword{multivariate}


