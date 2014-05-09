\name{row.summary}
\alias{row.summary}
\title{Summarize rows of a snp matrix}
\description{
This function calculates call rates and heterozygosity for each row of a an
object of class \code{"snp.matrix"}
}
\usage{
row.summary(object)
}
\arguments{
\item{object}{genotype data as a \code{\link{snp.matrix-class}} or \code{\link{X.snp.matrix-class}} object}
}
\value{
  A data frame with rows corresponding to rows of the input object
  and with columns/elements:
  \item{Call.rate }{Proportion of SNPs called}
  \item{Heterozygosity}{Proportion of called SNPs which are heterozygous}
}
\note{The current version does not deal with the X chromosome
  differently, so that males are counted as homozygous}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\examples{
data(testdata)
rs <- row.summary(Autosomes)
summary(rs)
rs <- row.summary(Xchromosome)
summary(rs)
}
\keyword{utilities}
