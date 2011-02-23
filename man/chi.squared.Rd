\name{chi.squared}
\alias{chi.squared}
\alias{p.value}
\title{Extract test statistics and p-values}
\description{
Generic functions to extract chi-squared values and p-values from SNP
association test objects
}
\usage{
chi.squared(x, df)
p.value(x, df)
}
\arguments{
  \item{x}{Currently an object of class \code{"snp.tests.single"} or
     \code{"snp.tests.single.score"}}
  \item{df}{Either the numeric value 1 or 2}
}
\details{
  Currently these function are only available for objects output by
  \code{\link{single.snp.tests}}. 
  Eventually versions will also be provided for output objects produced by
  \code{\link{snp.lhs.tests}} and  \code{\link{snp.rhs.tests}}.
}
\value{
  A numeric vector containing the
  chi-squared test statistics or p-values. The output vector has a \code{names}
  attribute.
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{snp.lhs.tests}}, \code{\link{snp.rhs.tests}},
  	\code{\link{snp.tests.single-class}},
    \code{\link{snp.tests.single.score-class}}
}
\keyword{utilities}

