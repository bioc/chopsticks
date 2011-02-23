\name{pool2}
\alias{pool2}
\title{Pool results of tests from two independent datasets}
\description{
Generic function to pool results of tests from two independent
datasets. It is not designed to be called directly, but is called
recursively by \code{\link{pool}}
}
\usage{
pool2(x, y, score)
}
\arguments{
  \item{x, y}{Objects holding the (extended) test results. 
    These must be of class
    \code{\link[=snp.tests.single.score-class]{snp.tests.single.score}}
    or
    \code{\link[=snp.tests.glm-class]{snp.tests.glm}}}
  \item{score}{Is extended score information to
  be returned in the output object?}
}
\value{
  An object of same class as the input objects (optionally without the
  \code{.score}) extension. Tests are produced for the \emph{union} of
  SNPs tested in all the input objects. 
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{pool}},
  \code{\link{snp.tests.single.score-class}},
  \code{\link{snp.tests.glm-class}},
  \code{\link{single.snp.tests}},
  \code{\link{snp.lhs.tests}},
  \code{\link{snp.rhs.tests}}
  }
\keyword{htest}
