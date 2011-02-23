\name{pool}
\alias{pool}
\title{Pool test results from several studies or sub-studies}
\description{
  Given the same set of "score" tests carried out in several studies or
  in several different sub-samples within a study, this function pools
  the evidence by summation of the score statistics and score
  variances. Currently it combines tests produced by
  \code{\link{single.snp.tests}} but it will be extended to deal with
  further tests.
}
\usage{
pool(..., score = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{Objects holding the (extended) test results. Currently
    these must be of class
    \code{\link[=snp.tests.single.score-class]{snp.tests.single.score}}}
  \item{score}{Is extended score information to
  be returned in the output object?}
}
\details{
This function works by recursive calls to the generic function 
\code{pool2} which pools the results of two studies. 
}
\value{
  An object of same class as the input objects (optionally without the
  \code{.score}) extension. Tests are produced for the \emph{union} of
  SNPs tested in all the input objects. 
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{snp.tests.single.score-class}}}
\examples{
# An artificial example which simply doubles the size of a study
library(snpMatrix)
data(testdata)
sst <- single.snp.tests(snp.data=Autosomes, cc, data=subject.data,
         score=TRUE)
sst2 <- pool(sst, sst)
summary(sst2)
}
\keyword{htest}
