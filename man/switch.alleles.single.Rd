\name{switch.alleles.single}
\alias{switch.alleles.single}
\title{Switch alleles in test results}
\description{
   In meta-analysis of two or more studies, tests are calculated
   separately for each study, saved, and then pooled across
   studies. However certain SNPs may have alleles reversed in one study,
   resulting in incorrect pooling of results. This function corrects the
   saved test statistics for such switches.
}
\usage{
switch.alleles.single(tests, switch)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tests}{An object of class \code{snp.tests.single.score}
    holding the test results}
  \item{switch}{A logical, integer, or character array identifying which
   SNPs should have alleles switched}
}
\details{
  It is sometimes more convenient to use this function rather than to
  use the function \code{switch.alleles} to correct the original data
  matrix. 
}
\value{
  An object of class \code{snp.tests.single.score} holding a copy of the
  input matrix in which the required allele switches have been made.
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{snp.tests.single.score-class}},
  \code{\link{switch.alleles}}
}
\keyword{utilities}
