\name{switch.alleles}
\alias{switch.alleles}
\title{Switch alleles in columns of a snp.matrix or in test results}
\description{
  This is a generic function which can be applied to objects of class
  \code{"snp.matrix"} or  \code{"X.snp.matrix"} (which hold SNP genotype
  data), or to objects of class 
  \code{"snp.tests.single.score"} or \code{"snp.tests.glm"} (which hold
  association test results). In the former case, specified SNPs can be
  recoded as if the alleles were switched (so that \emph{AA} genotypes
  become \emph{BB} and 
  vice-versa while \emph{AB} remain unchanged). In the latter case,
  test results are modified \emph{as if} alleles had been switched. 
}
\usage{
switch.alleles(x, snps)
}
\arguments{
  \item{x}{The input object, of class  \code{"snp.matrix"},
    \code{"X.snp.matrix"}, \code{"snp.tests.single.score"}, or
    \code{"snp.tests.glm"}}
  \item{snps}{A vector of type integer, character or logical specifying
      the SNP to have its alleles switched}
}
\value{
  An object of the same class as the input object
}
\note{
  Switching alleles for SNPs has no effect on test results. These
  functions are required when carrying out meta-analysis, bringing
  together several sets of results. It is then important that alleles
  line up in the datasets to be combined. It is often more convenient
  (and faster) to apply this process to the test result objects rather
  than to the genotype data themselves.
  }
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{  \code{\link{snp.matrix-class}},
  \code{\link{X.snp.matrix-class}},
  \code{\link{snp.tests.single-class}},
  \code{\link{snp.tests.glm-class}}
}
\examples{
library(snpMatrix)
data(testdata)
which <-  c("173774", "173811")
Asw <- switch.alleles(Autosomes, which)
summary(Autosomes)[which,]
summary(Asw)[which,]
}
\keyword{utilities}

