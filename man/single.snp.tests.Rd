\name{single.snp.tests}
\alias{single.snp.tests}
\title{1-df and 2-df tests for genetic associations with SNPs}
\description{
  This function carries out tests for association between phenotype and
  a series of single nucleotide polymorphisms (SNPs), within strata
  defined by a possibly confounding factor. SNPs are
  considered one at a time and both 1-df and 2-df tests are
  calculated. For a binary phenotype, the
  1-df test is the Cochran-Armitage test (or, when
  stratified, the Mantel-extension test).
}
\usage{
single.snp.tests(phenotype, stratum, data = sys.parent(), snp.data, subset, snp.subset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phenotype}{A vector containing the values of the phenotype}
  \item{stratum}{Optionally, a factor defining strata for the analysis}
  \item{data}{A dataframe containing the \code{phenotype} and
    \code{stratum} data. The row names of this are linked with the row
    names of the \code{snps} argument to establish correspondence of
    phenotype and genotype data. If this argument is not supplied,
    \code{phenotype} and \code{stratum} are evaluated in the calling
    environment and should be in the same order as rows of \code{snps}}
  \item{snp.data}{An object of class \code{"snp.matrix"} containing the SNP
    genotypes to be tested}
  \item{subset}{A vector or expression describing the subset of subjects
    to be used in teh analysis. This is evaluated in the same
    environment as the \code{phenotype} and \code{stratum} arguments}
  \item{snp.subset}{A vector describing the subset of SNPs to be
    considered. Default action is to test all SNPs.} 
}
\details{
  Formally, the test statistics are score tests for generalized linear
  models with canonical link. That is, they are inner products between
  genotype indicators and the deviations of phenotypes from their
  stratum means. Variances (and covariances) are those of the
  permutation distribution obtained by randomly permuting phenotype
  within stratum.

  The \code{subset} argument can either be a logical vector of length
  equal to the length of the vector of  phenotypes, an integer vector
  specifying positions in the \code{data} frame, or a character vector
  containing names of the selected rows in the \code{data}
  frame. Similarly, the \code{snp.subset} argument can be a logical,
  integer, or character vector.
}
\value{
  A dataframe, with columns
  \item{chi2.1df}{Cochran-Armitage type test for additive genetic component}
  \item{chi2.2df}{Chi-squared test for both additive and dominance
    components}
  \item{N}{The number of valid data points used}
}
\references{Clayton (2008) Testing for association on the X chromosome
  \emph{Biostatistics} (In press)}
\note{
  The behaviour of this function for objects of class
  \code{X.snp.matrix} is as described by Clayton (2008). Males are
  treated as homozygous females and corrected variance estimates are
  used.
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{snp.lhs.tests}}, \code{\link{snp.rhs.tests}}}
\examples{
data(testdata)
results <- single.snp.tests(cc, stratum=region, data=subject.data,
   snp.data=Autosomes, snp.subset=1:10)
summary(results)
# QQ plot - see \code{\link{qq.chisq}}
qq.chisq(results$chi2.1df)
qq.chisq(results$chi2.2df)
}
\keyword{htest}

