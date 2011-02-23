\name{read.snps.chiamo}
\alias{read.snps.chiamo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read genotype data from the output of Chiamo}
\description{
  This function reads data from the raw output of Chiamo
}
\usage{
read.snps.chiamo(filename, sample.list, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{List of file names of output from Chiamo ; 
  the outcome is the concatenation from runs of Chiamo, e.g. on blocks of SNPs,
  which is often done for practical reasons}
  \item{sample.list}{A character vector giving the sample list}
  \item{threshold}{Cut-off for the posterior probability for a no-call}
}
\details{
  The raw output of Chiamo consists of the first 5 columns of
  \code{\link{read.wtccc.signals}}, followed by triplets of posterior
  probabilities of calling A-A, A-B, or B-B.

  The sample list can typically be obtained using
  \code{\link{wtccc.sample.list}}, from one of the (smaller) signal
  files, which are the inputs to Chiamo.
}
\value{
  The result is a list of two items:
  \item{snp.data}{The genotype data as a \code{\link{snp.matrix-class}} object.}
  \item{snp.support}{The information from the first 5 columns of \code{\link{read.wtccc.signals}}.}
}
\references{To obtain a copy of the Chiamo software please email
  Jonathan L. Marchini <marchini@stats.ox.ac.uk>.
}
\author{Hin-Tak Leung \email{htl10@users.sourceforge.net}}
\note{
}
\seealso{\code{\link{wtccc.sample.list}}, \code{\link{read.wtccc.signals}}}
\examples{
#
}
\keyword{IO}
\keyword{file}
