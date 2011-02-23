\name{pair.result.ld.snp}
\alias{pair.result.ld.snp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to calculate the pairwise  D', r-squared, LOD of
  a pair of specified SNPs}
\description{
  \code{pair.result.ld.snp.Rd} calculates the pairwise D', $r^2$, LOD of
    a pair of specified SNPs in a \code{snp.matrix} object. This is used
    mainly for debugging.
}
\usage{
pair.result.ld.snp(snpdata, loc.snpA, loc.snpB)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snpdata}{An object of \code{snp.matrix} class with M samples of
      N snps}
  \item{loc.snpA}{index of the first snp; should be between 1 and N}
  \item{loc.snpB}{index of the second snp; should be between 1 and N}
}
\value{
  Returns nothing. Results are displayed in stdout/console.  
}
\references{
    Clayton, D.G. and Leung, Hin-Tak (2007) An R package for analysis of
  whole-genome association studies. 
  \emph{Human Heredity} \bold{64}:45-51.\cr
  GSL (GNU Scientific Library) \url{http://www.gnu.org/software/gsl/}
}
\author{Hin-Tak Leung \email{htl10@users.sourceforge.net}}
\note{
  Not really recommended for daily usage; the result isn't saved
  anywhere and this routine is primarily for debugging the details
  and correctness of the calculation.
}
\seealso{\code{\link{snp.matrix-class}}}
\examples{
data(testdata)
pair.result.ld.snp(Autosomes, 1, 2)
}
\keyword{htest}
\keyword{models}
\keyword{dplot}
