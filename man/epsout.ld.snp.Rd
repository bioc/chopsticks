\name{epsout.ld.snp}
\alias{epsout.ld.snp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to write an eps file directly to visualize LD }
\description{
  \code{epsout.ld.snp} takes an object of \code{snp.matrix} class and a given
  snp range and depth, draw a eps file to visualize the LD in the same color scheme
  as haploview's default view. It was the first prototype of
  this bunch of software. Also, it does not keep any pair-wise data
  in memory at all, and maybe more suitable where the actual pair-wise
  LD data is not needed.
}
\usage{
epsout.ld.snp(snpdata, filename, start, end, depth, do.notes=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snpdata}{An object of \code{snp.matrix} class with M samples of N snps}
  \item{filename}{The file name of the output, preferably ending with
    ".eps", but this rule not enforced}
  \item{start}{The index of the start of the range of interest. Should
    be between 1 and (N-1)}
  \item{end}{The index of the end of the range of interest. Should be
    between 2 and N.}
  \item{depth}{The depth or lag of pair-wise calculation. Should be
    between 1 and N-1}
  \item{do.notes}{Boolean for whether to generate pdf annotation-related
  code}
}
\details{
  The functionality of this routine has since been split into a two-stage
  processes involving \code{\link{ld.snp}} which generates a
  \code{\link[=snp.dprime-class]{snp.dprime}}
  object which contains the result of the
  pairwise LD calculation, and \code{\link{plot.snp.dprime}} (or the
    \code{plot} method of a \code{snp.dprime} object) which does the drawing.
}
\value{
  return nothing
}
\references{
    Clayton, D.G. and Leung, Hin-Tak (2007) An R package for analysis of
  whole-genome association studies.
  \emph{Human Heredity} \bold{64}:45-51.\cr
  GSL (GNU Scientific Library) \url{http://www.gnu.org/software/gsl/}\cr
  The postscript language reference manual:
  \url{http://www.adobe.com/products/postscript/pdfs/PLRM.pdf}\cr
  The pdf specification:
  \url{http://partners.adobe.com/public/developer/en/pdf/PDFReference16.pdf}
}
\author{Hin-Tak Leung \email{htl10@users.sourceforge.net}}
\seealso{\code{\link{snp.dprime-class}}, \code{\link{ld.snp}},
    \code{\link{plot.snp.dprime}}}
\examples{
#
data(testdata)
epsout.ld.snp(Autosomes, start=1, end=500, depth=50, filename="test.eps")
}
\keyword{dplot}
\keyword{hplot}
\keyword{models}
\keyword{htest}
