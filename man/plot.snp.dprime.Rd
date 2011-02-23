\name{plot.snp.dprime}
\alias{plot.snp.dprime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to draw the pairwise D' in a eps file }
\description{
  \code{plot.snp.dprime} takes a
  \code{\link[=snp.dprime-class]{snp.dprime}} object and
  draw an eps file to visualize the pairwise D', $r^2$ and LOD.
}
\usage{
\method{plot}{snp.dprime}(x, filename, scheme = "standard", do.notes = FALSE,
metric=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{\link[=snp.dprime-class]{snp.dprime}}}
  \item{filename}{The output file name, preferably ending with ".eps" (not enforced)}
  \item{scheme}{The colour scheme used. Valid values are "standard" for
    the Haploview default, and "rsq" for grayscale $r^2$. More may come later}
  \item{do.notes}{Boolean for whether to generate pdf annotation-related
    code}
  \item{metric}{An integer vector, detailing the chromosome position
    of the SNP, to drawa scaled metric of the location of the SNP. If
    NULL, no metric would be drawn}
  \item{...}{place holder}
}
\details{
  Annotation is a little used pdf features where certain part of a pdf
  file are hot spots where one can get pop-up balloons containing extra
  information, which doesn't appear in print. This is written to imitate
  the extra information one can get from right-clicking in Haploview's
  GUI.
}
\value{
  return nothing. Write a file as a result. And if \code{do.notes} is
  specified, Will also suggest user to execute \code{ps2pdf -dEPSCrop
    <filename>} to get a suitable pdf. 
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
\note{
  Unfortunately, there are two problems with annotations: only Acrobat Reader (out of all
  the pdf viewers, e.g. xpdf, kpdf, evince, various ghostscript based
  viewers) implements the feature, and a few thousand annotations can
  really make Acrobat Reader crawl.

  Also, Acrobat Reader has an implementation limit of 200 inches of the
  widest dimension of a document. This translates to 1200 snps in the
  current implementation of the drawing code, hence a warning is
  emitted that pdf written this way is not viewable by Acrobat
  Reader.(but viewable by xpdf, etc). A work around is possible based
  on LaTeX pdfpage, or eps can be included with scaling in another
  document, to stay inside 200 inches.

  The eps output is curently scaled to fit the long edge to A4
  Landscape. Users of gsview (up to version 4.9, a popular postcript
  file viewer on windows) should configure "A4" and "Rotate Media" under
  "Media" for correct display. Other useful options for running
  ghostscript are "-dEPSCrop" (for clipping again the bounding box) and
  "-dAutoRotatePages=/None" (to disable the automatic rotation based on
  text-run direction, i.e. not to rotate for SNP names running
  the right-way up).
  
  There is a Google Summer of code \url{http://code.google.com/soc/}
  2006 project to improve kpdf's annotation support.
  \url{http://wiki.kde.org/tiki-index.php?page=KDE\%20Google\%20SoC\%202006\%20ideas#id60851}
  I am involved.
}
\seealso{\code{\link{snp.dprime-class}}}
\examples{
data(testdata)
# As for ld.snp example ...
data(testdata)
ldinfo <- ld.snp(Autosomes, start=1, end=500, depth=50)
# Now plot to an eps file
plot.snp.dprime(ldinfo, filename="test.eps")
}
\keyword{htest}
\keyword{models}
\keyword{models}
\keyword{hplot}

