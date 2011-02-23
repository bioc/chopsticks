\name{ld.snp}
\alias{ld.snp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to calculate pairwise D', r-squared}
\description{
  \code{ld.snp} takes an object of \code{snp.matrix} class and suitable
  range and depth and calculation the pairwise D', $r^2$, LOD and return
  the result as a \code{\link[=snp.dprime-class]{snp.dprime}} object.
}
\usage{
ld.snp(snpdata, depth = 100, start = 1, end = dim(snpdata)[2], signed.r=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snpdata}{An object of \code{snp.matrix} class with M samples of
    N snps}
  \item{depth}{The depth or lag of pair-wise calculation. Should be
      between 1 and N-1; default 100. Using 0 (an invalid value) is the
      same as picking the maximum}
    \item{start}{The index of the start of the range of interest. Should
      be between 1 and (N-1); default 1}
    \item{end}{The index of the end of the range of interest. Should
      be between 2 and N. default N.}
    \item{signed.r}{Boolean for whether to returned signed $r$ values
      instead of $r^2$}
}
\details{
  The cubic equation and quadratic equation solver code is borrowed from GSL (GNU Scientific Library).
}
\value{
  return a \code{\link[=snp.dprime-class]{snp.dprime}}
  object, which is a list of 3 named
  matrices \code{dprime}, \code{rsq2} (or \code{r} depending on the input), \code{lod}, and an attribute
  \code{snp.names} for the list of snps involved. (Note that if $x$ snps
  are involved, the row numbers of the 3 matrices are $(x-1)$).
  Only one of \code{rsq2} or \code{r} is present.

  \item{dprime}{D'}
  \item{rsq2}{$r^2$}
  \item{r}{signed $r$}
  \item{lod}{Log of Odd's}
  All the matrices are defined such that the ($n, m$)th entry is the
  pair-wise value between the ($n$)th snp and the $(n+m)$th snp. Hence
  the lower right triangles are always filled with zeros. (See example
  section for the actual layout)
  
  Invalid values are represented by an out-of-range value - currently we
  use -1 for D', $r^2$ (both of which are between 0 and 1), and -2 for
  $r$ (valid values are between -1 and +1). lod is set to zero in
  most of these invalid cases. (lod can be any value so it is not indicative).

}
\references{
    Clayton, D.G. and Leung, Hin-Tak (2007) An R package for analysis of
  whole-genome association studies. 
  \emph{Human Heredity} \bold{64}:45-51.\cr
  GSL (GNU Scientific Library) \url{http://www.gnu.org/software/gsl/}
}
\author{Hin-Tak Leung \email{htl10@users.sourceforge.net}}
\note{
  The output \code{\link[=snp.dprime-class]{snp.dprime}} object is suitable
  for input to \code{\link{plot.snp.dprime}} for drawing.

  The speed of ``ld.snp'' LD calculation, on a single-processor
  opteron 2.2GHz box:

  unsigned $r^2$, 13191 snps, depth 100   = 36.4 s (~ 1.3 mil pairs)

  signed  r    , 13191 snps, depth 100   = 40.94s (~ 1.3 mil pairs)

  signed  r    , 13191 snps, depth 1500  = 582s   (~ 18.5 mil pairs)

  For depth=1500, it uses 500MB just for the three matrices. So I
  actually cannot do the full depth at ~13,000; full depth should be
  under 50 minutes for 87 mil pairs, even in the signed-r version.

  The LD code can be ran outside of R - mainly for debugging:

  \preformatted{
    gcc -DWITHOUT_R -o /tmp/hello pairwise_linkage.c solve_cubic.c \
       solve_quadratic.c -lm
  }

  When used in this form, it takes 9 numbers:

  \preformatted{
  $/tmp/hello 4 0 0 0 30 0 0 0 23
  case 3               <- internal code for which cases it falls in
  root count 1         <- how many roots
  trying 1.000000
  p = 1.000000
  4      0       0       6.333333        0.000000        0.000000
  0      30      0       0.000000        25.333333       0.000000
  0      0       23      0.000000        0.000000        25.333333
  57 8 38.000000 38 38
  8 0 0 46 30, 38 38 76 76
  0.333333 0.000000 0.000000 0.666667
  d' = 1.000000 , r2 = 1.000000, lod= 22.482643
  }
}
\seealso{\code{\link{snp.dprime-class}},
  \code{\link{plot.snp.dprime}}, \code{\link{ld.with}}}
\examples{
# LD stats between 500 SNPs at a depth of 50
data(testdata)
ldinfo <- ld.snp(Autosomes, start=1, end=500, depth=50)
}
\keyword{dplot}
\keyword{htest}
\keyword{models}
