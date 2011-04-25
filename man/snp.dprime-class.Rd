\name{snp.dprime-class}
\alias{snp.dprime-class}
\alias{print.snp.dprime}
\alias{niceprint}
\docType{class}
\title{Class "snp.dprime" for Results of LD calculation}
\description{
  The \code{snp.dprime} class encapsulates results returned by
  \code{\link{ld.snp}} (--- routine to calculate D', $r^2$ and LOD of a
  \code{\link{snp.matrix-class}} object, given
  a range and a depth) and is based on a list of three
  named matrices.

  The lower right triangle of the snp.dprime object returned by \code{\link{ld.snp}}
  always consists zeros. This is deliberate. The associated plotting routine
  would not normally access those elements either.
}

\value{
  The \code{snp.dprime} class is a list of 3 named
  matrices \code{dprime}, \code{rsq2} or \code{r}, \code{lod}, and an attribute
  \code{snp.names} for the list of snps involved. (Note that if $x$ snps
  are involved, the row numbers of the 3 matrices are $(x-1)$).
  Only one of \code{r} or \code{rsq2} is present.

  \item{dprime}{D'}
  \item{rsq2}{$r^2$}
  \item{r}{signed $r^2$}
  \item{lod}{Log of Odd's}
  \item{attr(*, class)}{"snp.dprime"}
  \item{attr(*, snp.names)}{character vectors of the snp names involved}

  All the matrices are defined such that the ($n, m$)th entry is the
  pair-wise value between the ($n$)th snp and the $(n+m)$th snp. Hence
  the lower right triangles are always filled with zeros.

  Invalid values are represented by an out-of-range value - currently we
  use -1 for D', $r^2$ (both of which are between 0 and 1), and -2 for
  $r$ (valid values are between -1 and +1). lod is set to zero in
  most of these invalid cases. (lod can be any value so it is not
  indicative).

}
\section{Methods}{
  See \code{\link{plot.snp.dprime}}.
}
\source{
  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
  ~~ possibly secondary sources and usages ~~
}
\note{
  TODO: Need a subsetting operator.

  TODO: an assemble operator
}
\examples{
data(testdata)
snps20.20 <- Autosomes[11:20,11:20]
obj.snp.dprime <- ld.snp(snps20.20)
class(obj.snp.dprime)
summary(obj.snp.dprime)
\dontrun{
# The following isn't executable-as-is example, so these illustrations
# are commented out to stop R CMD check from complaining:

> d<- ld.snp(all, 3, 10, 15)
rows = 48, cols = 132
... Done
> d
$dprime
     [,1] [,2] [,3]
[1,]    1    1    1
[2,]    1    1    1
[3,]    1    1    1
[4,]    1    1    0
[5,]    1    0    0

$rsq2
          [,1]      [,2]      [,3]
[1,] 1.0000000 0.9323467 1.0000000
[2,] 0.9285714 1.0000000 0.1540670
[3,] 0.9357278 0.1854481 0.9357278
[4,] 0.1694915 1.0000000 0.0000000
[5,] 0.1694915 0.0000000 0.0000000

$lod
          [,1]      [,2]      [,3]
[1,] 16.793677 11.909686 16.407120
[2,] 10.625650 15.117962  2.042668
[3,] 12.589586  2.144780 12.589586
[4,]  2.706318 16.781859  0.000000
[5,]  2.706318  0.000000  0.000000

attr(,"class")
[1] "snp.dprime"
attr(,"snp.names")
[1] "dil118"  "dil119"  "dil5904" "dil121"  "dil5905" "dil5906"

}
}
\author{Hin-Tak Leung \email{htl10@users.sourceforge.net}}
\keyword{classes}
