\name{ld.with}
\alias{ld.with}
\alias{ld.with,snp.matrix,character-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function to calculate the LD measures of specific SNPs against
other SNPs}
\description{
  This function calculates the LD measures ($r^2$, D', LOD) of specific
  SNPs against other SNPs.
}
\usage{
ld.with(data, snps, include.itself = as.logical(length(snps) - 1), signed.r = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{either a \code{\link{snp.dprime-class}} object or a
      \code{\link{snp.matrix-class}} object}
  \item{snps}{A list of snps, some of which are found in \code{data}}
  \item{include.itself}{Whether to include LD measures of SNPs against
    itself - it is FALSE for one SNP, since in that case, the
    result is known and trivial; but otherwise TRUE}
  \item{signed.r}{Logical, whether to output signed r or $r^2$}
}
\details{
  Not all combinations of the \code{include.itself} and \code{signed.r}
  make sense, nor fully operational.
}
\value{
  The returned value is somewhat similar to a
  \code{\link[=snp.dprime-class]{snp.dprime}}
  object, but not the same. It is a list of 3 named
  matrices \code{dprime}, \code{rsq2} (or \code{r} depending
  on the input), \code{lod}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Hin-Tak Leung \email{htl10@users.sourceforge.net}}
\note{
  The \code{ld.with} function started its life as an extractor
  function to take the output of \code{\link{ld.snp}}, a
  \code{\link{snp.dprime-class}} object, to rearrange it
  in a more convenient form to focus on the LD's against specific SNPs,
  but then evolved to take a \code{\link{snp.matrix-class}} object
  alternatively and perform the same task directly and more efficiently.  
}
\section{Warning}{
  Because this is really two functions rolled into one, depending on the
  class of \code{data}, not all combinations of the
  \code{include.itself} and \code{signed.r} make sense, nor fully operational.

  Also, the two versions have slightly different idea about invalid
  values, e.g. the LOD value for a SNPs against itself, or $r^2$ for two
  monomorphic snps (such as one against itself).

}
\seealso{\code{\link{ld.snp}}, \code{\link{snp.dprime-class}}}
\examples{
data(testdata)
snps10 <- Autosomes[1:10,1:10]
obj.snp.dprime <- ld.snp(snps10)

# result1 and result2 should be almost identical
# except where noted in the warning section above:
result1 <- ld.with(obj.snp.dprime, colnames(snps10))
result2 <- ld.with(snps10, colnames(snps10))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dplot}
\keyword{manip}
\keyword{models}
