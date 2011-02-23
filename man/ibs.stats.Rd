\name{ibs.stats}
\alias{ibs.stats}
\title{function to calculate the identity-by-state stats of a group
  of samples}
\description{
  Given a \link{snp.matrix-class} or a \link{X.snp.matrix-class} object
  with $N$ samples, calculates some statistics about the relatedness of
  every pair of samples within.
}
\usage{
ibs.stats(x)
}
\arguments{
  \item{x}{a \link{snp.matrix-class} or a \link{X.snp.matrix-class}
    object containing $N$ samples}
}
\details{
  No-calls are excluded from consideration here.
}
\value{
  A data.frame containing $N (N-1)/2$ rows, where the row names are the
  sample name pairs separated by a comma, and the columns are: 
  
  \item{Count}{count of identical calls, exclusing no-calls}
  \item{Fraction}{fraction of identical calls comparied to actual
    calls being made in both samples}
}
\author{Hin-Tak Leung \email{htl10@users.sourceforge.net}}
\note{ 
  This is mostly written to find mislabelled and/or duplicate samples.

  Illumina indexes their SNPs in alphabetical order so the
  mitochondria SNPs comes first - for most purpose it is undesirable
  to use these SNPs for IBS purposes.

  TODO: Worst-case S4 subsetting seems to make 2 copies of a large object,
  so one might want to subset before rbind(), etc; a future version
  of this routine may contain a built-in subsetting facility to work
  around that limitation.
}  
\section{Warning}{
  In some applications, it may be preferable to
  subset a (random) selection of SNPs first - the
  calculation time increases as $N (N-1) M /2$ . Typically for N = 800
  samples and M = 3000 SNPs, the calculation time is about 1 minute. A
  full GWA scan could take hours, and quite unnecessary for simple
  applications such as checking for duplicate or related samples.}
\examples{
data(testdata)
result <- ibs.stats(Autosomes[11:20,])
summary(result)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
