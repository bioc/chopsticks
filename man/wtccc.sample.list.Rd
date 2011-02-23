\name{wtccc.sample.list}
\alias{wtccc.sample.list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{read the sample list from the header of the WTCCC signal
  file format}
\description{
  This is a convenience function for constructing the sample list from
  the header of a WTCCC signal file.
}
\usage{
wtccc.sample.list(infile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infile}{One of the signal files in a set of 23 (it is advisable
    to use the smaller ones such as number 22, although it shouldn't matter).}
}
\details{
  The header of a WTCCC signal file is like this:

\preformatted{
  AFFYID RSID pos  AlleleA AlleleB 12999A2_A 12999A2_B ...
}

The first 5 fields are discarded. There after, every other token
is retained, with the "\_A" or "\_B" part removed to give the sample list.
  
See also \code{\link{read.wtccc.signals}} for more details. 
}
\value{
  The value returned is a character vector contain the sample names
  or the plate-well names as appropriate.
}
\references{\url{http://www.wtccc.org.uk}}
\author{Hin-Tak Leung \email{htl10@users.sourceforge.net}}
\seealso{\code{\link{read.wtccc.signals}}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
\keyword{file}% __ONLY ONE__ keyword per line
