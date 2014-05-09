\name{read.pedfile.map}
\alias{read.pedfile.map}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function to read the accompanying map file of a LINKAGE ped file}
\description{
  This function read the accompanying map file of a LINKAGE ped file,
  for the SNP names, position and chromosome.
}
\usage{
read.pedfile.map(file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{A Plink map file}
}
\details{
  One such map file is the one accompanying the sample ped file of
  Haploview.
}
\value{
  A data frame with columns "snp.names", "position", "chromosome".
}
\references{See the documentation and description of ped files in
  Haploview (\url{http://www.broad.mit.edu/mpg/haploview/})}
\author{Hin-Tak Leung \email{htl10@users.sourceforge.net}}
\note{
  This is used internally by \code{\link{read.snps.pedfile}}
  to read an accompanying map file.
}
\seealso{\code{\link{read.snps.pedfile}}}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
\keyword{IO}
\keyword{file}

