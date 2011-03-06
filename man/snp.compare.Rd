\name{snp.compare}
\alias{snp.compare}
\alias{GenTrain1}
\alias{GenTrain2}
\alias{Genotypes.GenTrain1}
\alias{Genotypes.GenTrain2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to do compare two snp.matrix class object
}
\description{
Currently this function just counts how many genotypes differ, for each SNP,
and output named integer counts, and the signed call differences.
}
\usage{
snp.compare(obj1, obj2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj1}{a \code{\link{snp.matrix-class}} object
}
  \item{obj2}{another \code{\link{snp.matrix-class}} object
}
}
\details{
The two \code{\link{snp.matrix-class}} objects must have the same row
and column names; e.g. from the same raw data clustered with
two different calling algorithms.
}
\value{
named integer vector of the counts of differed genotypes, and the signed
call differences, for each SNP
}
\author{Hin-Tak Leung \email{htl10@users.sourceforge.net}}
\examples{
# See "clustering comparison vignette" for detailed description of these two data files:
load(system.file("data/Genotypes.GenTrain1.RData",package="chopsticks"))
load(system.file("data/Genotypes.GenTrain2.RData",package="chopsticks"))
a<- snp.compare(GenTrain1, GenTrain2)
hist(a$count,breaks=50,col='black')
}
\keyword{manip}

