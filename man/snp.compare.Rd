\name{snp.compre}
\alias{snp.compare}
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
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
load(system.file("data/Genotypes.GenTrain1.Rda",package="snpMatrix"))
load(system.file("data/Genotypes.GenTrain2.Rda",package="snpMatrix"))
a<- snp.compare(GenTrain1, GenTrain2)
hist(a$count,breaks=50,col='black')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
