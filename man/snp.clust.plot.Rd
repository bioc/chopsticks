\name{snp.clust.plot}
\alias{snp.clust.plot}
\alias{example-new}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  function to do colorized cluster plots
}
\description{
  colorised cluster plots, typically for the purpose of
  verifying the validity of clustering
}
\usage{
snp.clust.plot(cluster, gtype, title = "test")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cluster}{one member from output of \code{read.wtccc.signals}
}
  \item{gtype}{one column of a \code{snp.matrix-class} object.
}
  \item{title}{title of the plot
}
}
\details{
See vignette for usage.
}
\value{None
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
\code{\link{read.wtccc.signals}},\code{\link{snp.matrix-class}}
}
\examples{
load(system.file("data/Genotypes.GenTrain1.Rda",package="snpMatrix"))
ab.signals <- read.wtccc.signals(system.file("data/example-new.txt",
package="snpMatrix"),"rs1")
snp.clust.plot(ab.signals[['rs1']], GenTrain1[,'rs1'], title='test')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
