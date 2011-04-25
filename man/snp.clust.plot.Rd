\name{snp.clust.plot}
\alias{snp.clust.plot}
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
\author{Hin-Tak Leung \email{htl10@users.sourceforge.net}}
\seealso{
\code{\link{read.wtccc.signals}}, \code{\link{snp.matrix-class}}
}
\examples{
load(system.file("data/Genotypes.GenTrain1.RData",package="chopsticks"))
ab.signals <- read.wtccc.signals(system.file("extdata/example-new.txt",
package="chopsticks"),"rs1")
snp.clust.plot(ab.signals[['rs1']], GenTrain1[,'rs1'], title='test')
}
\keyword{hplot}

