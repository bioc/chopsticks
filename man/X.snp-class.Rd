\name{X.snp-class}
\docType{class}
\alias{X.snp-class}
\alias{coerce,X.snp,character-method}
\alias{coerce,X.snp,numeric-method}
\alias{coerce,X.snp,genotype-method}
\alias{show,X.snp-method}

\title{Class "X.snp"}
\description{
  Compact representation of data concerning single nucleotide
  polymorphisms (SNPs) on the X chromosome
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("snp", ...,
    Female=...)} or 
  by subset selection from an object of class \code{"X.snp.matrix"}.
  Holds one row or column of an object of class \code{"X.snp.matrix"}
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{The genotype data coded as 0, 1, 2, or 3. For
      males are coded as homozygious females}
    \item{\code{Female}:}{A logical array giving the sex of the sample(s)}
  }
}
\section{Extends}{
Class \code{"snp"}, directly.
Class \code{"raw"}, by class \code{"snp"}.
Class \code{"vector"}, by class \code{"snp"}.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "X.snp", to = "character")}:
      map to codes "A/A", "A/B", "B/B", "A/Y", "B/Y", or "" } 
    \item{coerce}{\code{signature(from = "X.snp", to = "numeric")}:
      map to codes 0, 1, 2, or NA}
    \item{coerce}{\code{signature(from = "X.snp", to = "genotype")}: Yet
    to be implemented}
    \item{show}{\code{signature(object = "X.snp")}: shows character
      representation of the object}
  }
}
\references{\url{http://www-gene.cimr.cam.ac.uk/clayton}}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{X.snp.matrix-class}},
  \code{\link{snp.matrix-class}}, \code{\link{snp-class}}
}
\examples{
data(testdata)
s <- Xchromosome[,1]
class(s)
s
}
\keyword{classes}
