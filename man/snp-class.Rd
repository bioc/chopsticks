\name{snp-class}
\docType{class}
\alias{snp-class}
\alias{coerce,snp,character-method}
\alias{coerce,snp,numeric-method}
\alias{coerce,snp,genotype-method}
\alias{show,snp-method}

\title{Class "snp"}
\description{Compact representation of data concerning single nucleotide
  polymorphisms (SNPs)}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("snp", ...)} or
  by subset selection from an object of class \code{"snp.matrix"}.
  Holds one row or column of an object of class \code{"snp.matrix"}
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{The genotype data coded as 0, 1, 2, or 3}
  }
}
\section{Methods}{
  \describe{
   \item{coerce}{\code{signature(from = "snp", to =
	"character")}: map to codes "A/A", "A/B", "B/B", or ""}
    \item{coerce}{\code{signature(from = "snp", to = "numeric")}:
      map to numeric values 0, 1, 2 or, for uncertain assignments, to
      the posterior expectation of the 0, 1, 2 code}
    \item{coerce}{\code{signature(from = "snp", to = "genotype")}: maps
      a single SNP to an object of class \code{"genotype"}. See
      the \code{"genetics"} package}.
    \item{show}{\code{signature(object = "snp")}: shows character
      representation of the object}
    \item{is.na}{\code{signature(x = "snp")}: returns a logical
      vector of missing call indicators}
  }
}
\references{\url{http://www-gene.cimr.cam.ac.uk/clayton}}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{snp.matrix-class}},
  \code{\link{X.snp.matrix-class}}, \code{\link{X.snp-class}}
}
\examples{
## data(testdata)
## s <- autosomes[,1]
## class(s)
## s
}
\keyword{classes}
