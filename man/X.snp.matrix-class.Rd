\name{X.snp.matrix-class}
\docType{class}
\alias{X.snp.matrix-class}
\alias{initialize,X.snp.matrix-method}
\alias{[,X.snp.matrix-method}
\alias{[<-,X.snp.matrix,ANY,ANY,X.snp.matrix-method}
\alias{coerce,X.snp.matrix,character-method}
\alias{coerce,snp.matrix,X.snp.matrix-method}
\alias{show,X.snp.matrix-method}
\alias{summary,X.snp.matrix-method}

\title{Class "X.snp.matrix"}
\description{This class extends the  \code{\link{snp.matrix-class}} to
    deal with SNPs on the X chromosome.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("X.snp.matrix", x,
  Female)}.
Such objects have an additional \code{slot} to objects of class
\code{"snp.matrix"}
consisting of  a logical array of the same length as the number of
rows. This array indicates whether the sample corresponding to that row
came from a female (\code{TRUE}) or a male (\code{FALSE}). 
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"} and storage mode
      \code{"raw"}} 
    \item{\code{Female}:}{Object of class \code{"logical"} indicating
      sex of samples}
  }
}
\section{Extends}{
Class \code{"snp.matrix"}, directly, with explicit coerce.
Class \code{"matrix"}, by class \code{"snp.matrix"}.
Class \code{"structure"}, by class \code{"snp.matrix"}.
Class \code{"array"}, by class \code{"snp.matrix"}.
Class \code{"vector"}, by class "snp.matrix", with explicit coerce.
Class \code{"vector"}, by class "snp.matrix", with explicit coerce.
}
\section{Methods}{
  \describe{
    \item{[]}{\code{signature(x = "X.snp.matrix")}: subset
      operations. Currently rather slow owing to excessive copying}
    \item{[<-}{\code{signature(x = "X.snp.matrix")}: subset
      assignment operation to replace part of an object}
    \item{coerce}{\code{signature(from = "X.snp.matrix", to =
	"character")}:  map to codes 0, 1, 2, or NA}
    \item{coerce}{\code{signature(from = "snp.matrix", to =
	"X.snp.matrix")}:
      maps a snp.matrix to an X.snp.matrix. Sex is inferred from the
      genotype data since males should not be heterozygous at any locus.
      After inferring sex, heterozygous calls for males are set to
      \code{NA}}
     \item{show}{\code{signature(object = "X.snp.matrix")}: map to codes
      "A/A", "A/B", "B/B", "A/Y", "B/Y" or ""}
    \item{summary}{\code{signature(object = "X.snp.matrix")}: calculate
      call rates, allele frequencies, genotype frequencies,
      and chi-square tests for
      Hardy-Weinberg equilibrium. Genotype frequencies are calculated for
      males and females separately and Hardy-Weinberg equilibrium tests
      use only the female data. Allele frequencies are calculated using
      data from both males and females. Results are returned as a
      dataframe with
      column names \code{Calls}, \code{Call.rate}, \code{MAF}, \code{P.AA},
      \code{P.AB}, \code{P.BB},  \code{P.AY}, \code{P.BY}, and
      \code{z.HWE}
    }
  }
}
\references{\url{http://www-gene.cimr.cam.ac.uk/clayton}}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{\code{\link{X.snp-class}}, \code{\link{snp.matrix-class}},
    \code{\link{snp-class}}
}
\examples{
data(testdata)
summary(summary(Xchromosome))
}
\keyword{classes}
