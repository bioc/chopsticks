\name{snp.tests.glm-class}
\docType{class}
\alias{snp.tests.glm-class}
\alias{[,snp.tests.glm,ANY,missing,missing-method}
\alias{chi.squared,snp.tests.glm,missing-method}
\alias{deg.freedom,snp.tests.glm-method}
\alias{effect.sign,snp.tests.glm,logical-method}
\alias{names,snp.tests.glm-method}
\alias{p.value,snp.tests.glm,missing-method}
\alias{sample.size,snp.tests.glm-method}
\alias{show,snp.tests.glm-method}
\alias{summary,snp.tests.glm-method}

\alias{snp.tests.glm.score-class}
\alias{[,snp.tests.glm.score,ANY,missing,missing-method}
\alias{pool2,snp.tests.glm.score,snp.tests.glm.score,logical-method}
\alias{switch.alleles,snp.tests.glm.score,character-method}
  
\title{Classes "snp.tests.glm" and "snp.tests.glm.score"}
\description{Classes of objects created by \code{\link{snp.lhs.tests}}
  and  \code{\link{snp.rhs.tests}}. The class
    \code{"snp.tests.glm.score"} extends the class
    \code{"snp.tests.glm"} and is invoked by setting the argument
    \code{score=TRUE} when calling testing functions
    in order to save the scores and their variances
    (and covariances)
}
\section{Objects from the Class}{Objects of class \code{"snp.tests.glm"}
  have four slots:
  \describe{
    \item{test.names}{A character vector of names}
    \item{chisq}{A numerical vector of chi-squared test values}
    \item{df}{An integer vector of degrees of freedom for the tests}
    \item{N}{A integer vector of the number of samples contributing to
      each test}
  }
  The \code{"snp.tests.glm.score"} class extends this, adding a slot
  \code{score} containing a list with elements which are themselves
  lists with two elements:
  \describe{
    \item{U}{The vector of efficient scores}
    \item{V}{The upper triangle of the variance-covariance matrix of
      \code{U}, stored as a vector}
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "snp.tests.glm", i = "ANY", j =
	"missing", drop = "missing")}: Subsetting operator}
    \item{chi.squared}{\code{signature(x = "snp.tests.glm", df =
	"missing")}: Extract chi-squared test values}
    \item{deg.freedom}{\code{signature(x = "snp.tests.glm")}: Extract
      degrees of freedom for tests}
    \item{names}{\code{signature(x="snp.tests.glm")}: Extract names of
      test values (\code{test.names} slot)}
    \item{p.value}{\code{signature(x = "snp.tests.glm", df =
	"missing")}: Extract \emph{p}-values}
    \item{sample.size}{\code{signature(object = "snp.tests.glm")}:
      Extract sample sizes for tests}
    \item{show}{\code{signature(object = "snp.tests.glm")}: Show method}
    \item{summary}{\code{signature(object = "snp.tests.glm")}: Summary
      method}

    \item{[}{\code{signature(x = "snp.tests.glm.score", i = "ANY", j =
	"missing", drop = "missing")}: Subsetting operator}
    \item{effect.sign}{\code{signature(x = "snp.tests.glm.score", simplify =
	"logical")}: Extract signs of associations. If \code{simpify} is
      \code{TRUE} then a simple vector is returned if all tests are on 1df}
    \item{pool2}{\code{signature(x = "snp.tests.glm.score", y =
	"snp.tests.glm.score", score = "missing")}: Combine results from two
      sets of tests}
    \item{switch.alleles}{\code{signature(x = "snp.tests.glm.score", snps =
	"character")}: Emulate, in the score vector and its (co)variances,
      the effect of switching of the alleles of specified SNPs}
  }
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\note{Most of the methods for this class are shared with
the \code{\linkS4class{snp.tests.single}} and
\code{\linkS4class{snp.tests.single.score}} classes
}
\seealso{\code{\link{snp.lhs.tests}},\code{\link{snp.rhs.tests}},
   \code{\linkS4class{snp.tests.single}},
   \code{\linkS4class{snp.tests.single.score}}
}
\examples{
showClass("snp.tests.glm")
}
\keyword{classes}
