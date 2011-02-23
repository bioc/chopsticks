\name{snp.tests.glm-class}
\docType{class}
\alias{snp.tests.glm-class}
\alias{[,snp.tests.glm,ANY,missing,missing-method}
\alias{chi.squared,snp.tests.glm,missing-method}
\alias{deg.freedom,snp.tests.glm-method}
\alias{effect.sign,snp.tests.glm,logical-method}
\alias{pool2,snp.tests.glm,snp.tests.glm,logical-method}
\alias{p.value,snp.tests.glm,missing-method}
\alias{show,snp.tests.glm-method}
\alias{summary,snp.tests.glm-method}
\alias{switch.alleles,snp.tests.glm,character-method}
  
\title{Class "snp.tests.glm"}
\description{The class of objects created by \code{\link{snp.lhs.tests}}
  and  \code{\link{snp.rhs.tests}}
}
\section{Objects from the Class}{Objects are lists of
  test results. Each result is itself a list 
  with the following elements:
  \describe{
    \item{parameters}{The names of the parameters tested}
    \item{chi.squared}{The chi-squared test statistic}
    \item{df}{The degrees of freedom for the test}
    \item{U}{The score vector}
    \item{V}{The upper triangle of the variance-covariance matrix of the
      score, stored as a 1-D array}
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "snp.tests.glm", i = "ANY", j =
	"missing", drop = "missing")}: Subsetting operator}
    \item{chi.squared}{\code{signature(x = "snp.tests.glm", df =
	"missing")}: Extract chi-squared test values}
    \item{deg.freedom}{\code{signature(x = "snp.tests.glm")}: Extract
      degrees of freedom for tests}
    \item{effect.sign}{\code{signature(x = "snp.tests.glm", simplify =
	"logical")}: Extract signs of associations. If \code{simpify} is
    \code{TRUE} then a simple vector is returned if all tests are on 1df}
    \item{pool2}{\code{signature(x = "snp.tests.glm", y =
	"snp.tests.glm", score = "missing")}: Combine results from two
      sets of tests}
    \item{p.value}{\code{signature(x = "snp.tests.glm", df =
	"missing")}: Extract \emph{p}-values}
    \item{show}{\code{signature(object = "snp.tests.glm")}: Show method}
    \item{summary}{\code{signature(object = "snp.tests.glm")}: Summary
      method}
    \item{switch.alleles}{\code{signature(x = "snp.tests.glm", snps =
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
