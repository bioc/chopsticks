\name{snp.tests.single-class}
\docType{class}
\alias{snp.tests.single-class}
\alias{[,snp.tests.single,ANY,missing,missing-method}
\alias{chi.squared,snp.tests.single,numeric-method}
\alias{p.value,snp.tests.single,numeric-method}
\alias{sample.size,snp.tests.single-method}
\alias{show,snp.tests.single-method}
\alias{summary,snp.tests.single-method}

\alias{snp.tests.single.score-class}
\alias{[,snp.tests.single.score,ANY,missing,missing-method}
\alias{chi.squared,snp.tests.single.score,numeric-method}
\alias{effect.sign,snp.tests.single.score,missing-method}
\alias{pool2,snp.tests.single.score,snp.tests.single.score,logical-method}
\alias{p.value,snp.tests.single.score,numeric-method}
\alias{switch.alleles,snp.tests.single.score,ANY-method}


\title{Classes "snp.tests.single" and "snp.tests.single.score"}
\description{These are classes to hold the objects created by
\code{\link{single.snp.tests}} and provide methods for extracting key
elements. The class \code{"snp.tests.single.score"} extends class
\code{"snp.tests.single"}
to  include the score and score variance statistics in order to
provide methods for pooling results from several studies or parts of a study
}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("snp.tests.single", ...)} and
  \code{new("snp.tests.single.score", ...)} but, more usually, will be created
  by calls to \code{\link{single.snp.tests}}
}
\section{Slots}{
  \describe{
    \item{\code{snp.names}:}{The names of the SNPs tested, as they
      appear as column names in the original \code{snp.matrix}}
    \item{\code{chisq}:}{A two-column
      matrix holding the 1 and 2 df association tests}
    \item{\code{N}:}{The numbers of
      observations included in each test}
    \item{\code{N.r2}:}{For tests on imputed SNPs, the product of \code{N} and
	the imputation \eqn{r^2}. Otherwise a zero-length object}
    \item{\code{U}:}{(class \code{"snp.tests.single.score"}) Score statistics}
    \item{\code{V}:}{(class \code{"snp.tests.single.score"}) Score variances}
  }
}
\section{Methods}{
  \describe{
    \item{[]}{\code{signature(x =  "snp.tests.single", i = "ANY")}:
      Subsetting operator}
    \item{[]}{\code{signature(x =  "snp.tests.single.score", i = "ANY")}:
      Subsetting operator}
    \item{chi.squared}{\code{signature(x = "snp.tests.single", df =
	"numeric")}:
      Extract 1- and 2-df chi-squared test values}
    \item{chi.squared}{\code{signature(x = "snp.tests.single.score", df =
	"numeric")}:
      Extract 1- and 2-df chi-squared test values}
    \item{effect.sign}{\code{signature(x = "snp.tests.single.score", simplify =
	"missing")}: Extract signs of associations tested by the 1df tests}
    \item{p.value}{\code{signature(x = "snp.tests.single", df =
	"numeric")}:
      Evaluate 1- and 2-df test p-values}
    \item{p.value}{\code{signature(x = "snp.tests.single.score", df =
	"numeric")}:
      Evaluate 1- and 2-df test p-values}
    \item{show}{\code{signature(object = "snp.tests.single")}:
      List all tests and p-values}
    \item{sample.size}{\code{signature(object = "snp.tests.single")}:
      Extract sample sizes for tests}
    \item{summary}{\code{signature(object = "snp.tests.single")}:
      Summarize all tests and p-values}
    \item{pool2}{\code{signature(x = "snp.tests.single.score",
	y = "snp.tests.single.score", score = "logical")}: Combine two
      sets of test results. Used recursively by {\code{\link{pool}}}}
    \item{switch.alleles}{\code{signature(x = "snp.tests.single.score", snps =
	"ANY")}: Emulate, in the score vector and its (co)variances,
      the effect of switching of the alleles for the specified tests}
  }
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{
\code{\link{single.snp.tests}}, \code{\link{pool}}
}
\examples{
showClass("snp.tests.single")
showClass("snp.tests.single.score")
}
\keyword{classes}
