\name{snp.reg.imputation-class}
\docType{class}
\alias{snp.reg.imputation-class}
\alias{show,snp.reg.imputation-method}
\alias{summary,snp.reg.imputation-method}
\alias{plot,snp.reg.imputation,missing-method}
\alias{[,snp.reg.imputation,ANY,missing,missing-method}

\title{Class "snp.reg.imputation"}
\description{A class defining a list "rules" for imputation of
  SNPs. Rules are either linear regression equations or estimated
  haplotype probabilities for a target SNP and one or more predictor SNPs}
\section{Objects from the Class}{Objects are lists of
  \emph{rules}. Rules are named list elements each describing imputation
  of a SNP by a linear regression equation. Each element is itself a
  list with the following elements:
  \describe{
    \item{maf}{The minor allele frequency of the imputed SNP} 
    \item{r.squared}{The squared Pearson correlation coefficient between
    observed and predicted SNP duration derivation of the rule.}
    \item{snps}{The names of the SNPs to be included in the regression.}
    \item{coefficients}{A numeric array containing the regression
      equation 
      intercept followed by the regression coefficients for the SNPs
      listed in \code{snps}, OR}
    \item{hap.probs}{A numerical array containing estimated
      probabilities for haplotypes of the SNP to be imputed and all the
      predictor SNPs}
  }
  If any target SNP is monomorphic, the corresponding rule is
  returned as \code{NULL}. An object of class \code{snp.reg.imputation}
  has an attribute, \code{Max.predictors}, which gives the maximum
  number of predictors used for any imputation.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "snp.reg.imputation")}: prints
      an abreviated listing of the rules}
    \item{summary}{\code{signature(object = "snp.reg.imputation")}:
      returns a table which shows the distribution of r-squared values
      achieved against the number of snps used for imputation}
    \item{plot}{\code{signature(x="snp.reg.imputation", y="missing")}:
      plots the distribution of r-squared values as a stacked bar chart}
    \item{[]}{\code{signature(x = "snp.reg.imputation", i = "ANY")}: subset
      operations}
  }
}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\seealso{
  \code{\link{snp.imputation}}, \code{\link{impute.snps}},
  \code{\link{single.snp.tests}}
}
\examples{
showClass("snp.reg.imputation")
}
\keyword{classes}
