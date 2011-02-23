\name{snp.rhs.tests}
\alias{snp.rhs.tests}
\title{ Score tests with SNP genotypes as independent variable}
\description{
This function fits a generalized linear model with phenotype as
dependent variable and, optionally, one or more potential confounders of
a phenotype-genotype association as independent variable. A series of
SNPs (or small groups of SNPs) are then tested for additional
association with phenotype. In order to protect against misspecification
of the variance function,  "robust" tests may be selected.
}
\usage{
snp.rhs.tests(formula, family = "binomial", link, weights, subset, data = parent.frame(),
   snp.data, tests=NULL, robust = FALSE,
   control=glm.test.control(maxit=20, epsilon=1.e-4, R2Max=0.98), allow.missing=0.01)
}
\arguments{
  \item{formula}{The base model formula, with phenotype as dependent variable}
  \item{family}{A string defining the generalized linear model
    family. This currently should (partially)  match one of
    \code{"binomial"}, \code{"Poisson"}, \code{"Gaussian"} or
      \code{"gamma"} (case-insensitive)}
  \item{link}{A string defining the link function for the GLM. This
  currently should (partially) match one of \code{"logit"},
  \code{"log"}, \code{"identity"} or \code{"inverse"}. The
  default action is to use the "canonical" link for the family selected}
  \item{data}{The dataframe in which the base model is to be fitted}
  \item{snp.data}{An object of class \code{"snp.matrix"} or
    \code{"X.snp.matrix"} containing the SNP data} 
  \item{tests}{Either a vector of column names or  numbers for the SNPs
    to be tested, or a list of short vectors defining groups of SNPs to be
    tested (again by name or number). The default action is to carry out
    \emph{all} single SNP tests, but
    \code{\link{single.snp.tests}} will often achieve the same
    result much faster}
  \item{weights}{"Prior" weights in the generalized linear model}
  \item{subset}{Array defining the subset of rows of  \code{data} to use}
  \item{robust}{If \code{TRUE}, robust tests will be carried out}
  \item{control}{An object giving parameters for the IRLS algorithm
    fitting of the base model and for the acceptable aliasing amongst
    new terms to be tested. See\ code{\link{glm.test.control}}}
  \item{allow.missing}{The maximum proportion of SNP genotype that can
    be missing before it becomes necessary to refit the base model}
}
\details{
  The tests used are asymptotic chi-squared tests based on the vector of
  first and second derivatives of the log-likelihood with respect to the
  parameters of the additional model. The "robust" form is a generalized
  score test in the sense discussed by Boos(1992). The "base" model is
  first fitted, and a score test is performed for addition of one or
  more SNP genotypes to the model. Homozygous SNP genotypes are coded 0
  or 2 and heterozygous genotypes are coded 1. For SNPs on the X
  chromosome, males are coded as homozygous females. For X SNPs, it will
  often be
  appropriate to include sex of subject in the base model (this is
  not done automatically).
  
  If a \code{data} argument is supplied, the \code{snp.data} and
  \code{data} objects are aligned by rowname. Otherwise all variables in
  the model formulae are assumed to be stored in the same order as the
  columns of the \code{snp.data} object.
}
\value{
  A data frame containing, for each SNP,
  \item{Chi.squared}{The value of the  chi-squared test statistic}
  \item{Df}{The corresponding degrees of freedom}
  \item{Df.residual}{The residual degrees of freedom for the base model;
  \emph{i.e.} the number of observations minus the number of parameters
  fitted}
  For the binomial family model, the base model can, in some circumstances,
  lead to perfect prediction of some observations (\emph{i.e.} fitted
  probabilities of 0 or 1). These observations are ignored in
  subsequent calculations; in particular they are not counted
  in the residual degrees of freedom. Similarly for Poisson means fitted
  exactly to zero.
}
\references{Boos, Dennis D. (1992) On generalized score tests. \emph{The
  American Statistician}, \strong{46}:327-333.}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\note{
  A factor (or
  several factors) may be included as arguments to the function
  \code{strata(...)} in the \code{formula}. This fits all
  interactions of the factors so included, but leads to faster
  computation  than fitting these in the normal way. Additionally, a
  \code{cluster(...)} call may be included in the base model
  formula. This identifies clusters of potentially correlated
  observations (e.g. for members of the same family); in this case, an
  appropriate robust estimate of the variance of the score test is used.
}
\seealso{\code{\link{single.snp.tests}}, \code{\link{snp.lhs.tests}},
  \code{\link{snp.matrix-class}}, \code{\link{X.snp.matrix-class}}}
\examples{
data(testdata)
slt3 <- snp.rhs.tests(cc~strata(region), family="binomial",
   data=subject.data, snp.data= Autosomes, tests=1:10)
print(slt3)
}
\keyword{htest}
