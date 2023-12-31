\name{snp.lhs.tests}
\alias{snp.lhs.tests}
\title{Score tests with SNP genotypes as dependent variable}
\description{
  Under the assumption of Hardy-Weinberg equilibrium, a SNP genotype is
  a binomial variate with two trials for an autosomal SNP or with one or
  two trials (depending on sex) for a SNP on the X chromosome.
  With each SNP in an input
  \code{"snp.matrix"} as dependent variable, this function first fits a
  "base" logistic regression model and then carries out a score test for
  the addition of further term(s). The Hardy-Weinberg
  assumption can be relaxed by use of a "robust" option.
}
\usage{
snp.lhs.tests(snp.data, base.formula, add.formula, subset, snp.subset,
              data = sys.parent(), robust = FALSE,
              control=glm.test.control(maxit=20, epsilon=1.e-4, R2Max=0.98))
}
\arguments{
  \item{snp.data}{The SNP data, as an object of class
    \code{"snp.matrix"} or \code{"X.snp.matrix"} }
  \item{base.formula}{A \code{formula} object describing the base model,
    with dependent variable omitted }
  \item{add.formula}{A \code{formula} object describing the additional
    terms to be tested, also with dependent variable omitted}
  \item{subset}{An array describing the subset of observations to be
    considered}
  \item{snp.subset}{An array describing the subset of SNPs to be
    considered. Default action is to test all SNPs.}
  \item{data}{The data frame in which \code{base.formula},
    \code{add.formula} and \code{subset} are to be evaluated}
  \item{robust}{If \code{TRUE}, a test which does not assume
    Hardy-Weinberg equilibrium will be used}
  \item{control}{An object giving parameters for the IRLS algorithm
    fitting of the base model and for the acceptable aliasing amongst
    new terms to be tested. See\ code{\link{glm.test.control}}}
}
\details{
  The tests used are asymptotic chi-squared tests based on the vector of
  first and second derivatives of the log-likelihood with respect to the
  parameters of the additional model. The "robust" form is a generalized
  score test in the sense discussed by Boos(1992).
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
  For the logistic model, the base model can, in some circumstances,
  lead to perfect prediction of some observations (\emph{i.e.} fitted
  probabilities of 0 or 1). These observations are ignored in
  subsequent calculations; in particular they are not counted
  in the residual degrees of freedom.
}
\references{Boos, Dennis D. (1992) On generalized score tests. \emph{The
  American Statistician}, \strong{46}:327-333.}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\note{
  A factor (or
  several factors) may be included as arguments to the function
  \code{strata(...)} in the \code{base.formula}. This fits all
  interactions of the factors so included, but leads to faster
  computation  than fitting these in the normal way. Additionally, a
  \code{cluster(...)} call may be included in the base model
  formula. This identifies clusters of potentially correlated
  observations (e.g. for members of the same family); in this case, an
  appropriate robust estimate of the variance of the score test is used.
}
\seealso{\code{\link{glm.test.control}},\code{\link{snp.rhs.tests}}
    \code{\link{single.snp.tests}}, \code{\link{snp.matrix-class}},
    \code{\link{X.snp.matrix-class}}}
\examples{
data(testdata)
library(survival)
slt1 <- snp.lhs.tests(Autosomes[,1:10], ~cc, ~region, data=subject.data)
print(slt1)
slt2 <- snp.lhs.tests(Autosomes[,1:10], ~strata(region), ~cc,
                      data=subject.data)
print(slt2)
}
\keyword{htest}
