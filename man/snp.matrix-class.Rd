\name{snp.matrix-class}
\docType{class}
\alias{snp.matrix-class}
\alias{initialize,snp.matrix-method} 
\alias{[,snp.matrix-method}
\alias{[,snp.matrix,ANY,ANY,ANY-method}
\alias{rbind2,snp.matrix,snp.matrix-method}
\alias{cbind2,snp.matrix,snp.matrix-method}
\alias{coerce,snp.matrix,numeric-method}
\alias{coerce,snp.matrix,character-method}
\alias{coerce,matrix,snp.matrix-method}
%% Not having a duplicate coerce alias here as one should read
%% X.snp.matrix when the method is inquired.
%\alias{coerce,snp.matrix,X.snp.matrix-method}
\alias{is.na,snp.matrix-method}
\alias{cbind,snp.matrix-method}
\alias{rbind,snp.matrix-method}
\alias{show,snp.matrix-method}
\alias{summary,snp.matrix-method}
\alias{switch.alleles,snp.matrix,ANY-method}


\title{Class "snp.matrix"}
\description{ This class defines objects holding large arrays of single
  nucleotide polymorphism (SNP) genotypes generated using array
  technologies.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("snp.matrix", x)}
  where \code{x} is a matrix with storage mode \code{"raw"}.
	 Chips (usually corresponding to samples or
	 subjects) define rows of the matrix while polymorphisms (loci) define
	 columns. Rows and columns will usually have names which can be
	 used to link the data to further data concerning samples and SNPs
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"} and storage
      mode \code{raw}  
    Internally, missing data are coded 00 and SNP genotypes are coded
    01, 02 or 03. }
  }
}
\section{Extends}{
Class \code{"matrix"}, from data part.
Class \code{"structure"}, by class \code{"matrix"}.
Class \code{"array"}, by class \code{"matrix"}.
Class \code{"vector"}, by class "matrix", with explicit coerce.
Class \code{"vector"}, by class "matrix", with explicit coerce.
}
\section{Methods}{
  \describe{
    \item{[]}{\code{signature(x = "snp.matrix")}: subset
      operations}
    \item{cbind2}{\code{signature(x = "snp.matrix", y = "snp.matrix")}:
      S4 generic function to provide cbind() for two or
      more matrices together by column. Row names must match and column
      names must not coincide. If the matrices are of the derived class
      \code{\link{X.snp.matrix-class}}, the \code{Female} slot values must also
      agree}
    \item{coerce}{\code{signature(from = "snp.matrix", to = "numeric")}:
      map to codes 0, 1, 2, or NA}
    \item{coerce}{\code{signature(from = "snp.matrix", to =
	"character")}: map to codes "A/A", "A/B", "B/B", ""}
    \item{coerce}{\code{signature(from = "matrix", to = "snp.matrix")}:
      maps numeric matrix (coded 0, 1, 2 or NA) to a snp.matrix}
    \item{coerce}{\code{signature(from = "snp.matrix", to =
	"X.snp.matrix")}:
      maps a snp.matrix to an X.snp.matrix. Sex is inferred from the
      genotype data since males should not be heterozygous at any locus.
      After inferring sex, heterozygous calls for males are set to
      \code{NA}}
    \item{is.na}{\code{signature(x = "snp.matrix")}: returns a logical
      matrix indicating whether each element is \code{NA}}
    \item{rbind2}{\code{signature(x = "snp.matrix", y = "snp.matrix")}:
      S4 generic function to provide rbind() for two or
      more matrices by row. Column names must match and duplicated row
      names prompt warnings}
    \item{show}{\code{signature(object = "snp.matrix")}: shows the size
      of the matrix (since most objects will be too large to show in full)}
    \item{summary}{\code{signature(object = "snp.matrix")}: returns
      summaries of the data frames returned by
      \code{\link{row.summary}} and \code{\link{col.summary}}}
    \item{is.na}{\code{signature(x = "snp.matrix")}: returns a logical
      matrix of missing call indicators}
    \item{switch.alleles}{\code{signature(x = "snp.matrix", snps
	="ANY")}: Recode specified columns of of the matrix to reflect
      allele switches}

	 }
}
\references{\url{http://sourceforge.net/projects/outmodedbonsai/}}
\author{David Clayton \email{david.clayton@cimr.cam.ac.uk}}
\note{This class requires at least version 2.3 of R}
\seealso{\code{\link{snp-class}}, \code{\link{X.snp-class}},
  \code{\link{X.snp.matrix-class}}
}
\examples{
data(testdata)
summary(Autosomes)

# Just making it up - 3-10 will be made into NA during conversion
snps.class<-new("snp.matrix", matrix(1:10))
snps.class
if(!isS4(snps.class)) stop("constructor is not working")

pretend.X <- as(Autosomes, 'X.snp.matrix')
if(!isS4(pretend.X)) stop("coersion to derived class is not S4")
if(class(pretend.X) != 'X.snp.matrix') stop("coersion to derived class is not working")

pretend.A <- as(Xchromosome, 'snp.matrix')
if(!isS4(pretend.A)) stop("coersion to base class is not S4")
if(class(pretend.A) != 'snp.matrix') stop("coersion to base class is not working")

# display the first 10 snps of the first 10 samples
print(as(Autosomes[1:10,1:10], 'character'))

# convert the empty strings (no-calls) explicitly to "NC" before
# writing to an (anonymous and temporary) csv file
csvfile <- tempfile()
write.csv(file=csvfile, gsub ('^$', 'NC',
                              as(Autosomes[1:10,1:10], 'character')
                             ), quote=FALSE)
unlink(csvfile)
}
\keyword{classes}
