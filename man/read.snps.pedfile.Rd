\name{read.snps.pedfile}
\alias{read.snps.pedfile}
\title{Read genotype data from  a LINKAGE "pedfile"}
\description{
  This function reads data arranged as a LINKAGE "pedfile" with some
  restrictions and returns a list of three objects: a data frame
  containing the initial 6 fields giving pedigree structure, sex and
  disease status, a vector or a data frame containing snp assignment
  and possibly other snp information, and an object of class
  \code{"snp.matrix"} or \code{"X.snp.matrix"} containing the genotype data
}
\usage{
read.snps.pedfile(file, snp.names=NULL, assign=NULL, missing=NULL, X=FALSE, sep=".", low.mem = FALSE)
}
\arguments{
  \item{file}{The file name for the input pedfile}
  \item{snp.names}{A character vector giving the SNP names. If an
    accompanying map file or an info file is present, it will be read and the
    information used for the SNP names, and also the information
    merged with the result. If absent,
    the SNPs will be named numerically ("1", "2", ...)}
  \item{assign}{A list of named mappings for which letter maps to which
    Allele; planned for the future, not currently used
  }
  \item{missing}{Meant to be a single character giving the code recorded for alleles
    of missing genotypes ; not used in the current code}
  \item{X}{If \code{TRUE} the pedfile is assumed to describe loci on the
    X chromosome}
  \item{sep}{The character separating the family and member identifiers
    in the constructed row names; not used}
  \item{low.mem}{Switch over to input with a routine which requires less memory to run,
    but takes a little longer. This option also has the disadvantage that
    assignment of A/B genotype is somewhat non-deterministic and depends
    the listed order of samples.}
}
\details{
  Input variables are assumed to take the usual codes.
  Genotype should be coded as pairs of single character
  allele codes (which can be alphanumeric or numeric), from either 'A',
  'C', 'G', 'T' or '1', '2', '3', '4', with 'N', '-' and '0' denoting a
  missing; everything else is considered invalid and would invalidate
  the whole snp; also more than 2 alleles also cause the snp to be
  marked invalid.

  Row names of the output objects are constructed by
  concatenation of the pedigree and member identifiers, "Family",
  "Individual" joined by ".", e.g. "Family.Adams.Individual.0".

  It has been called to the authors' attention that there are LINKAGE
  ped files out there of case-control-study origin which encode the entire
  collection as one big family. This is wrong because in a case-control
  study evey sample is supposed to be unrelated, and the PED file should be
  a few thousand families, each with 1 member, instead of 1 family with
  a few thousand members.

  To fix such faulty ped files, run this perl snipplet below, which
  joins the family field with the member field with an underscore as the
  new family field, and put "1" in the member field:

  \preformatted{
cat input.ped | perl -n -e \
 '($a, $b, $c) = split /\s/,$_, 3; print $a, "_", $b, " 1 ", $c;' > output.ped
  }

}
\value{
  \item{snps}{The output \code{"snp.matrix"} or \code{"X.snp.matrix"}}
  \item{subject.support}{A data frame containing the first six fields of
  the pedfile}
}
\author{Hin-Tak Leung}
\seealso{\code{\link{snp.matrix-class}}, \code{\link{X.snp.matrix-class}},
  \code{\link{read.snps.long}}, \code{\link{read.HapMap.data}},
  \code{\link{read.pedfile.info}}, \code{\link{read.pedfile.map}}}
\keyword{IO}
\keyword{file}
