\name{read.HapMap.data}
\alias{read.HapMap.data}
\title{ function to import HapMap genotype data as snp.matrix}
\description{
  Given a URL for HapMap genotype data, \code{read.HapMap.data},
  download and convert the genotype data into a \code{snp.matrix} class
  object, and saving snp support infomation into an associated data.frame.
}
\usage{
read.HapMap.data(url, verbose=FALSE, save=NULL, ...)
}
\arguments{
  \item{url}{URL for HapMap data. Web data is to be specified with
    prefix "http://", ftp data with prefix "ftp://", and local file as "file://"}
  \item{verbose}{Where the dnSNPalleles annotation is ambiguous, output
    more details information about how/why assignment is made. See
    Notes below.}
  \item{save}{filename to save the download - if unspecified, a
    temporary file will be created but removed afterwards.}
  \item{\dots}{Place-holder for further switches - currently
    ignored.}
}
\details{
  During the conversion, if the dbSNPAlleles entry is exactly of the
  form "X/Y", where X, Y = A or C or G or T, then it is used directly
  for assigning allele 1 and allele 2.

  However, about 1 in 1000 entries are more complicated e.g. may involving
  deletion, e.g. "-/A/G" or "-/A/AGT/G/T". Some heuristics are
  used in such cases, in which the observed genotypes in the specific
  snp of the current batch are examined in two passes.  The first time
  to see which bases are present, excluding "N".

  If more than 2 bases are observed in the batch specified in the url,
  the routine aborts, but so far this possibility has not arisen in tests.
  If there is exactly two, then allele 1 and 2 are assigned in
  alphabatical order (dbSNPAlleles entries seems to be always in
  dictionary order, so the assignment made should agree with a shorten
  version of the dbSNPAlleles entry). Likewise, if only "A" or "T"
  is observed, then we know automatically it is the first (assigned as
  "A/.") or the last allele (assigned as "./T") of a hypothetical
  pair, without looking at the dbSNPAlleles entry.  For other observed
  cases of 1 base, the
  routine goes further and look at the dnSNPAlleles
  entry and see if it begins with "-/X/" or ends with "/X", as a single
  base, and compare it with the single base observed to see if it should be
  allele 1 (same as the beginning, or different from the end) and
  allele 2 (same as the end, or different from the beginning).
  If no decision can be made for a particular snp entry, the routine
  aborts with an appropriate message.
  (For zero observed bases, assignment is "./.", and of course,
  all observed genotypes of that snp are therefore converted to the
  equivalent of NA.)

  (This heuristics does not cover all grounds, but practically it seems
  to work. See Notes below.)
}
\value{
  Returns a list containing these two items when successful, otherwise
  returns NULL:
  \item{snp.data}{A \code{snp.matrix-class} object containing the snp data}
  \item{snp.support}{A data.frame, containing the \code{dbSNPalleles},
        \code{Chromosome}, \code{Position}, \code{Strand} entries from
	the hapmap genotype file, together with the actual
	\code{Assignment} used for allele 1 and allele 2 during the
	conversion (See Details above and Note below).}
}
\references{\url{http://www.hapmap.org/genotypes}}
\author{Hin-Tak Leung \email{htl10@users.sourceforge.net}}
\note{

  Using both "file://" for \code{url} and \code{save} duplicates the file.
  (i.e. by default, the routine make a copy of the url in any case, but
  tidy up afterwards if run without \code{save}).

  Sometimes the assignment may not be unique e.g. dnSNPAlleles entry
  "A/C/T" and only "C" is observed - this can be assigned "A/C" or
  "C/T". (currently it does the former). One needs to be especially
  careful when joining two sets of snp data and it is imperative to compare
  the assigment supplementary data to see they are compatible. (e.g. for an
  "A/C/T" entry, one data set may have "C" only and thus have assignment
  "A/C" and have all of it assigned Allele 2 homozygotes, whereas
  another data set contains both "C" and "T" and thus the
  first set needs to be modified before joining).

  A typical run, chromosome 1 for CEU, contains about ~400,000 snps and
  ~100 samples, and the snp.matrix object is about ~60MB (40 million
  bytes for snps plus overhead) and similiar for the support data
  (i.e. ~ 2x), takes about 30 seconds, and at peak memory usage requires
  ~ 4x . The actual download is ~20MB, which is compressed from
  ~200MB.

}
\seealso{\code{\link{snp.matrix-class}}}
\examples{
\dontrun{

## ** Please be aware that the HapMap project generates new builds from
## ** time to time and the build number in the URL changes.

> library(snpMatrix)
> testurl <- paste0("http://www.hapmap.org/genotypes/latest/fwd_strand/",
                    "non-redundant/genotypes_chr1_CEU_r21_nr_fwd.txt.gz")
> result1 <- read.HapMap.data(testurl)
> sum1 <- summary(result1$snp.data)

> head(sum1[is.finite(sum1$z.HWE),], n=10)
           Calls Call.rate         MAF      P.AA       P.AB      P.BB      z.HWE
rs1933024     87 0.9666667 0.005747126 0.0000000 0.01149425 0.9885057 0.05391549
rs11497407    89 0.9888889 0.005617978 0.0000000 0.01123596 0.9887640 0.05329933
rs12565286    88 0.9777778 0.056818182 0.0000000 0.11363636 0.8863636 0.56511033
rs11804171    83 0.9222222 0.030120482 0.0000000 0.06024096 0.9397590 0.28293272
rs2977656     90 1.0000000 0.005555556 0.9888889 0.01111111 0.0000000 0.05299907
rs12138618    89 0.9888889 0.050561798 0.0000000 0.10112360 0.8988764 0.50240136
rs3094315     88 0.9777778 0.136363636 0.7272727 0.27272727 0.0000000 1.48118392
rs17160906    89 0.9888889 0.106741573 0.0000000 0.21348315 0.7865169 1.12733108
rs2519016     85 0.9444444 0.047058824 0.0000000 0.09411765 0.9058824 0.45528615
rs12562034    90 1.0000000 0.088888889 0.0000000 0.17777778 0.8222222 0.92554468

## ** Please be aware that the HapMap project generates new builds from
## ** to time and the build number in the URL changes.

## This URL is broken up into two to fit the width of
## the paper. There is no need in actual usage:
> testurl2 <- paste0("http://www.hapmap.org/genotypes/latest/",
         "fwd_strand/non-redundant/genotypes_chr1_JPT_r21_nr_fwd.txt.gz")
> result2 <- read.HapMap.data(testurl2)

> head(result2$snp.support)
           dbSNPalleles Assignment Chromosome Position Strand
rs10399749          C/T        C/T       chr1    45162      +
rs2949420           A/T        A/T       chr1    45257      +
rs4030303           A/G        A/G       chr1    72434      +
rs4030300           A/C        A/C       chr1    72515      +
rs3855952           A/G        A/G       chr1    77689      +
rs940550            C/T        C/T       chr1    78032      +
}
}
\keyword{manip}
\keyword{IO}
\keyword{file}
