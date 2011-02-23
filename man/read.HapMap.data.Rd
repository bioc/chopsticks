\name{read.HapMap.data}
\alias{read.HapMap.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function to import HapMap genotype data as snp.matrix}
\description{
  Given a URL for HapMap genotype data, \code{read.HapMap.data}, 
  download and convert the genotype data into a \code{snp.matrix} class
  object, and saving snp support infomation into an associated data.frame.
}
\usage{
read.HapMap.data(url, verbose=FALSE, save=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
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
  snp of the current batch are examined in two passes. The first time to
  see which bases are present, excluding "N". Deletions and Insertions
  ('D' and 'I') are currently treated as no-calls.

  If more than 2 bases are observed in the batch specified in the url,
  the routine aborts, but so far this possibility has not arisen in tests.
  If there is exactly two, then allele 1 and 2 are assigned in
  alphabetical order (dbSNPAlleles entries seems to be always in
  dictionary order, so the assignment made should agree with a shorten
  version of the dbSNPAlleles entry). Likewise, if only "A" or "T"
  is observed, then we know automatically it is the first (assigned as
  "A/.") or the last allele (assigned as "./T") of a hypothetical
  pair, without looking at the dbSNPAlleles entry. For other observed cases of 1 base, the
  routine goes further and look at the dnSNPAlleles
  entry and see if it begins with "-/X/" or ends with "/X", as a single
  base, and compare it with the single base observed to see if it should be
  allele 1 (same as the beginning, or different from the end) and
  allele 2 (same as the end, or different from the beginning).
  If no decision can be made for a particular snp entry, the routine
  aborts with an appropriate message.
  (for zero observed bases, assignment is "./.", and of course,
  all observed genotypes of that snp are therefore converted to the equivalent of NA)

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
  the assignment supplementary data to see they are compatible. (e.g. for an
  "A/C/T" entry, one data set may have "C" only and thus have assignment
  "A/C" and have all of it assigned Allele 2 homozygotes, whereas
  another data set contains both "C" and "T" and thus the
  first set needs to be modified before joining).

  A typical run, chromosome 1 for CEU, contains about ~400,000 snps and
  ~100 samples, and the snp.matrix object is about ~60MB (40 million
  bytes for snps plus overhead) and similar for the support data
  (i.e. ~ 2x), takes about 30 seconds, and at peak memory usage requires
  ~ 4x . The actual download is ~20MB, which is compressed from
  ~200MB.  
  
}
\seealso{\code{\link{snp.matrix-class}}}
\examples{
\dontrun{

## ** Please be aware that the HapMap project generates new builds from
## ** to time and the build number in the URL changes.
## ** The follow is valid as of Aug 2009:

> library(snpMatrix)
> testurl <- "http://ftp.hapmap.org/genotypes/latest/forward/non-redundant/genotypes_chr1_CEU_r27_nr.b36_fwd.txt.gz"
> result1 <- read.HapMap.data(testurl)
trying URL
> 'http://ftp.hapmap.org/genotypes/latest/forward/non-redundant/genotypes_chr1_CEU_r27_nr.b36_fwd.txt.gz'
Content type 'application/x-gzip' length 20059083 bytes (19.1 Mb)
opened URL
==================================================
downloaded 19.1 Mb
 
Reading 174 samples
...
EOF reached after 314024 snps
 ...conversion complete...
> sum1 <- col.summary(result1$snp.data)

> head(sum1[is.finite(sum1$z.HWE),], n=10)
           Calls Call.rate         MAF      P.AA       P.AB       P.BB        z.HWE
rs6650104    164 0.9425287 0.012195122 0.9756098 0.02439024 0.00000000   0.15810183
rs9629043     88 0.5057471 0.085227273 0.8295455 0.17045455 0.00000000   0.87399051
rs11510103   161 0.9252874 0.006211180 0.9937888 0.00000000 0.00621118 -12.68857754
rs11497407    89 0.5114943 0.005617978 0.0000000 0.01123596 0.98876404   0.05329933
rs12565286   160 0.9195402 0.028125000 0.0000000 0.05625000 0.94375000   0.36605143
rs11804171    83 0.4770115 0.030120482 0.0000000 0.06024096 0.93975904   0.28293272
rs2977670     85 0.4885057 0.058823529 0.8823529 0.11764706 0.00000000   0.57622153
rs2977656     90 0.5172414 0.005555556 0.9888889 0.01111111 0.00000000   0.05299907
rs12138618    89 0.5114943 0.050561798 0.0000000 0.10112360 0.89887640   0.50240136
rs3094315    163 0.9367816 0.153374233 0.7116564 0.26993865 0.01840491   0.50328457
## ** Please be aware that the HapMap project generates new builds from
## ** to time and the build number in the URL changes.
## ** The follow is valid as of Aug 2009:
 
## This URL is broken up into two to fit the width of
## the paper. There is no need in actual usage:

testurl2 <- paste("http://ftp.hapmap.org/genotypes/latest/",
                  "forward/non-redundant/genotypes_chr1_JPT_r27_nr.b36_fwd.txt.gz",
                  sep="")
> result2 <- read.HapMap.data(testurl2)
...
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
