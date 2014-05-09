\name{read.wtccc.signals}
\alias{read.wtccc.signals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ read normalized signals in the WTCCC signal
  file format }
\description{
  \code{read.wtccc.signals} takes a file and a list of snp ids (either
  Affymetrix ProbeSet IDs or rs numbers), and extract the entries
  into a form suitable for plotting and further analysis  
}
\usage{
read.wtccc.signals(file, snp.list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{\code{file} contains the signals. There is no need to gunzip.}
  \item{snp.list}{ A list of snp id's. Some Affymetrix
    SNPs don't have rsnumbers both rsnumbers and Affymetrix ProbeSet IDs
    are accepted}
}

\details{

  Do not specify both rs number and Affymetrix Probe Set ID in the input;
  one of them is enough.

  The signal file is formatted as follows, with the first 5 columns
  being the Affymetrix Probe Set ID, rs number, chromosome position,
  AlleleA and AlleleB. The rest of the header containing the sample id
  appended with "\_A" and "\_B". 

\preformatted{  
  AFFYID         RSID       pos   AlleleA AlleleB 12999A2_A 12999A2_B ...
  SNP_A-4295769  rs915677   14433758  C     T     0.318183  0.002809 
  SNP_A-1781681  rs9617528  14441016  A     G     1.540461  0.468571 
  SNP_A-1928576  rs11705026 14490036  G     T     0.179653  2.261650
} 

  The routine matches the input list against the first and the 2nd column.
  
(some early signal files, have the first "AFFYID" missing - this routine
  can cope with that also)
  
}
\value{
  The routine returns a list of named matrices, one for each input SNP
  (\code{NULL} if the SNP is not found); the row names are sample IDs
  and columns are "A", "B" signals.
}
\references{\url{http://www.wtccc.org.uk}}
\author{Hin-Tak Leung \email{htl10@users.sourceforge.net}}
\note{
  TODO: There is a built-in limit to the input line buffer (65535)
  which should be sufficient for 2000 samples and 30 characters each.
  May want to seek backwards, re-read and dynamically expand if
  the buffer is too small. 
}
\seealso{}
\examples{
\dontrun{
answer <-
  read.wtccc.signals("NBS_22_signals.txt.gz", c("SNP_A-4284341","rs4239845"))
> summary(answer)
              Length Class  Mode
SNP_A-4284341 2970   -none- numeric
rs4239845     2970   -none- numeric

> head(a$"SNP_A-4284341")
               A        B
12999A2 1.446261 0.831480
12999A3 1.500956 0.551987
12999A4 1.283652 0.722847
12999A5 1.549140 0.604957
12999A6 1.213645 0.966151
12999A8 1.439892 0.509547
>
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
\keyword{file}

