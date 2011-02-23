\name{snpMatrix-package}
\alias{snpMatrix-package}
\alias{snpMatrix}
\docType{package}
\title{
The snp.matrix and X.snp.matrix classes
}
\description{
Implements classes and methods for large-scale SNP association studies
}
\details{
\tabular{ll}{
Package: \tab snpMatrix\cr
Version: \tab 1.5.5\cr
Date: \tab 2008-10-09\cr
Depends: \tab R(>= 2.3.0), survival, methods\cr
Suggests: \tab hexbin\cr
Enhances: \tab genetics\cr
License: \tab GNU General Public Licence (GPLv3)\cr
URL: \tab http://www-gene.cimr.cam.ac.uk/clayton/software/\cr
Collate: \tab ld.with.R ss.R contingency.table.R glm\_test.R ibs.stats.R imputation.R indata.R ld.snp.R ld.with.R.eml.R misc.R outdata.R plink.R qc.R qq\_chisq.R read.chiamo.R read.HapMap.R read.snps.pedfile.R single.R tdt\_single.R structure.R wtccc.sample.list.R wtccc.signals.R xstuff.R zzz.R\cr
LazyLoad: \tab yes\cr
biocViews: \tab Microarray, SNPsAndGeneticVariability\cr
Packaged: \tab Thu Oct  9 15:07:21 BST 2008\cr
Built: \tab R 2.7.1; x86\_64-redhat-linux-gnu; 2008-10-13 09:00:35; unix\cr
}

Index:
\preformatted{
X.snp-class             Class "X.snp"
X.snp.matrix-class      Class "X.snp.matrix"
chi.squared             Extract test statistics and p-values
epsout.ld.snp           Function to write an eps file directly to
                        visualize LD
families                Test data for family association tests
filter.rules            Filter a set of imputation rules
for.exercise            Data for exercise in use of the snpMatrix
                        package
genotype-class          snpMatrix-internal
glm.test.control        Set up control object for GLM tests
ibs.stats               function to calculate the identity-by-state
                        stats of a group of samples
ibsCount                Count alleles identical by state
ibsDist                 Distance matrix based on identity by state
                        (IBS)
imputation.maf          Extract MAF and R-squared values from
                        imputation rules
impute.snps             Impute snps
ld.snp                  Function to calculate pairwise D', r-squared
ld.with                 function to calculate the LD measures of
                        specific SNPs against other SNPs
misinherits             Find non-Mendelian inheritances in family data
pair.result.ld.snp      Function to calculate the pairwise D',
                        r-squared, LOD of a pair of specified SNPs
plot.snp.dprime         Function to draw the pairwise D' in a eps file
pool                    Pool test results from several studies or
                        sub-studies
pool2                   Pool results of tests from two independent
                        datasets
qq.chisq                Quantile-quantile plot for chi-squared tests
read.HapMap.data        function to import HapMap genotype data as
                        snp.matrix
read.pedfile.info       function to read the accompanying info file of
                        a LINKAGE ped file
read.pedfile.map        function to read the accompanying map file of a
                        LINKAGE ped file
read.plink              Read a PLINK binary data file as a snp.matrix
read.snps.chiamo        Read genotype data from the output of Chiamo
read.snps.long          Read SNP data in long format
read.snps.long.old      Read SNP input data in "long" format (old
                        version)
read.snps.pedfile       Read genotype data from a LINKAGE "pedfile"
read.wtccc.signals      read normalized signals in the WTCCC signal
                        file format
row.summary             Summarize rows of a snp matrix
single.snp.tests        1-df and 2-df tests for genetic associations
                        with SNPs (or imputed SNPs)
snp-class               Class "snp"
snp.cbind               Bind together two or more snp.matrix objects
snp.cor                 Correlations with columns of a snp.matrix
snp.dprime-class        Class "snp.dprime" for Results of LD
                        calculation
snp.imputation          Calculate regression-based imputation rules
snp.lhs.tests           Score tests with SNP genotypes as dependent
                        variable
snp.matrix-class        Class "snp.matrix"
snp.pre                 Pre- or post-multiply a snp.matrix object by a
                        general matrix
snp.reg.imputation-class
                        Class "snp.reg.imputation"
snp.rhs.tests           Score tests with SNP genotypes as independent
                        variable
snp.tests.glm-class     Class "snp.tests.glm"
snp.tests.single-class
                        Classes "snp.tests.single" and
                        "snp.tests.single.score"
snpMatrix-package       The snp.matrix and X.snp.matrix classes
switch.alleles          Switch alleles in columns of a snp.matrix or in
                        test results
tdt.snp                 1-df and 2-df tests for genetic associations
                        with SNPs (or imputed SNPs) in family data
test.allele.switch      Test for switch of alleles between two
                        collections
testdata                Test data for the snpMatrix package
write.snp.matrix        Write a snp.matrix object as a text file
wtccc.sample.list       read the sample list from the header of the
                        WTCCC signal file format
xxt                     X.X-transpose for a standardized snp.matrix
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{imputation-vignette} \tab imputation (source, pdf)\cr
\code{snpMatrix-vignette} \tab snpMatrix (source, pdf)\cr
}
}
\author{
David Clayton <david.clayton@cimr.cam.ac.uk> and
Hin-Tak Leung <htl10@users.sourceforge.net>

Maintainer: David Clayton <david.clayton@cimr.cam.ac.uk>
}
\keyword{ package }
