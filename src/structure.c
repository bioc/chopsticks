/* 
  If A is a snp.matrix this routine calculates the matrix B.B-transpose, 
  where B is derived from A by normalising columns to have zero mean and 
  unit standard deviation under HWE. That is, if p is the allele frequency 
  for one column of A, the corresponding column of B (if elements are
  coded 0, 1, or 2) is obtained by subtracting the mean, 2*p and dividing by 
  the SD, sqrt(2*p*(1-p)). For male samples and the X chromosome, codes are 
  0 and 2 so that the mean is again 2*p, but the SD is now 2*sqrt(p*(1-p)).
  Missing genotypes score zero (equivalent to replacing missing values by 
  the mean in the original matrix.

  Missing data treatment:

  If Correct_for_missing is FALSE, missing genotypes are replaced by 
  their (marginal) expectations - i.e. twice the allele frequency

  If TRUE, contributions to the output matrix are weighted using inverse
  probability weights. The (small) probability that locus k is missing in 
  subject i is assumed to be mu*alpha_i*beta_k where alpha and beta are 
  are vectors with mean 1. Then the probability that subject i is observed 
  at locus k is 1-mu*alpha_i*beta_k, and the probability that locus k 
  is observed in both subject i and subject j is 
  (1-mu*alpha_i*beta_k)*(1-mu*alpha_j*beta_k) 

  alpha_i and beta_k are estimated from the observed numbers of missing calls 
  as follows:

  T_ik = 1 if call (i,k) is missing, 0 otherwise

  alpha_i = N*T_i./T_..
  beta_k = M*T_.k/T..
  mu = T_../(N*M)

  so that

  mu*alpha_i*beta_k = T_i.*T_.k/T_..

  

*/

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Rmissing.h"

SEXP xxt(const SEXP Snps, const SEXP Correct_for_missing, 
	 const SEXP Lower_only) {
  
  if (TYPEOF(Correct_for_missing)!=LGLSXP)
    error("Argument error - Correct_for_missing wrong type");
  const int correct = *LOGICAL(Correct_for_missing);

  if (TYPEOF(Lower_only)!=LGLSXP)
    error("Argument error - Lower_only wrong type");
  const int lower = *LOGICAL(Lower_only);

  int *ifFemale = NULL;
  SEXP cl = GET_CLASS(Snps);
  if (TYPEOF(cl) != STRSXP) {
    cl = R_data_class(Snps, FALSE); /* S4 way of getting class attribute */
  }
  if (!strcmp(CHAR(STRING_ELT(cl, 0)), "X.snp.matrix")) {
    SEXP Female = R_do_slot(Snps, mkString("Female"));
    if (TYPEOF(Female)!=LGLSXP)
      error("Argument error -  Female slot wrong type");
    ifFemale = LOGICAL(Female);
  }
  else if (strcmp(CHAR(STRING_ELT(cl, 0)), "snp.matrix")) {
    error("Argument error - Snps wrong type");
  }    

  const unsigned char *snps = RAW(Snps);
  int N, M;
  int *dim = INTEGER(getAttrib(Snps, R_DimSymbol));
  N = dim[0];
  M = dim[1];

  /* Weights for missing data correction */

  int *Ti=NULL, *Tk=NULL, T=0;
  if (correct) {
    warning("With correct.for.missing option set, result may not be a positive semi-definite matrix");
    T = 0;
    Ti = (int *) calloc(N, sizeof(int));
    memset(Ti, 0x00, N*sizeof(int));
    Tk = (int *) calloc(M, sizeof(int));
    memset(Tk, 0x00, M*sizeof(int));
    for (int k=0, ik=0; k<M; k++) {
      for (int i=0; i<N; i++) {
	int sik = (int) snps[ik++];
	if (!sik){
	  T++;
	  Ti[i]++;
	  Tk[k]++;
	}
      }
    } 
  }

  /* Result matrix */

  SEXP Result;
  PROTECT(Result = allocMatrix(REALSXP, N, N));
  double *result = REAL(Result);
  memset(result, 0x00, N*N*sizeof(double));

  /* Update result matrix for each locus in turn */

  for (int ik=0, k=0; k<M; k++) {

    /* Calculate allele frequency */

    int s1=0, s2=0;
    for (int ki=ik, i=0; i<N; i++) {
      int w = (int) snps[ki++];
      if (w) {
	w--;
	if (ifFemale && !ifFemale[i]) {
	  s1++;
	  s2 += w/2;
	}
	else {
	  s1 += 2;
	  s2 += w;
	}
      }
    }
    
    /* If polymorphic, add contribution */

    if ((s1>0) && (s2>0) && (s2<s1)) {

      double afk = ((double) s2) / ((double) s1);
      double mean = 2.0*afk + 1.0;
      double sd2 = sqrt(2.0*afk*(1.0-afk));
      double sd1 = 2.0*sqrt(afk*(1.0-afk));
      double tk=0.0, ipw=0.0;
      if (correct)
	tk = T? (double) Tk[k] / (double) T: 0.0;

      /* Update X.X-transpose matrix */

      for (int i=0, ij=0; i<N; i++, ik++) {
	int sik = (int) snps[ik];
	if (sik) {
	  double xik;
	  if (ifFemale && !ifFemale[i])
	    xik = ((double) sik - mean)/sd1;
	  else
	    xik = ((double) sik - mean)/sd2;
	  if (correct) {
	    ipw = 1.0/(1.0 - tk*(double)Ti[i]); /* IPW for diagonal */
	  }
	  ij += i;
	  for (int jk=ik, j=i; j<N; j++, ij++) {
	    int sjk = (int) snps[jk++]; 
	    if (sjk) {
	      double xjk;
	      if (ifFemale && !ifFemale[j])
		xjk =  ((double) sjk - mean)/sd1;
	      else
		xjk =  ((double) sjk - mean)/sd2;
	      if (correct) 
		result[ij] += xik*xjk*((i==j)?ipw: ipw/(1.0-tk*(double)Ti[j]));
	      else 
		result[ij] += xik*xjk;
	    }
	  }
	}
	else {
	  ij += N;
	}
      }
    }
    else {
      ik += N;
    }
  }

  /* Copy lower triangle to upper triangle */

  if (!lower) {
    for (int i=0, ij=0; i<N; i++) {
      ij += (i+1);
      for (int j=i+1, ji=ij-1+N; j<N; j++, ij++, ji+=N) {
	result[ji] = result[ij];
      }
    }
  }

  /* Return work space */

  if (correct) {
    free(Tk);
    free(Ti);
  }

      
  UNPROTECT(1);
  return(Result);
}
  

/* 
   Correlations between columns of a snpmatrix and columns of a normal matrix
*/

 
SEXP corsm(const SEXP Snps, const SEXP X) { 

  if (!inherits(Snps, "snp.matrix"))
    error("Argument error - Snps wrong type");
  const unsigned char *snps = RAW(Snps);
  int *dim = INTEGER(getAttrib(Snps, R_DimSymbol));
  int N = dim[0];
  int M = dim[1];

  if (TYPEOF(X)!=REALSXP)
    error("Argument error - X wrong type");
  if (X == R_NilValue) {
    error("Argument error - X = NULL");
  }
  const double *x = REAL(X);
  dim = INTEGER(getAttrib(X, R_DimSymbol));
  if (dim[0] != N) {
    error("Unequal numbers of rows");
  }
  int P = dim[1];


  SEXP Result;
  PROTECT(Result = allocMatrix(REALSXP, M, P));
  double *result = REAL(Result);

  for (int j=0, ij=0, jks=0; j<P; j++, jks+=N) {
    for (int i=0, ik=0; i<M; i++, ij++) {
      double sg=0.0, sgg=0, sx=0.0, sxx=0.0, sgx=0.0;
      int s=0;
      for (int k=0, jk=jks; k<N; k++) {
	int g = (int) snps[ik++];
	double xk = x[jk++];
	if (g && !ISNA(xk)) {
	  s++;
	  sg += g;
	  sgg += g*g;
	  sx += xk;
	  sxx += xk*xk;
	  sgx += g*xk;
	}
      }
      if (s) {
	sgg -= sg*sg/(double)s;
	sxx -= sx*sx/(double)s;
	sgx -= sg*sx/(double)s;
	if (sgg>0.0 && sxx>0.0)
	  result[ij] = sgx/sqrt(sgg*sxx);
	else
	  result[ij] = NA_REAL;
      }
      else {
	result[ij] = NA_REAL;
      }
    }
  }
  
  UNPROTECT(1);
  return(Result);
}

/*
  IBS matrix.

  Calculates NxN matrix of type integer, whose upper triangle counts number 
  non-missing pairs of chromosomes and whose lower triangle counts number 
  IBS. The diagonal counts non-missing calls for each subject

  For autosomes, each pair of non-missing SNPs add 4 above the diagonal, and 
  0, 2, or 4 below the diagonal, according to IBS state:
                      AA AB BB
		   AA  4  2  0
                   AB  2  2  2
                   BB  0  2  4
  For an X SNP, pairs of females count similarly but pairs of males count 1 
  above the diagonal and a male/femal pair counts 2. Below the diagonal we 
  have:
       AY BY               AA AB BB
    AY  1  0      or    AY  2  1  0
    BY  0  1            BY  0  1  2

  However this version does not currently does not deal with X.snp.matrices

*/

SEXP ibs_count(const SEXP Snps) { 

  int *ifFemale = NULL;
  SEXP cl = GET_CLASS(Snps);
  if (TYPEOF(cl) != STRSXP) {
    cl = R_data_class(Snps, FALSE); /* S4 way of getting class attribute */
  }
  if (!strcmp(CHAR(STRING_ELT(cl, 0)), "X.snp.matrix")) {
    SEXP Female = R_do_slot(Snps, mkString("Female"));
    if (TYPEOF(Female)!=LGLSXP)
      error("Argument error -  Female slot wrong type");
    ifFemale = LOGICAL(Female);
  }
  else if (strcmp(CHAR(STRING_ELT(cl, 0)), "snp.matrix")) {
    error("Argument error - Snps wrong type");
  }    
  SEXP names = getAttrib(Snps, R_DimNamesSymbol);
  if (names == R_NilValue) {
    error("Argument error - Snps object has no names");
  }
  SEXP rowNames = VECTOR_ELT(names, 0);
  if (rowNames == R_NilValue) {
    error("Argument error - Snps object has no row names");
  }

  const unsigned char *snps = RAW(Snps);
  int N, M;
  int *dim = INTEGER(getAttrib(Snps, R_DimSymbol));
  N = dim[0];
  M = dim[1];


  /* Result matrix */

  SEXP Result, dimnames;
  PROTECT(Result = allocMatrix(INTSXP, N, N));
  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 0, duplicate(rowNames));    
  SET_VECTOR_ELT(dimnames, 1, duplicate(rowNames));    
  setAttrib(Result, R_DimNamesSymbol, dimnames);    
  int *result = INTEGER(Result);
  memset(result, 0x00, N*N*sizeof(int));

  /* Update result matrix for each locus in turn */

  for (int ik=0, k=0; k<M; k++) {

    /* Update IBS matrix */

    int N1 = N+1;
    for (int i=0, ii=0; i<N; i++, ii+=N1) {
      int base_div;
      if (ifFemale && !ifFemale[i])
	base_div = 2;
      else
	base_div = 1;
      int sik = (int) snps[ik++];
      if (sik) {
	result[ii]++;
	for (int j=i+1, jk=ik, ji=ii+1, ij=ii+N; 
	     j<N; j++, ji++, ij+=N) {
	  int div = base_div ;
	  if (ifFemale && !ifFemale[j])
	    div *= 2;
	  int sjk = (int) snps[jk++]; 
	  if (sjk) {
   	    int add = sjk==sik?
	      (sik==2 ? 2: 4) : (4-2*abs(sik-sjk));
	    result[ij] += add/div;
	    result[ji] += 4/div;
	  }
	}
      }
    }
  }
  
  UNPROTECT(2);
  return(Result);
}

/*
  Distance matrix based on IBS counts
*/

SEXP ibs_dist(const SEXP Ibsc) {
  
  if (!IS_INTEGER(Ibsc))
    error("Input object is not an integer array");

  const int *ibsc = INTEGER(Ibsc);
  int N, M;
  int *dim = INTEGER(getAttrib(Ibsc, R_DimSymbol));
  N = dim[0];
  M = dim[1];
  if (!N || N!=M) 
    error("Input object is not a square matrix");
  SEXP names = getAttrib(Ibsc, R_DimNamesSymbol);
  if (names == R_NilValue) {
    error("Argument error - no names");
  }
  SEXP rowNames = VECTOR_ELT(names, 0);
  if (rowNames == R_NilValue) {
    error("Argument error - no sample identifiers");
  }

  /* Result matrix */

  R_len_t Nout = (N*(N-1))/2;
  SEXP Result, Size, Class;
  PROTECT(Result = allocVector(REALSXP, Nout));
  PROTECT(Size = allocVector(INTSXP, 1));
  INTEGER(Size)[0] = N;
  PROTECT(Class = allocVector(STRSXP, 1));
  SET_STRING_ELT(Class, 0, mkChar("dist"));
  setAttrib(Result, install("Size"), Size);
  setAttrib(Result, install("Labels"), duplicate(rowNames));
  classgets(Result, Class);
  double *result = REAL(Result);
  memset(result, 0x00, Nout*sizeof(double));
  for (int i=1, ii=0, k=0; i<N; i++, ii+=(N+1)) {
    for (int j=i, ji=ii+1, ij=ii+N; j<N; j++, ji++, ij+=N){
      result[k++] = (double) (ibsc[ji]-ibsc[ij])/(double) ibsc[ji];
    }
  }

  UNPROTECT(3);
  return(Result);
}


