/* 
  If A is a snp.matrix these routines calculate the matrices B.M or M.B, 
  where M is a general matrix and B is derived from A by normalising 
  columns to have zero mean and unit standard deviation under HWE. 
  That is, if p is the allele frequency for one column of A, the 
  corresponding column of B (if elements are coded 0, 1, or 2) is 
  obtained by subtracting the mean, 2*p and dividing by the SD, 
  sqrt(2*p*(1-p)). For male samples and the X chromosome, codes are 
  0 and 2 so that the mean is again 2*p, but the SD is now 2*sqrt(p*(1-p)).
  Missing genotypes score zero (equivalent to replacing missing values by 
  the mean in the original matrix. Missing genotypes are replaced by 
  their (marginal) expectations - i.e. twice the allele frequency

  Allele frequencies may be taken from the data or supplied as an argument
*/

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Rmissing.h"

SEXP snp_pre(const SEXP Snps, const SEXP Mat, const SEXP Frequency) {
  
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
  SEXP Snpnames = VECTOR_ELT(getAttrib(Snps, R_DimNamesSymbol), 1);

  cl = GET_CLASS(Mat);
  if (TYPEOF(cl) != STRSXP) {
    cl = R_data_class(Mat, FALSE); /* S4 way of getting class attribute */
  }
  if (strcmp(CHAR(STRING_ELT(cl, 0)), "matrix"))
    error("Argument error - Mat wrong type");
  dim = INTEGER(getAttrib(Mat, R_DimSymbol));
  if (dim[1]!=N)
    error("non-conformable arguments");
  int P = dim[0];
  double *mat = REAL(Mat);
  SEXP Rownames = GetRowNames(Mat);

  /* Allele frequencies */

  double *frequency = NULL;
  if (TYPEOF(Frequency) == REALSXP) {
    if (LENGTH(Frequency)!=M)
      error("incorrect length for allele frequency vector");
    frequency = REAL(Frequency);
  }
  else if (TYPEOF(Frequency) != NILSXP)
    error("Argument error: Frequency is wrong type");
     
  /* Result matrix */

  SEXP Result, Dimnames;
  PROTECT(Result = allocMatrix(REALSXP, P, M));
  PROTECT(Dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(Dimnames, 0, duplicate(Rownames));
  SET_VECTOR_ELT(Dimnames, 1, duplicate(Snpnames));
    
  double *result = REAL(Result);
  memset(result, 0x00, P*M*sizeof(double));

  /* Update result matrix for each locus in turn */

  for (int  k=0, ik=0, jks=0; k<M; k++, jks+=P) {

    /* Get allele frequency */

    double afk = NA_REAL;
    if (frequency) 
      afk = frequency[k];
    else {     
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
      if (s1)
	afk = ((double) s2) / ((double) s1);
    }
    
    /* If polymorphic, add contribution */

    if (afk!=NA_REAL && afk>0.0 && afk<1.0) {
      double mean = 2.0*afk + 1.0;
      double sd2 = sqrt(2.0*afk*(1.0-afk));
      double sd1 = 2.0*sqrt(afk*(1.0-afk));

      for (int i=0, jis=0; i<N; i++, ik++, jis+=P) {
	int sik = (int) snps[ik];
	if (sik) {
	  double xik;
	  if (ifFemale && !ifFemale[i])
	    xik = ((double) sik - mean)/sd1;
	  else
	    xik = ((double) sik - mean)/sd2;
	  for (int j=0, jk=jks, ji=jis; j<P; j++, jk++, ji++) {
	    result[jk] += xik*mat[ji];
	  }
	}
      }
    }
    else {
      ik += N;
    }
  }
      
  UNPROTECT(2);
  return(Result);
}

SEXP snp_post(const SEXP Snps, const SEXP Mat, const SEXP Frequency) {
  
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
  SEXP Sampnames = VECTOR_ELT(getAttrib(Snps, R_DimNamesSymbol), 0);

  cl = GET_CLASS(Mat);
  if (TYPEOF(cl) != STRSXP) {
    cl = R_data_class(Mat, FALSE); /* S4 way of getting class attribute */
  }
  if (strcmp(CHAR(STRING_ELT(cl, 0)), "matrix"))
    error("Argument error - Mat wrong type");
  dim = INTEGER(getAttrib(Mat, R_DimSymbol));
  if (dim[0]!=M)
    error("non-conformable arguments");
  int P = dim[1];
  double *mat = REAL(Mat);
  SEXP Colnames = GetColNames(Mat);
     
  /* Allele frequencies */

  double *frequency = NULL;
  if (TYPEOF(Frequency) == REALSXP) {
    if (LENGTH(Frequency)!=M)
      error("incorrect length for allele frequency vector");
    frequency = REAL(Frequency);
  }
  else if (TYPEOF(Frequency) != NILSXP)
    error("Argument error: Frequency is wrong type");
     
  /* Result matrix */

  SEXP Result, Dimnames;
  PROTECT(Result = allocMatrix(REALSXP, N, P));
  PROTECT(Dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(Dimnames, 0, duplicate(Sampnames));
  SET_VECTOR_ELT(Dimnames, 1, duplicate(Colnames));
    
  double *result = REAL(Result);
  memset(result, 0x00, N*P*sizeof(double));

  /* Update result matrix for each locus  in turn */

  for (int  k=0, ik=0, jks=0; k<M; k++, jks+=P) {

      /* Get allele frequency */

    double afk = NA_REAL;
    if (frequency) 
      afk = frequency[k];
    else {     
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
      if (s1)
	afk = ((double) s2) / ((double) s1);
    }
    
    /* If polymorphic, add contribution */

    if (afk!=NA_REAL && afk>0.0 && afk<1.0) {
      double mean = 2.0*afk + 1.0;
      double sd2 = sqrt(2.0*afk*(1.0-afk));
      double sd1 = 2.0*sqrt(afk*(1.0-afk));

      for (int i=0; i<N; i++, ik++) {
	int sik = (int) snps[ik];
	if (sik) {
	  double xik;
	  if (ifFemale && !ifFemale[i])
	    xik = ((double) sik - mean)/sd1;
	  else
	    xik = ((double) sik - mean)/sd2;
	  for (int j=0, ij=i, kj=k; j<P; j++, ij+=N, kj+=M) {
	    result[ij] += xik*mat[kj];
	  }
	}
      }
    }
    else {
      ik += N;
    }
  }
      
  UNPROTECT(2);
  return(Result);
}
  

