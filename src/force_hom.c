#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "Rmissing.h"

SEXP force_hom(const SEXP Xsnps, const SEXP Female) {

  int *female = LOGICAL(Female);
  int N, M;
  int *dim = INTEGER(getAttrib(Xsnps, R_DimSymbol));
  N = dim[0];
  M = dim[1];


  SEXP Result;
  PROTECT(Result = duplicate(Xsnps));
  unsigned char *result = RAW(Result);
  for (int i=0; i<N; i++) {
    int male = !female[i];
    for (int j=0, ij=i; j<M; j++, ij+=N) {
      if (male && result[ij]==2)
	result[ij] = 0;
    }
  }
  UNPROTECT(1);
  return(Result);
}
  
