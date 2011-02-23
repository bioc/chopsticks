/*
 *  Copyright (C) 2007  Hin-Tak Leung
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available via WWW at
 *  http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
 *  writing to the Free Software Foundation, Inc., 51 Franklin Street
 *  Fifth Floor, Boston, MA 02110-1301  USA.
 */

#include <R.h>
#include <Rinternals.h>

#include "pairwise_linkage.h"

SEXP ld_with(SEXP x, SEXP snps, SEXP signed_r) {
  SEXP dims = R_NilValue;
  SEXP in_names = R_NilValue, cnames = R_NilValue; 
  SEXP dprime = R_NilValue, rmisc = R_NilValue, lod = R_NilValue, ans = R_NilValue;
  SEXP dimnames = R_NilValue, snpnames = R_NilValue;
  SEXP ans_name = R_NilValue;
  int rows = 0 , cols = 0;
  int need_signed_r = 0;
  int i=0,j=0;

  int n_select = LENGTH(snps); 

  if(TYPEOF(x) != RAWSXP)
    error(" input snp.data wrong type\n");
  if(TYPEOF(snps) != INTSXP)
    error(" input snps wrong type\n");
  if(TYPEOF(signed_r) != LGLSXP)
    error(" input signed_r wrong type\n");
  
  PROTECT(dims = getAttrib(x, R_DimSymbol));
  if (length(dims) == 2) {
    rows = INTEGER(dims)[0];
    cols = INTEGER(dims)[1];
    Rprintf("Information: The input contains %i samples with %i snps\n", rows, cols);
  } else {
    error("The input does not seem to have two dimensions\n");
  }
  
  if(LOGICAL(signed_r)[0]) {
    need_signed_r = 1;
  }
  /* we'll need the column names later */
  in_names = getAttrib(x, R_DimNamesSymbol);
  cnames = GetColNames(in_names);

  /* finished playing with the input, now do some real work */
  
  PROTECT(dprime = allocMatrix(REALSXP, cols, n_select));
  PROTECT(rmisc  = allocMatrix(REALSXP, cols, n_select));
  PROTECT(lod    = allocMatrix(REALSXP, cols, n_select));

  PROTECT(dimnames = allocVector(VECSXP, 2));
  PROTECT(snpnames = allocVector(STRSXP, n_select));
  for(i = 0; i < n_select; i++)
    SET_STRING_ELT(snpnames, i, mkChar(CHAR(STRING_ELT(cnames, INTEGER(snps)[i]))));
  SET_VECTOR_ELT(dimnames, 0, duplicate(cnames));
  SET_VECTOR_ELT(dimnames, 1, duplicate(snpnames));
  setAttrib(dprime, R_DimNamesSymbol, dimnames);
  setAttrib(rmisc,  R_DimNamesSymbol, duplicate(dimnames));
  setAttrib(lod,    R_DimNamesSymbol, duplicate(dimnames));

  memset(REAL(dprime), 0, cols * n_select * sizeof(double));
  memset(REAL(rmisc) , 0, cols * n_select * sizeof(double));
  memset(REAL(lod)   , 0, cols * n_select * sizeof(double));

  for (i = 0 ; i < cols ; i++) {
    for(j = 0 ; j < n_select ; j++) {
      int offset = i + j * cols;
      geno_cptr res = get_geno_count(RAW(x) + i * rows,
				     RAW(x) + INTEGER(snps)[j] * rows ,
				     rows);
      REAL(dprime)[offset] = res->dprime;
      if (need_signed_r) {
	if (!ISNA(res->rsq2)) {
	  REAL(rmisc)[offset] = res->sign_of_r * sqrt(res->rsq2);
	} else {
	  REAL(rmisc)[offset] = NA_REAL;
	}
      } else {
	REAL(rmisc)[offset] = res->rsq2;
      }
      REAL(lod)[offset] = res->lod;
      free(res->expt);
      free(res);
    }
  }

  /* constructing the rest of the result */
  PROTECT(ans    = allocVector(VECSXP, 3));
  PROTECT(ans_name = allocVector(STRSXP, 3));
  SET_STRING_ELT(ans_name, 0, mkChar("dprime"));
  if (need_signed_r)
    {
      SET_STRING_ELT(ans_name, 1, mkChar("r"));
    }
  else
    {
      SET_STRING_ELT(ans_name, 1, mkChar("rsq2"));
    }
  SET_STRING_ELT(ans_name, 2, mkChar("lod"));
  setAttrib(ans, R_NamesSymbol, ans_name);

  SET_VECTOR_ELT(ans, 0, dprime);
  SET_VECTOR_ELT(ans, 1, rmisc);
  SET_VECTOR_ELT(ans, 2, lod);

  UNPROTECT(8);
  return ans;
}
