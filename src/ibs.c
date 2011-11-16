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

/* Identity-by-state; we only count calls and exclude no calls */

#include <R.h>
#include <Rinternals.h>

SEXP do_ibs(SEXP x)
{
  SEXP dims = R_NilValue;
  SEXP names = R_NilValue, rnames = R_NilValue;
  int rows = 0, cols =0;
  int i_first =0, i_second = 0;
  
  if(TYPEOF(x) != RAWSXP) {
    Rprintf(" input x wrong type\n");
    return R_NilValue;
  }
  
  PROTECT(dims = getAttrib(x, R_DimSymbol));
  if (length(dims) == 2)
    {
      rows = INTEGER(dims)[0];
      cols = INTEGER(dims)[1];
      Rprintf("Information: samples = %i, snps = %i\n", rows, cols);
    }
  else
    {
      Rprintf("wrong dims\n");
      UNPROTECT(1);
      return R_NilValue;
    }
  
  names = getAttrib(x, R_DimNamesSymbol);
  rnames = GetRowNames(names);

  int out_size = (rows * (rows - 1)) /2;
  SEXP out_count, out_fract, out_names;
  PROTECT(out_count   = allocVector(INTSXP, out_size));
  PROTECT(out_fract   = allocVector(REALSXP, out_size));
  PROTECT(out_names = allocVector(STRSXP, out_size));

  int i_element = 0;
  for (i_first = 0; i_first < rows - 1; i_first++) {
    for (i_second = i_first + 1; i_second < rows ; i_second++) {
      int j; 
      int count = 0, calls_valid = 0;
      char buffer[256];
      for(j = 0 ; j < cols ; j++) {
	int offset = j * rows;
	Rbyte char_first = RAW(x)[offset + i_first];
	Rbyte char_second = RAW(x)[offset + i_second];
	if ((char_first) && (char_second)) {
	  calls_valid++;
	  if (char_first  == char_second) {
	    count++;
	  }
	}
      }
      INTEGER(out_count)[i_element] = count;
      REAL(out_fract)[i_element] = ((double) count)/calls_valid;
      snprintf(buffer, 256, "%s,%s", CHAR(STRING_ELT(rnames, i_first)), CHAR(STRING_ELT(rnames, i_second)));
      SET_STRING_ELT(out_names, i_element, mkChar(buffer));
      i_element++;
    }
  }
  SEXP Result, Names;
  PROTECT(Result = allocVector(VECSXP, 2));
  PROTECT(Names = allocVector(STRSXP, 2));
  SET_VECTOR_ELT(Result, 0, out_count);
  SET_VECTOR_ELT(Result, 1, out_fract);
  SET_STRING_ELT(Names, 0, mkChar("Count"));
  SET_STRING_ELT(Names, 1, mkChar("Fraction"));
  setAttrib(Result, R_NamesSymbol, Names);
  setAttrib(Result, R_RowNamesSymbol, out_names);
  SEXP dfClass;
  PROTECT(dfClass = allocVector(STRSXP, 1));
  SET_STRING_ELT(dfClass, 0, mkChar("data.frame"));
  setAttrib(Result, R_ClassSymbol, dfClass);

  UNPROTECT(7);
  return Result;
}
