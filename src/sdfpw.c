/*
 *  Copyright (C) 2006  Hin-Tak Leung
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

/* see sdfpw.h for description of these routines */

#include <R.h>
#include <Rinternals.h>

#include <stdio.h>

#include "pairwise_linkage.h"
#include "ld_graphic_eps.h"
#include "hash_index.h"

#include "sdfpw.h"

/* "getListElement" was taken from R extension manual */
/* static so as to hide from the rest of the library */

/* get the list element named str, or return NULL */
static SEXP getListElement(SEXP list, char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

SEXP snp_pairwise(SEXP v, SEXP i, SEXP j)
{
  SEXP dims = R_NilValue;
  int rows = 0, cols =0;
  int ii, jj;
  geno_cptr res = NULL;
  

  if(TYPEOF(v) != RAWSXP)
    Rprintf(" input v wrong type\n");

  /* literals come in as numeric/double */
  PROTECT(i = coerceVector(i, INTSXP));
  PROTECT(j = coerceVector(j, INTSXP));
  if(TYPEOF(i) != INTSXP)
    Rprintf(" input i wrong type\n");
  if(TYPEOF(j) != INTSXP)
    Rprintf(" input j wrong type\n");

  PROTECT(v);

  PROTECT(dims = getAttrib(v, R_DimSymbol));
  if (length(dims) == 2)
    {
      rows = INTEGER(dims)[0];
      cols = INTEGER(dims)[1];
      Rprintf("Information: samples = %i, snps = %i\n", rows, cols);
    }
  else
    {
      Rprintf("wrong size\n");
    }

  ii = INTEGER(i)[0];
  jj = INTEGER(j)[0];

  Rprintf(" ii = %i, jj = %i\n", ii, jj);

  res = get_geno_count(RAW(v) + (ii - 1) * rows, RAW(v) + (jj - 1) * rows, rows);
  
  Rprintf(" %i\t%i\t%i\t%f\t%f\t%f\n", res->count[0], res->count[1], res->count[2], res->expt[0], res->expt[1], res->expt[2]);
  Rprintf(" %i\t%i\t%i\t%f\t%f\t%f\n", res->count[3], res->count[4], res->count[5], res->expt[3], res->expt[4], res->expt[5]);
  Rprintf(" %i\t%i\t%i\t%f\t%f\t%f\n", res->count[6], res->count[7], res->count[8], res->expt[6], res->expt[7], res->expt[8]);

  Rprintf("d\' = %f , r^2 = %f, lod= %f\n", res->dprime, res->rsq2, res->lod);
  free(res->expt);
  free(res);
  UNPROTECT(4);
  return R_NilValue;
}

SEXP snp_pair_graphics(SEXP v, SEXP fileoutput, SEXP i, SEXP j, SEXP depth, SEXP do_notes)
{
  const char *filename = NULL;
  SEXP dims = R_NilValue;
  SEXP names = R_NilValue, cnames = R_NilValue;
  int rows = 0 , cols = 0;
  int i_notes = 0;
  int ii, jj;
  int i_depth;
  int idx_i, idx_j;
  ld_graphics_ptr graphic_handle = NULL;

  if(TYPEOF(v) != RAWSXP)
    Rprintf(" input v wrong type\n");

  /* literals come in as numeric/double */
  PROTECT(i = coerceVector(i, INTSXP));
  PROTECT(j = coerceVector(j, INTSXP));
  PROTECT(fileoutput = coerceVector(fileoutput, STRSXP));
  PROTECT(depth = coerceVector(depth, INTSXP));
  PROTECT(do_notes = coerceVector(do_notes, INTSXP));
  if(TYPEOF(i) != INTSXP)
    Rprintf(" input i wrong type\n");
  if(TYPEOF(j) != INTSXP)
    Rprintf(" input j wrong type\n");
  if(TYPEOF(fileoutput) != STRSXP)
    Rprintf(" input filename wrong type\n");
  if(TYPEOF(depth) != INTSXP)
    Rprintf(" input depth wrong type\n");

  PROTECT(v);

  PROTECT(dims = getAttrib(v, R_DimSymbol));
  if (length(dims) == 2)
    {
      rows = INTEGER(dims)[0];
      cols = INTEGER(dims)[1];
      Rprintf("Information: samples = %i, snps = %i\n", rows, cols);
    }
  else
    {
      Rprintf("wrong size\n");
    }

  ii = INTEGER(i)[0];
  jj = INTEGER(j)[0];
  if (INTEGER(do_notes)[0] != 0)
    {
      i_notes = 1;
    }
  if (jj > cols)
    {
      jj = cols;
    }
  i_depth = INTEGER(depth)[0];


  names = getAttrib(v, R_DimNamesSymbol);
  cnames = GetColNames(names);

  filename = CHAR(STRING_ELT(fileoutput, 0));
  Rprintf("Writing to %s ...", filename);
  graphic_handle = graphic_init(filename, ii, jj, i_depth, 0, 0);
  if (!graphic_handle) {
    Rprintf("Cannot open %s for writing\n", filename);
    return R_NilValue;
  }

  for (idx_i = ii -1 ; idx_i <= jj -1 ; idx_i ++)
    {
      graphic_do_name(graphic_handle, idx_i, CHAR(STRING_ELT(cnames, idx_i)));
    }

  for(idx_j = 0; idx_j < i_depth; idx_j++)
    {
      graphic_scan_line_begin(graphic_handle, idx_j);
      for (idx_i = ii - 1 ; idx_i <= jj - 2 - idx_j; idx_i++)
	{	  
	  if (0) Rprintf(" start = %i, end = %i\n", idx_i, idx_i + idx_j + 1);
	  graphic_do_pair(graphic_handle, RAW(v) + idx_i * rows, RAW(v) + (idx_i + idx_j + 1) * rows , idx_i, idx_j, rows, i_notes);
	}
      graphic_scan_line_end(graphic_handle, idx_j);
    }
  UNPROTECT(7);
  graphic_close(graphic_handle);
  Rprintf("... Done\n");  
  return R_NilValue;
}

SEXP snp_dprime_draw(SEXP list_in, SEXP fileoutput, SEXP color_scheme, SEXP do_notes, SEXP metric)
{
  SEXP dprime = R_NilValue;
  SEXP rmisc = R_NilValue, lod = R_NilValue; /* we are storing r or rsq2 in rmisc */
  SEXP snp_names = R_NilValue;
  SEXP dims = R_NilValue;
  int rows = 0 , cols = 0;
  int i_notes = 0;
  int idx_i, idx_j;
  const char *filename = NULL;
  ld_graphics_ptr graphic_handle = NULL;
  index_db dist_db = NULL;
  int is_r = 0;
  int do_metric = 0;

  PROTECT(list_in);
  PROTECT(do_notes = coerceVector(do_notes, INTSXP));
  PROTECT(color_scheme = coerceVector(color_scheme, INTSXP));
  if(TYPEOF(list_in) != VECSXP)
    Rprintf("list in wrong type\n");
  if(TYPEOF(fileoutput) != STRSXP)
    Rprintf("filename in wrong type\n");
  if ((metric != R_NilValue) && (TYPEOF(metric) != INTSXP)) {
    Rprintf("metric in wrong type\n");
  }

  dprime = getListElement(list_in, "dprime");
  rmisc   = getListElement(list_in, "rsq2");
  lod    = getListElement(list_in, "lod");

  if (rmisc == R_NilValue)
    {
      rmisc = getListElement(list_in, "r");
      is_r = 1;
    }

  if ((TYPEOF(dprime) != REALSXP) || (TYPEOF(rmisc) != REALSXP) || (TYPEOF(lod) != REALSXP))
    {
      Rprintf("filename in wrong type\n");
      return R_NilValue;
    }
  
  PROTECT(dims = getAttrib(dprime, R_DimSymbol));
  if (length(dims) == 2)
    {
      rows = INTEGER(dims)[0];
      cols = INTEGER(dims)[1];
      Rprintf("Information: range = %i, depth = %i\n", rows, cols);
    }
  else
    {
      Rprintf("wrong size\n");
    }

  if (INTEGER(do_notes)[0] != 0)
    {
      i_notes = 1;
    }

#define MAX(x,y) ((x) > (y) ? (x):(y))
#define MIN(x,y) ((x) < (y) ? (x):(y))

  filename = CHAR(STRING_ELT(fileoutput, 0));
  Rprintf("Writing to %s ...", filename);

  if (metric != R_NilValue) {
    do_metric = 1;
  }
  graphic_handle = graphic_init(filename, 1, rows + 1, cols, INTEGER(color_scheme)[0] , do_metric);

  if (do_metric) {
    int max_metric = INT_MIN;
    int min_metric = INT_MAX;
    int i_metric;
    SEXP metric_names = getAttrib(metric, R_NamesSymbol);
    dist_db = index_create(LENGTH(metric));
    do_metric = 1;
    for(i_metric = 0; i_metric < LENGTH(metric) ; i_metric++) {
      if (INTEGER(metric)[i_metric] != NA_INTEGER) {
	min_metric = MIN(INTEGER(metric)[i_metric], min_metric);
	max_metric = MAX(INTEGER(metric)[i_metric], max_metric);
	index_insert(dist_db, CHAR(STRING_ELT(metric_names, i_metric)), INTEGER(metric)[i_metric]);
      }
    }
    if (max_metric> min_metric) {
      graphic_add_metric(graphic_handle, min_metric, max_metric - min_metric);
    }
  }
  
  snp_names = getAttrib(list_in, install("snp.names"));
  
  if (snp_names != R_NilValue)
    {
      if (length(snp_names) == rows + 1)
	{
	  
	  for (idx_i = 0 ; idx_i < rows + 1 ; idx_i ++)
	    {
	      int dist_x = -1;
	      graphic_do_name(graphic_handle, idx_i, CHAR(STRING_ELT(snp_names, idx_i)));
	      if ((do_metric) && ((dist_x = index_lookup(dist_db, CHAR(STRING_ELT(snp_names, idx_i)))) >= 0)) {	
		graphic_do_metric(graphic_handle, idx_i, dist_x);
	      }
	    }
	}
      else
	{
	  Rprintf("size of snp.names doesn't agree with size of dprime data, not doing names");
	}

    }
  for (idx_j = 0 ; idx_j < cols; idx_j++)
    {	  
      graphic_scan_line_begin(graphic_handle, idx_j);
      for(idx_i = 0; idx_i < rows - idx_j; idx_i++)
	{
	  geno_cptr res = (geno_cptr)calloc(1, sizeof(geno_count));
	  double rmaybe    = REAL(rmisc)[idx_i + rows * idx_j];
	  res->dprime = REAL(dprime)[idx_i + rows * idx_j];
	  if (is_r) 
	    {
	      if (rmaybe < -1.1) /* invalid value of r is -2, but is less than -1 anyway */
		{
		  res->rsq2 = -1; /* invalid value */
		}
	      else
		{
		  res->rsq2 = rmaybe * rmaybe;
		}
	    }
	  else
	    {
	      res->rsq2 = rmaybe;
	    }
	  res->lod     = REAL(lod)[idx_i + rows * idx_j];
	  graphic_draw_pair(graphic_handle, res, idx_i, idx_j, i_notes);
	  free(res);
	}
      graphic_scan_line_end(graphic_handle, idx_j);
    }
  UNPROTECT(4);
  graphic_close(graphic_handle);
  Rprintf("... Done\n");  
  return R_NilValue;
}

SEXP snp_pair_range(SEXP v, SEXP i, SEXP j, SEXP depth, SEXP signed_r)
{
  SEXP dims = R_NilValue;
  SEXP in_names = R_NilValue, cnames = R_NilValue;
  /* rmisc is either r or rsq2 depending on context */
  SEXP dprime = R_NilValue, rmisc = R_NilValue, lod = R_NilValue, ans = R_NilValue;
  SEXP ans_name = R_NilValue;
  SEXP snp_names = R_NilValue;
  SEXP class_name = R_NilValue;
  int rows = 0 , cols = 0;
  int ii, jj;
  int i_depth;
  int idx_i, idx_j;
  int width = 0; /* width = range -1 */
  int need_signed_r = 0;

  if(TYPEOF(v) != RAWSXP)
    Rprintf(" input v wrong type\n");

  /* literals come in as numeric/double */
  PROTECT(i = coerceVector(i, INTSXP));
  PROTECT(j = coerceVector(j, INTSXP));
  PROTECT(depth = coerceVector(depth, INTSXP));
  PROTECT(signed_r = coerceVector(signed_r, LGLSXP));
  if(TYPEOF(i) != INTSXP)
    Rprintf(" input i wrong type\n");
  if(TYPEOF(j) != INTSXP)
    Rprintf(" input j wrong type\n");
  if(TYPEOF(depth) != INTSXP)
    Rprintf(" input depth wrong type\n");
  if(TYPEOF(signed_r) != LGLSXP)
    Rprintf(" input signed_r wrong type\n");

  PROTECT(v);

  PROTECT(dims = getAttrib(v, R_DimSymbol));
  if (length(dims) == 2)
    {
      rows = INTEGER(dims)[0];
      cols = INTEGER(dims)[1];
      Rprintf("Information: The input contains %i samples with %i snps\n", rows, cols);
    }
  else
    {
      Rprintf("wrong size\n");
    }

  ii = INTEGER(i)[0];
  jj = INTEGER(j)[0];
  if (jj > cols)
    {
      jj = cols;
    }
  i_depth = INTEGER(depth)[0];
  width = jj - ii;

  if(LOGICAL(signed_r)[0])
    {
      need_signed_r = 1;
    }

  PROTECT(dprime = allocMatrix(REALSXP, width, i_depth));
  PROTECT(rmisc  = allocMatrix(REALSXP, width, i_depth));
  PROTECT(lod    = allocMatrix(REALSXP, width, i_depth));

  memset(REAL(dprime), 0, width * i_depth * sizeof(double));
  memset(REAL(rmisc) , 0, width * i_depth * sizeof(double));
  memset(REAL(lod)   , 0, width * i_depth * sizeof(double));

  PROTECT(ans    = allocVector(VECSXP, 3)); 

  in_names = getAttrib(v, R_DimNamesSymbol);
  cnames = GetColNames(in_names);

  PROTECT(snp_names = allocVector(STRSXP, width + 1));
  idx_j = 0;
  for (idx_i = ii -1 ; idx_i <= jj -1 ; idx_i ++)
    {
      SET_STRING_ELT(snp_names, idx_j, STRING_ELT(cnames, idx_i));
      idx_j++;
    }
  
  for(idx_j = 0; idx_j < i_depth; idx_j++)
    {
      for (idx_i = ii - 1 ; idx_i <= jj - 2 - idx_j; idx_i++)
	{	  
	  geno_cptr res = get_geno_count(RAW(v) + idx_i * rows, 
					 RAW(v) + (idx_i + idx_j + 1) * rows , 
					 rows);
	  int offset = idx_i +1 - ii;
	  REAL(dprime)[offset + idx_j * width] = res->dprime;
	  if (need_signed_r)
	    {
	      if (res->rsq2 > 0)
		{
		  REAL(rmisc)[offset + idx_j * width] = res->sign_of_r * sqrt(res->rsq2);
		}
	      else
		{
		  REAL(rmisc)[offset + idx_j * width] = -2;
		}
	    }
	  else
	    {
	      REAL(rmisc)[offset + idx_j * width] = res->rsq2;
	    }
	  REAL(lod)[offset + idx_j * width] = res->lod;
	  free(res->expt);
	  free(res);
	}
    }

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

  PROTECT(class_name = allocVector(STRSXP, 1));
  SET_STRING_ELT(class_name, 0, mkChar("snp.dprime"));
  classgets(ans, class_name);  
  setAttrib(ans, install("snp.names"), snp_names);

  UNPROTECT(13);
  Rprintf("... Done\n");  
  return ans;
}
