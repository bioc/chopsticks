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

#include <R.h>
#include <Rinternals.h>
#include "zlib.h"
#include "hash_index.h"

#include "read_tokens.h"

/* taken from R/src/main/bind.c */
static void SetRowNames(SEXP dimnames, SEXP x)
{
  if (TYPEOF(dimnames) == VECSXP)
    SET_VECTOR_ELT(dimnames, 0, x);
  else if (TYPEOF(dimnames) == LISTSXP)
    SETCAR(dimnames, x);
}

static void SetColNames(SEXP dimnames, SEXP x)
{
  if (TYPEOF(dimnames) == VECSXP)
    SET_VECTOR_ELT(dimnames, 1, x);
  else if (TYPEOF(dimnames) == LISTSXP)
    SETCADR(dimnames, x);
}

SEXP read_signals(SEXP signalfile, SEXP snp_list) {
  SEXP ans = R_NilValue, ans_name =  R_NilValue, sample_names = R_NilValue;
  const char *filename = NULL;
  int buffersize = 200000; /* T1D - 2000 samples, requires 72428 */ 
  int i =0;
  
  char *line_buffer = (char *)malloc(buffersize);
  char *line_ptr = line_buffer;
  
  index_db snp_db = NULL;
  
  if(TYPEOF(signalfile) != STRSXP)
    Rprintf(" input filename wrong type\n");
  if(TYPEOF(snp_list) != STRSXP)
    Rprintf(" input snp.list wrong type\n");
  
  filename = CHAR(STRING_ELT(signalfile, 0));
  gzFile ginfile = gzopen(filename, "rb");
  if (!ginfile) {
    Rprintf("Cannot read %s\n", filename);
    return R_NilValue;
  }
  Rprintf("Reading %s ...\nCan take a while...\n", filename);
  
  int array_length = LENGTH(snp_list);

  PROTECT(ans    = allocVector(VECSXP, array_length));
  PROTECT(ans_name = duplicate(snp_list));
  setAttrib(ans, R_NamesSymbol, ans_name);

  snp_db = index_create(array_length);  

  for(i = 0 ; i < array_length ; i++) {
    index_insert(snp_db, CHAR(STRING_ELT(snp_list, i)), i);
    /* fills with NULLs first, just in case we can't find it later */
    SET_VECTOR_ELT(ans, i, R_NilValue);
  }
  gzgets(ginfile, line_buffer, buffersize);
  
  int space_count = 0;
  line_ptr = line_buffer;
  while (*line_ptr) {
    if ((*line_ptr == ' ') || (*line_ptr == '\t')) {      
      space_count++;
    }
    line_ptr++;
  }

  int pair_count = (space_count - 4)/2;
  int array_count_down = array_length;

  line_ptr = line_buffer;
  skip_5_tokens(line_ptr);

  PROTECT(sample_names = allocVector(STRSXP, pair_count));
  for (i = 0; i< pair_count ; i++) {
    char *ptr_start = line_ptr;
    goto_next_token(line_ptr); *(line_ptr-3) = '\0'; /* chop off the last two characters */
    SET_STRING_ELT(sample_names, i, mkChar(ptr_start));
    goto_next_token(line_ptr);
  }
  
  int line_read_count = 1;
  while(1) {
    int get1 = 0;
    int i = 0;
    /* found all, no need to read further */
    if (!array_count_down) {
      break;
    }
    
    if (gzeof(ginfile)) {
      break;
    }
    if (!(line_read_count % 200)) {
      Rprintf("Reading line %i\r",line_read_count);
    }
    
    /* zlib 1.1 don't have VERNUM nor gzungetc() */
#ifdef ZLIB_VERNUM
    /* gzeof() is sometimes unreliable until the next read */
    if ((get1 = gzgetc(ginfile)) != -1) {
      /* get1 was successful, put it back */
      if (gzungetc(get1, ginfile) != get1) {
	Rprintf("Unexpected file system error\n");
        break;
      }
    } else {
      if (gzeof(ginfile)) {
        break;
      }
    }
#endif

    gzgets(ginfile, line_buffer, buffersize);
    line_read_count ++;

    line_ptr = line_buffer;
    goto_next_token(line_ptr); *(line_ptr-1) = '\0';
    int i_part = index_lookup(snp_db, line_buffer);
    if(i_part <0) {
      char *ptr_start = line_ptr;
      goto_next_token(line_ptr); *(line_ptr-1) = '\0';
      i_part = index_lookup(snp_db, ptr_start);
    }

    if(i_part <0) {
      /* neither first or 2nd field matched */
      continue;
    }

    /* if we get here we found something */
    array_count_down--;
    SEXP AB = R_NilValue, part =  R_NilValue, part_name = R_NilValue;
    SEXP dims = R_NilValue;
    /* allocMatrix(READSXP...) seems to be sane enough
       to automatically filled with NA_REAL so that when sscanf() below 
       fails, that's what one luckily gets, in the Illumina platform...
    */
    PROTECT(AB = allocMatrix(REALSXP, pair_count, 2));
    double *A = REAL(AB);
    double *B = A + pair_count;

    line_ptr = line_buffer;
    skip_5_tokens(line_ptr);

    for (i = 0; i< pair_count ; i++) {
      char *ptr_start = line_ptr;
      goto_next_token(line_ptr); *(line_ptr-1) = '\0';
      sscanf(ptr_start, "%lf", A + i);
      ptr_start = line_ptr;
      goto_next_token(line_ptr); *(line_ptr-1) = '\0';
      sscanf(ptr_start, "%lf", B + i);
    }
    
    PROTECT(part = allocVector(VECSXP, 2));
    PROTECT(part_name = allocVector(STRSXP, 2));
    PROTECT(dims = allocVector(INTSXP, 2));
    INTEGER(dims)[0] = pair_count;
    INTEGER(dims)[1] = 2;
    SET_STRING_ELT(part_name, 0, mkChar("A"));
    SET_STRING_ELT(part_name, 1, mkChar("B"));
    SetRowNames(part, sample_names);
    SetColNames(part, part_name);
    setAttrib(AB, R_DimSymbol, dims);
    /* set DimSymbol has the side-effect of erasing DimNames */
    setAttrib(AB, R_DimNamesSymbol, part); 
    SET_VECTOR_ELT(ans, i_part, AB);
  }

  if(array_count_down) {
    warning("EOF reached but %d SNPs was not found\n", array_count_down);
  } else {
    /* extra linefeed to avoid cobbering from early informative input */
    Rprintf("\n...Done\n");
  }

  index_destroy(snp_db);
  gzclose(ginfile);

  /* TODO: do we need to do something about the original 
     copy of sample_names or do nothing and let it gets 
     garbage-collected eventually?
   */

  UNPROTECT(3 + (array_length - array_count_down) *4);
  return ans;
}
