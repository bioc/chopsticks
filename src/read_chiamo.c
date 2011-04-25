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

#include "zlib.h"

#include "read_tokens.h"

/* headerless, and spaces (not tabs):
SNP_A-2186786 rs3097718 341390 A G 0.9988 0.0012 0 1 0 0
SNP_A-4298665 --- 344301 A G 0 0 1 0.0019 0.9930 0 0 0
*/

#define MAX_ID_LENGTH 128

typedef struct linecontent {
  struct linecontent *next;
  char vendor_id[MAX_ID_LENGTH];
  char public_id[MAX_ID_LENGTH];
  int pos;
  char alleleA;
  char alleleB;
  char *snps;
} *linecontent_ptr;

static void contentlist_destroy(linecontent_ptr start) {
  linecontent_ptr cur_content_ptr = NULL;
  for(cur_content_ptr = start; cur_content_ptr ; /* advancing pointer inside */ ) {
    linecontent_ptr copy = cur_content_ptr; /* make a copy before moving on */
    cur_content_ptr = cur_content_ptr->next;
    free(copy->snps);
    free(copy);
  }
}

static char genotype_from_prob(float *array, float cutoff)
{
  char calls[3] = {0x01, 0x02, 0x03};
#define idx_max3(a,b,c) ( (a) < (b) ? ( (b) < (c) ? 2 : 1 )	\
			  : ( (a) < (c) ? 2 : 0 ) )
  int idx = idx_max3(array[0], array[1], array[2]);

  if (array[idx] > cutoff) {
    return calls[idx];
  } else {
  return 0x00; /* no call */
  }
}

SEXP read_chiamo(SEXP mcmc_file, SEXP sample_list, SEXP threshold) {
  int n_samples = 0;
  int n_snps = 0;
  float f_threshold = 0.9;
  int read_failed = 1;

  struct linecontent listhead;
  linecontent_ptr cur_content_ptr = &listhead;

  /* check inputs */
  if(TYPEOF(mcmc_file) != STRSXP)
    error("input file name wrong type\n");
  if(TYPEOF(sample_list) != STRSXP)
    error("input sample list wrong type\n");
  if(TYPEOF(threshold) != REALSXP)
    error("input threshold wrong type\n");

  for(int i_file = 0; i_file < LENGTH(mcmc_file); i_file++) {
  const char *filename = CHAR(STRING_ELT(mcmc_file, i_file));
  gzFile gz_fp = gzopen(filename, "rb");

  f_threshold = REAL(threshold)[0];

  if (!gz_fp) {
    Rprintf("opening input file %s failed.\n", filename);
    return R_NilValue;
  }

  /* work out line length and number of samples - don't really
     trust the length of the input sample_list */
  /* read 4k blocks until the first line ending, then rewind,
     to count number of fields, and line length */
#define BUF_SIZE 4192
  char buffer[BUF_SIZE+1]; /* add one so we'll can null terminate later */
  int read_done = 0; /* not done */
  int n_buf = 0; /* used for line buffer size below */
  int n_read =0;
  int n_fields = 0;
  int is_blank = 1; /* start with white space, so the first non-white counts as one */
  while((n_read = gzread(gz_fp, buffer, BUF_SIZE)) != 0) {
    char *b_ptr = buffer;
    n_buf++; /* increment the size of the line buffer needed later */
    /* terminate explicitly as the buffer may contain garbage from last
       read after the end of current read */
    buffer[n_read] = '\0';
    while(*b_ptr) {
      char a = *b_ptr;
      if (a == '\n') {
        read_done = 1;
        break;
      }
      if ((a == ' ') || (a == '\t') || (a == '\f') || (a == '\r') || (a == '\v')) {
        is_blank = 1;
      } else {
        if (is_blank) {
          /* counting blank to non-blank transitions */
          n_fields++;
        }
        is_blank = 0; /* is non-blank */
      }
      b_ptr++;
    }
    if (read_done)
      break;
  }

  if (!(read_done)) {
    Rprintf("Had read error after scanning for %i fields in header\n", n_fields);
    return R_NilValue;
  } else {
    int len_list = LENGTH(sample_list);
    n_samples = (n_fields - 5)/3 ;
    Rprintf("Found %i samples in the first line.\n", n_samples);
    if (len_list != n_samples) {
      Rprintf("Input sample list contains %i samples, not equal to %i from chiamo output.\n",
	      len_list, n_samples);
      return R_NilValue;
    }
  }
  gzrewind(gz_fp);

  /* chiamo inputs are widely varying in length, so we use the
     field length, and hopefully the decimal places are not too many*/
#define MAX(a,b) ((a) > (b) ? (a):(b))
  int line_buffer_size = MAX(n_buf * BUF_SIZE, n_fields * 10);
  char *current_line = malloc(line_buffer_size);

  read_failed = 1;
  while (1) {
    int get1 = 0;
    int i = 0;
    char *line_ptr=current_line;
    /* read until breaking out because of eof or error */
    if (gzeof(gz_fp)) {
      gzclose(gz_fp);
      Rprintf("current line [%d] : %.20s...\n", n_snps, current_line);
      Rprintf("EOF reached after %d samples\n", n_snps);
      read_failed = 0; /* read_failed = "no" */
      break;
    }

    /* gzeof() is not reliable for non-compressed files,
       trying alternative method for testing end of file*/
    if ((get1 = gzgetc(gz_fp)) != -1) {
      /* get1 was successful, put it back */
      if (gzungetc(get1, gz_fp) != get1) {
        Rprintf("read error (gzungetc) after %d snps\n", n_snps);
        break;
      }
    } else {
      if (gzeof(gz_fp)) {
        gzclose(gz_fp);
        Rprintf("last line [%d] : %.20s...\n", n_snps, current_line);
        Rprintf("EOF reached after %d snps\n", n_snps);
        read_failed = 0; /* read_failed = "no" */
        break;
      }
    }
    /* normal eof detection ends - below are either processing or errors */

    if (gzgets(gz_fp, current_line, line_buffer_size) == Z_NULL) {
      Rprintf("read error (gzgets) after %d snps\n", n_snps);
      break;
    }
    n_snps++;

    cur_content_ptr->next = (linecontent_ptr)calloc(sizeof(struct linecontent), 1);
    cur_content_ptr = cur_content_ptr->next;
    cur_content_ptr->snps = (char *)calloc(sizeof(char), n_samples);
    char *snps = cur_content_ptr->snps;

    sscanf(line_ptr, "%s %s %d %c %c",
	   cur_content_ptr->vendor_id,
	   cur_content_ptr->public_id,
	   &(cur_content_ptr->pos),
	   &(cur_content_ptr->alleleA),
	   &(cur_content_ptr->alleleB));
    skip_5_tokens(line_ptr);

    for(i = 0; i < n_samples ; i++) {
      float f_num[3];
      sscanf(line_ptr, "%f %f %f", f_num, (f_num+1), (f_num+2));
      *snps = genotype_from_prob(f_num, f_threshold);
      snps++;
      skip_3_tokens(line_ptr);
    }

  } /* while (1) */

  free(current_line);

  if (read_failed) {
    break; /* break out of the for loop */
  }
  }

  if(read_failed) {
    /* no message, as any relevant message should be above
       where the error occurs */
    return R_NilValue;
  } else {
    Rprintf("Read %i samples with %i snps from input, now converting...\n", n_samples, n_snps);
  }

  SEXP vendor_id = R_NilValue, public_id = R_NilValue, pos = R_NilValue;
  SEXP alleleA = R_NilValue, alleleB = R_NilValue;
  SEXP snp_data = R_NilValue;

  /* sample support data */
  PROTECT(vendor_id = allocVector(STRSXP, n_snps));
  PROTECT(public_id = allocVector(STRSXP, n_snps));
  PROTECT(pos       = allocVector(INTSXP, n_snps));
  PROTECT(alleleA   = allocVector(STRSXP, n_snps));
  PROTECT(alleleB   = allocVector(STRSXP, n_snps));

  PROTECT(snp_data  = allocMatrix(RAWSXP, n_samples, n_snps));
  int protected = 6;
  int i_snps = 0;

  for(cur_content_ptr = listhead.next; cur_content_ptr ; cur_content_ptr = cur_content_ptr->next) {
    char *snps = cur_content_ptr->snps;
    memcpy((RAW(snp_data) + n_samples * i_snps), snps, n_samples);
    char aA[2] = {cur_content_ptr->alleleA, 0x00};
    char aB[2] = {cur_content_ptr->alleleB, 0x00};
    SET_STRING_ELT(vendor_id, i_snps, mkChar(cur_content_ptr->vendor_id));
    SET_STRING_ELT(public_id, i_snps, mkChar(cur_content_ptr->public_id));
    INTEGER(pos)[i_snps] = cur_content_ptr->pos;
    SET_STRING_ELT(alleleA,   i_snps, mkChar(aA));
    SET_STRING_ELT(alleleB,   i_snps, mkChar(aB));
    i_snps++;
  }
  contentlist_destroy(listhead.next);

  /* the rest of the snp.matrix */
  SEXP snp_data_class = R_NilValue, snp_data_dimnames = R_NilValue;
  PROTECT(snp_data_dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(snp_data_dimnames, 0, duplicate(sample_list));
  SET_VECTOR_ELT(snp_data_dimnames, 1, duplicate(vendor_id));
  setAttrib(snp_data, R_DimNamesSymbol, snp_data_dimnames);
  PROTECT(snp_data_class = allocVector(STRSXP, 1));
  SET_STRING_ELT(snp_data_class, 0, mkChar("snp.matrix"));
  classgets(snp_data, snp_data_class);
  SET_S4_OBJECT(snp_data);
  protected +=2;

  /* the rest of the snps support data frame */
  SEXP snps_support_df, snps_support_class, snps_support_df_names;
  PROTECT(snps_support_df_names = allocVector(STRSXP, 5));
  PROTECT(snps_support_df = allocVector(VECSXP, 5));
  SET_VECTOR_ELT(snps_support_df, 0, vendor_id);
  SET_VECTOR_ELT(snps_support_df, 1, public_id);
  SET_VECTOR_ELT(snps_support_df, 2, pos);
  SET_VECTOR_ELT(snps_support_df, 3, alleleA);
  SET_VECTOR_ELT(snps_support_df, 4, alleleB);
  SET_STRING_ELT(snps_support_df_names, 0, mkChar("Vendor.ID"));
  SET_STRING_ELT(snps_support_df_names, 1, mkChar("Public.ID"));
  SET_STRING_ELT(snps_support_df_names, 2, mkChar("Position"));
  SET_STRING_ELT(snps_support_df_names, 3, mkChar("AlleleA"));
  SET_STRING_ELT(snps_support_df_names, 4, mkChar("AlleleB"));

  setAttrib(snps_support_df, R_RowNamesSymbol, duplicate(vendor_id));
  setAttrib(snps_support_df, R_NamesSymbol, snps_support_df_names);
  PROTECT(snps_support_class = allocVector(STRSXP, 1));
  SET_STRING_ELT(snps_support_class, 0, mkChar("data.frame"));
  classgets(snps_support_df, snps_support_class);
  protected += 3;

  /* final result */
  SEXP ans = R_NilValue, ans_names = R_NilValue;
  PROTECT(ans = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(ans, 0, snp_data);
  SET_VECTOR_ELT(ans, 1, snps_support_df);
  PROTECT(ans_names = allocVector(STRSXP, 2));
  SET_STRING_ELT(ans_names, 0, mkChar("snp.data"));
  SET_STRING_ELT(ans_names, 1, mkChar("snp.support"));
  setAttrib(ans, R_NamesSymbol, ans_names);
  protected += 2;

  UNPROTECT(protected);
  return ans;
}
