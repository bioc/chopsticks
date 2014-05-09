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
#include <Rversion.h>

#if (R_VERSION >= R_Version(2,3,0))
#define HAVE_FLUSHCONSOLE
#endif

#include "zlib.h"

#define HEADER_BEGINNING "rs# SNPalleles chrom pos strand genome_build center protLSID assayLSID panelLSID QC_code "
#define HEADER_BEGINNING_ALT "rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode "
#define SKIP_FIELDS 11 /* number of spaces in the string above */

#define LINE_BUF_SIZE 4096
#define MAX_SNP_LEN 64
#define MAX_SNP 1000000

typedef struct linecontent {
  struct linecontent *next;
  char name[MAX_SNP_LEN];
  char alleles[64]; /* These sizes of these fields are chosen to be big and align position to 8 byte */ 
  char chrom[8];
  char assignment[4];
  char strand[4];  
  unsigned int pos;     /* this one should be aligned to 8 byte boundary */
  char *snps;
} *linecontent_ptr;

typedef struct _char_trail {
  char base;
  int offset;
} char_trail;

static char_trail bases[] = {
  {'A', 0},
  {'C', 2},
  {'G', 6},
  {'T', 19}
};


static int valid_gtype(char a) {
  switch (a) {
  case 'A':
  case 'C':
  case 'T':
  case 'G':
    return 1;
    break;
  default: 
    return 0;
    break;
  }
  return 0;
}

static int expected_gtype(char a) {
  if (a == 'N') {
    return 1;
  } else {
    return valid_gtype(a);
  }
}

#define GMASK_INVALID 0x40
#define GMASK_1 0x01
#define GMASK_2 0x02

static void contentlist_destroy(linecontent_ptr start) {
  linecontent_ptr cur_content_ptr = NULL;
  for(cur_content_ptr = start; cur_content_ptr ; /* advancing pointer inside */ ) {
    linecontent_ptr copy = cur_content_ptr; /* make a copy before moving on */
    cur_content_ptr = cur_content_ptr->next;
    free(copy->snps);
    free(copy);
  }
}

#define MAKEUP_GTYPE '.'

static int guess_genotypes(char *line_ptr, char *alleles, int verbose) {
  char bin[26];
  int result = 0;
  int seen_genotype = 0;
  int batch_genotype =0;
  char seen[5];
  int a_offset = strlen(alleles) - 2;
  int i = 0;
  memset(bin, 0x00, 26);
  memset(seen, 0x00, 5);
  char *seenbuf = seen;
  while (*line_ptr) {
    char a = *line_ptr++;
    if (valid_gtype(a)) {
      bin[ a - 'A' ] = 0x01;
    } else if ((a != ' ') && (a != 'N') && (a != '\n')) {
      /* see anything strange, bail out */
      return -1; 
    }
  }

  for (i = 0; i < 4 ; i++) {
    if (bin[bases[i].offset]) {
      char base = bases[i].base;
      seen_genotype++;
      result <<= 8;
      result |= base;
      sprintf(seenbuf, "%c", base);
      seenbuf++;
    }
  }

  batch_genotype = seen_genotype;

  if (seen_genotype == 0) {
    result = ((MAKEUP_GTYPE << 8)  | (MAKEUP_GTYPE));
    seen_genotype +=2;
  } else if (seen_genotype == 1) {
    if (result == 'A') {
      result <<= 8;
      result |= MAKEUP_GTYPE;
      seen_genotype ++;
    } else if (result == 'T') {
      result |= (MAKEUP_GTYPE << 8);
      seen_genotype ++;
    } else if (!strncmp(alleles, "-/A/", 4) || !strncmp(alleles, "A/", 2) || !strncmp(alleles, "(LARGEINSERTION)/-/A/", 21)) {
      if (result != 'A') {
	result |= ('A' << 8);
      } else {
	result <<=8;
	result |= MAKEUP_GTYPE; /* same as first case */
      }
      seen_genotype ++;
    } else if (!strncmp(alleles, "-/C/", 4) || !strncmp(alleles, "C/", 2) || !strncmp(alleles, "(LARGEINSERTION)/-/G/", 21)) {
      if (result != 'C') {
	result |= ('C' << 8);
      } else {
	result <<=8;
	result |= MAKEUP_GTYPE;
      }
      seen_genotype ++;
    } else if (!strncmp(alleles, "-/G/", 4) || !strncmp(alleles, "G/", 2) || !strncmp(alleles, "(LARGEINSERTION)/-/G/", 21)) {
      if (result != 'G') {
	result |= ('G' << 8);
      } else {
	result <<=8;
	result |= MAKEUP_GTYPE;
      }
      seen_genotype ++;
    } else if (!strncmp(alleles + a_offset, "/C", 2)) {
      if (result != 'C') {
	result <<= 8;
	result |= 'C';
      } else {
	result |= (MAKEUP_GTYPE << 8);
      }
      seen_genotype ++;
    } else if (!strncmp(alleles + a_offset, "/G", 2)) {
      if (result != 'G') {
	result <<= 8;
	result |= 'G';
      } else {
	result |= (MAKEUP_GTYPE << 8);
      }
      seen_genotype ++;
    } else if (!strncmp(alleles + a_offset, "/T", 2)) {
      if (result != 'T') {
	result <<= 8;
	result |= 'T';
      } else {
	result |= (MAKEUP_GTYPE << 8);
      }
      seen_genotype ++;
    }    
  }
  
  if (seen_genotype == 2) {
    if ((verbose) && (batch_genotype != seen_genotype)) {
      Rprintf("Reducing dbSNPs %s\tagainst seen %s", alleles, seen);      
      Rprintf("\t: %c/%c\n", (result >> 8), (result & 0xff)); 
    }
    return result;
  } else {
    return -1;
  }
}

SEXP read_hapmap_data(SEXP downloaded_file, SEXP sexp_verbose) {
  SEXP dimnames = R_NilValue;
  SEXP ans = R_NilValue, ans_name = R_NilValue, snp_matrix_ans = R_NilValue;
  SEXP snp_names = R_NilValue;
  SEXP sample_names = R_NilValue;

  SEXP supp_df = R_NilValue, supp_colnames = R_NilValue, supp_snp_names = R_NilValue, supp_alleles = R_NilValue,
    supp_chrom = R_NilValue, supp_assignment = R_NilValue, supp_strand = R_NilValue,
    supp_pos = R_NilValue;

  int verbose = 0;
  const char *filename = NULL;
  char current_line[LINE_BUF_SIZE];
  char *line_ptr = NULL;
  gzFile gz_fp;
  int snp_count = 0; 
  int i_snps = 0;
  int read_failed = 1; /* = "yes", default */
  int n_samples = 1;
  char gmask[96]; /* only need to up to 'Z', which is 91 */
  struct linecontent listhead;
  linecontent_ptr cur_content_ptr = &listhead;
  int header_skip = strlen(HEADER_BEGINNING);
  int i;

  PROTECT(downloaded_file = coerceVector(downloaded_file, STRSXP));
  PROTECT(sexp_verbose = coerceVector(sexp_verbose, LGLSXP));
  if(TYPEOF(downloaded_file) != STRSXP)
    Rprintf(" input filename wrong type\n");
  if(TYPEOF(sexp_verbose) != LGLSXP)
    Rprintf(" verbose argument not of logical type\n");
  
  verbose = LOGICAL(sexp_verbose)[0];

  filename = CHAR(STRING_ELT(downloaded_file, 0));
  gz_fp = gzopen(filename, "rb");
  
  if (!gz_fp) {
    Rprintf("opening downloaded file failed.\n");
    UNPROTECT(2); 
    return R_NilValue;
  }

  do {
    /* do-while loops does the body at least once */
    if (gzgets(gz_fp, current_line, LINE_BUF_SIZE) == Z_NULL) {
      Rprintf("read header failed\n");
      UNPROTECT(2); 
      return R_NilValue;
    }
    
    if (strlen(current_line) +1 >= LINE_BUF_SIZE) {
      Rprintf("Header unexpectedly long\n");
      Rprintf("Please recompile with large LINE_BUFFER_SIZE\n");
      UNPROTECT(2); 
      return R_NilValue;
    }
    /* skip over initial comments starting with '#' */
  } while (current_line[0] == '#');

  if (strncmp(current_line, HEADER_BEGINNING, header_skip)) {
    /* doesn't match normal header, try alternative header */
    header_skip = strlen(HEADER_BEGINNING_ALT);
    if (strncmp(current_line, HEADER_BEGINNING_ALT, header_skip)) {
      Rprintf("Unexpected header: %s\n", current_line);
      UNPROTECT(2);
      return R_NilValue;
    }
  }

  /* n_samples starts at 1, since no space is 1 sample */
  for (line_ptr = current_line + header_skip; (*line_ptr) ; line_ptr++) {
    if (*line_ptr == ' ') { 
      n_samples++;
      *line_ptr = '\0'; 
    }
  }
  line_ptr--;
  if (*line_ptr == '\n') {
    /* debug Rprintf */
    /* Rprintf("removing end of line in header\n"); */
    *line_ptr = '\0';
  }

  Rprintf("Reading %d samples\n", n_samples);

  PROTECT(sample_names = allocVector(STRSXP, n_samples));
  
  line_ptr = current_line + header_skip;
  for (i = 0; i < n_samples; i ++) {    
    SET_STRING_ELT(sample_names, i, mkChar(line_ptr));
    
    while(*line_ptr) line_ptr++; /* advance to the next null */
    line_ptr++; /* next non-null */
  }

  while (1) {
    char *snps = NULL;
    char *alleles = NULL;
    int get1 = 0; 
    /* read until breaking out because of eof or error */
    if (gzeof(gz_fp)) {
      gzclose(gz_fp);
      Rprintf("current line [%d] : %.20s...\n", snp_count - 1, current_line);
      Rprintf("EOF reached after %d snps\n", snp_count);
      read_failed = 0; /* read_failed = "no" */
      break;
    }
    
    /* gzeof() is not reliable for non-compressed files,
       trying alternative method for testing end of file*/
    if ((get1 = gzgetc(gz_fp)) != -1) {
      /* get1 was successful, put it back */
      if (gzungetc(get1, gz_fp) != get1) {
	Rprintf("read error (gzungetc) after %d snps\n", snp_count);
	break;
      }
    } else {
      Rprintf("last line [%d] : %.20s...\n", snp_count - 1, current_line);
      if (gzeof(gz_fp)) {
        gzclose(gz_fp);
	Rprintf("EOF reached after %d snps\n", snp_count);
      } else {
	Rprintf("No more to read after %d snps\n", snp_count);	
      }
      read_failed = 0; /* read_failed = "no" */
      break;
    }

    if (gzgets(gz_fp, current_line, LINE_BUF_SIZE) == Z_NULL) {
      Rprintf("read error (gzgets) after %d snps\n", snp_count);
      break;
    }
    
    if (strlen(current_line) +1 >= LINE_BUF_SIZE) {
      Rprintf("Unexpected long lines after %d snps\n", snp_count);
      Rprintf("Please recompile with large LINE_BUFFER_SIZE\n");
      break;
    }

    cur_content_ptr->next = (linecontent_ptr)calloc(sizeof(struct linecontent), 1);
    cur_content_ptr = cur_content_ptr->next;
    cur_content_ptr->snps = (char *)calloc(sizeof(char), n_samples);
    alleles = cur_content_ptr->alleles;
    snps = cur_content_ptr->snps;
    int i = sscanf(current_line, "%s %s %5s %u %1s", cur_content_ptr->name, 
		   alleles,
		   cur_content_ptr->chrom,
		   &(cur_content_ptr->pos),
		   cur_content_ptr->strand
		   );
    if ((i != 5) || (i == EOF) || 
	(strlen(cur_content_ptr->chrom) > 5) || 
	(strlen(cur_content_ptr->strand) > 1) || 
        ( (cur_content_ptr->strand[0] != '+') && (cur_content_ptr->strand[0] != '-') )
	){
      Rprintf("mailformed input line: %s", current_line);
      break;
    }

    memset(gmask, GMASK_INVALID, 96);

    /* Advance to the snp data position */
    line_ptr = current_line;
    int skip_count;
    for (skip_count = 0; (*line_ptr) && (skip_count < SKIP_FIELDS) ; line_ptr++) {
      if (*line_ptr == ' ') skip_count++;
    }

    if ((alleles[1] != '/') || (alleles[0] == alleles[2]) || !valid_gtype(alleles[0]) || !valid_gtype(alleles[2]) 
	|| (strlen(alleles) != 3)) {
      int guesses = guess_genotypes(line_ptr, alleles, verbose);
      if (guesses < 0) {	
	Rprintf("snp %s allele field unexpected %s\n", cur_content_ptr->name, alleles);
	Rprintf(" calls: %s", line_ptr); /* no need to do \n, since it is included */ 
	break;
      } else {
	unsigned char a = guesses >> 8;
	unsigned char b = guesses & 0xff;
	gmask[a] = GMASK_1 ;
	gmask[b] = GMASK_2 ;
	sprintf(cur_content_ptr->assignment, "%c/%c", a, b);
	
      }
    } else {
      gmask[(int) alleles[0]] = GMASK_1 ;
      gmask[(int) alleles[2]] = GMASK_2 ;
      strcpy(cur_content_ptr->assignment, alleles); /* we know at this point strlen(alleles) is 3 */
    }

    while (*line_ptr) {
      char c;
      unsigned char a = *line_ptr++;
      unsigned char b = *line_ptr++;
      if (!expected_gtype(a) || !expected_gtype(b)) {
	Rprintf("unexpected input buffer %c%c%20s\n", a, b, line_ptr);
	break;
      }
      c =  gmask[a] | gmask[b];
      switch (c) {
      case 1: /* 1 | 1 */
	*snps = 1;
	break;
      case 2: /* 2 | 2 */
	*snps = 3;
	break;
      case 3: /* 1 | 2 or 2 | 1 */
	*snps = 2;
	break;
      default:
	*snps = 0;
	break;
      }
      snps ++;
      c = *line_ptr++; 
      if ((snps - cur_content_ptr->snps) == n_samples) {
	while (c == ' ') {
	  c = *line_ptr++;
	}
	if ((c != '\n')) {	  
	  Rprintf(" unexpected end of line %02X\n", c);
	}
      } else {
	if (c != ' ') {
	  Rprintf(" unexpected middle of line %02X\n", c);
	}
      }
    } 

    if (!(snp_count % 20000)) {
      Rprintf("current line [%d] : %.20s...\n", snp_count, current_line);
      /* windows and aqua needs the flush */
#ifdef HAVE_FLUSHCONSOLE
      R_FlushConsole();
#endif
    }

    snp_count++;
  }

  /* read_failed defaults to yes, unless EOF reached */
  if (read_failed) {
    /* something went wrong, tidy up and return nil */
    contentlist_destroy(listhead.next);
    UNPROTECT(3);
    return R_NilValue;
  }

  PROTECT(snp_matrix_ans = allocMatrix(RAWSXP, n_samples, snp_count));
  PROTECT(snp_names = allocVector(STRSXP,  snp_count));

  /* support data  - snp_names is duplicated later */
  PROTECT(supp_df         = allocVector(VECSXP, 5));
  PROTECT(supp_colnames   = allocVector(STRSXP, 5));
  PROTECT(supp_alleles    = allocVector(STRSXP,  snp_count));
  PROTECT(supp_assignment = allocVector(STRSXP,  snp_count));
  PROTECT(supp_chrom      = allocVector(STRSXP,  snp_count));
  PROTECT(supp_strand     = allocVector(STRSXP,  snp_count));
  PROTECT(supp_pos        = allocVector(INTSXP,  snp_count));
  
  i_snps = 0;
  for(cur_content_ptr = listhead.next; cur_content_ptr ; cur_content_ptr = cur_content_ptr->next) {
    char *snps = cur_content_ptr->snps;
    memcpy((RAW(snp_matrix_ans) + n_samples * i_snps) , snps, n_samples);
    SET_STRING_ELT(snp_names, i_snps, mkChar(cur_content_ptr->name));
    /* copying the support data */
    SET_STRING_ELT(supp_alleles,    i_snps, mkChar(cur_content_ptr->alleles));
    SET_STRING_ELT(supp_assignment, i_snps, mkChar(cur_content_ptr->assignment));
    SET_STRING_ELT(supp_chrom,      i_snps, mkChar(cur_content_ptr->chrom));
    SET_STRING_ELT(supp_strand,     i_snps, mkChar(cur_content_ptr->strand));
    INTEGER(supp_pos)[i_snps] = cur_content_ptr->pos;
    /* done copying support data */
    i_snps++;
  }

  contentlist_destroy(listhead.next);
  
  /* make another copy of the snp names for the support data frame */
  PROTECT(supp_snp_names = duplicate(snp_names));
  /* constructing the support data frame */
  SET_VECTOR_ELT(supp_df, 0, supp_alleles);
  SET_VECTOR_ELT(supp_df, 1, supp_assignment);
  SET_VECTOR_ELT(supp_df, 2, supp_chrom);
  SET_VECTOR_ELT(supp_df, 3, supp_pos);
  SET_VECTOR_ELT(supp_df, 4, supp_strand);
  SET_STRING_ELT(supp_colnames, 0 , mkChar("dbSNPalleles"));
  SET_STRING_ELT(supp_colnames, 1 , mkChar("Assignment"));
  SET_STRING_ELT(supp_colnames, 2 , mkChar("Chromosome"));
  SET_STRING_ELT(supp_colnames, 3 , mkChar("Position"));
  SET_STRING_ELT(supp_colnames, 4 , mkChar("Strand"));
  setAttrib(supp_df, R_NamesSymbol, supp_colnames);
  setAttrib(supp_df, R_RowNamesSymbol, supp_snp_names);
  setAttrib(supp_df, R_ClassSymbol, mkString("data.frame")); /* data.frame is not built-in */

  /* construct the rest of the snp.data */
  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 1, snp_names);    /* colnames */
  SET_VECTOR_ELT(dimnames, 0, sample_names); /* rownames */
  setAttrib(snp_matrix_ans, R_DimNamesSymbol, dimnames);
  setAttrib(snp_matrix_ans, R_ClassSymbol, mkString("snp.matrix"));
  SET_S4_OBJECT(snp_matrix_ans);

  /* now make the result list */
  PROTECT(ans      = allocVector(VECSXP, 2));
  PROTECT(ans_name = allocVector(STRSXP, 2));
  SET_VECTOR_ELT(ans, 0, snp_matrix_ans);
  SET_VECTOR_ELT(ans, 1, supp_df);
  SET_STRING_ELT(ans_name, 0, mkChar("snp.data"));
  SET_STRING_ELT(ans_name, 1, mkChar("snp.support"));
  setAttrib(ans, R_NamesSymbol, ans_name);

  UNPROTECT(16);
  return ans;
}
