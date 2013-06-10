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


/*
   The interface to insnp_long_unsorted() is deliberately similiar
   to insnp_long().
*/


/* Read genotype data in long format

filename    Name of input text file (4 fields per row)

nchip       Number of chips
chip_id     Array of chip identifier strings
nsnp        Number of snps
snp_id      Array of snp identifier strings

ifX         Indicator for X chromosome
female      If X chromosome, indicator array for female sex

ifconf      If +ve, high score on confidence score is good, if -ve high
            score is bad. If zero no score is present in input.
threshold   Quality threshold for inclusion of data

progress    1 if a progess report to stdout, else 0.

gtypes      Array of char's, length nchip*nsnp, to hold genotype data
counts      Various counts
iferror     Error indicator

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "zlib.h"

#include "hash_index.h"

#ifndef WITHOUT_R
#include "R.h"
#include "Rversion.h"
#define message_out Rprintf
#define FLUSH R_FlushConsole
#endif

#define LINE_BUFFER_SIZE 1024
#define MAX_GT 16
/* anything but 0 to 3 - cannot use 0xFF nor -1 - char can be signed or unsigned */
#define CONFLICT_GTYPE 0x7E
#define INVALID_GTYPE 0x7F

#if (R_VERSION >= R_Version(2,3,0))
#define HAVE_FLUSHCONSOLE
#endif

void insnp_long_unsorted(char **filename,
			 int *nchip, char **chip_id, int *nsnps, char **snp_id,
			 int *ifX, int *female,
			 char **codes, int *ifconf, double *threshold,
			 int *progress,
			 unsigned char *gtypes, int *counts,  int *iferror) {
  int conf = *ifconf;
  int report = *progress;
  int i=0, j=0;
  int read_ok = 0; /* set to not ok for a start */

  if (report)
    message_out("Filling matrix\n");

  char *fname = *filename;
  gzFile ginfile = NULL;
  ginfile = gzopen(fname, "rb");
  if (!ginfile)
    goto open_error;

  memset(counts, 0x00, sizeof(int) * 9);

  int *called     = counts;
  int *not_called = counts +1;
  int *not_found  = counts +2;
  int *nread      = counts +3;
  int *repeats    = counts +4;
  int *disagree   = counts +5;
  int *nskip      = counts +6;
  int *no_call    = counts +7;
  int *not_conf   = counts +8;

  char *code_aa = codes[0];
  char *code_ab = codes[1];
  char *code_bb = codes[2];
  char *code_ay = NULL, *code_by = NULL;
  if (*ifX) {
    code_ay = codes[3];
    code_by = codes[4];
  }
  char chip_in[MAX_ID], snp_in[MAX_ID], gt_in[MAX_GT];
  char line_buffer[LINE_BUFFER_SIZE];
  double thr_in;
  int full = (*nsnps)*(*nchip);
  int pcpnt = full/100; /* One percentage point */
  *snp_in = '\0';
  *chip_in = '\0';
  *gt_in = '\0';

  index_db chip_db = NULL, snp_db = NULL;
  chip_db = index_create(*nchip);
  snp_db = index_create(*nsnps);

  for (i=0; i<(*nchip); i++) {
    index_insert(chip_db, chip_id[i], i);
  }
  for (j=0; j<(*nsnps); j++) {
    index_insert(snp_db, snp_id[j], j);
  }

  memset(gtypes, INVALID_GTYPE, full);

  while(1) {
    int i_snp = -1, i_chip = -1;
    char gt = 0x00;
    int get1 = 0;

    /* gzeof() is sometimes unreliable until the next read */
    if ((get1 = gzgetc(ginfile)) != -1) {
      /* get1 was successful, put it back */
      if (gzungetc(get1, ginfile) != get1) {
	*iferror = *nread;
	goto read_error;
        break;
      }
    } else {
      if (gzeof(ginfile)) {
        gzclose(ginfile);
        break;
      }
    }

    /* test for EOF then read, to make sure we get the last item */
    if (gzeof(ginfile)) {
      gzclose(ginfile);
      break;
    }
    if (gzgets(ginfile, line_buffer, LINE_BUFFER_SIZE) == Z_NULL) {
      *iferror = *nread;
      goto read_error;
      break;
    }

    if (conf==0) {
      if (sscanf(line_buffer, " %s %s %s",
		 chip_in, snp_in, gt_in)!=3) {
	*iferror = *nread;
	goto read_error;
      }
    }
    else {
      if (sscanf(line_buffer, " %s %s %s %lf",
		 chip_in, snp_in, gt_in, &thr_in)!=4) {
	*iferror = *nread;
	goto read_error;
      }
    }
    (*nread)++;

    if((i_snp = index_lookup(snp_db, snp_in)) < 0) {
      (*nskip)++;
      goto next;
    }

    if((i_chip = index_lookup(chip_db, chip_in)) < 0) {
      (*nskip)++;
      goto next;
    }

    if (*ifX && !female[i_chip]) {
      if (!strcmp(code_ay, gt_in))
	gt = 1;
      else if (!strcmp(code_by, gt_in))
	gt = 3;
      else
	gt = 0;
    }
    else {
      if (!strcmp(code_aa, gt_in))
	gt = 1;
      else if (!strcmp(code_ab, gt_in))
	gt = 2;
      else if (!strcmp(code_bb, gt_in))
	gt = 3;
      else
	gt = 0;
    }

    if (gt == 0)
      (*no_call)++;
    else if (((conf>0) && (thr_in < *threshold)) ||
	     ((conf<0) && (thr_in > *threshold))) {
      (*not_conf)++;
      gt = 0;
    }

    unsigned char *gtype = gtypes + i_snp * (*nchip) + i_chip;

    if (*gtype != INVALID_GTYPE) {
      (*repeats) ++;
      if((*gtype) && (gt) && (*gtype != gt)) {
	(*disagree) ++;
	/* the follow code is to deal with the case where
	   we have more than 2 conflicts; the 3rd time
	   is always a conflict whatever it is, but there is
	   no further need to update the called nor the gtype */
	if (*gtype != CONFLICT_GTYPE) {
	  (*called) --; /* it was counted as a call in the first instance */
	  *gtype = CONFLICT_GTYPE;
	}
      } else {
	/* don't care the calls are the same, or when the 2nd
	   is a no-call, so only deal with the case
	where the 2nd is a call and the first was no call */
	if (!(*gtype) && (gt)) {
	  (*called)++;
	  *gtype = gt;
	  (*not_called)--; /* it was no call the first time */
	}
      }
    } else {
      if (gt > 0)
	(*called) ++;
      else if (gt == 0)
	(*not_called) ++;
      *gtype = gt;
    }

  next:
    if (report && pcpnt && !((*called)%pcpnt)) {
      message_out("%d%%\r", (*called)/pcpnt);
#ifdef HAVE_FLUSHCONSOLE
      FLUSH();
#endif
    }
  }

  /* Finish normally  */

  int idx = 0;
  unsigned char *cur_gtype = gtypes;
  for (idx = 0 ; idx < full ; cur_gtype++, idx++) {
    if (*cur_gtype == INVALID_GTYPE) {
      *cur_gtype = 0x00;
      (*not_found)++;
    } else if (*cur_gtype == CONFLICT_GTYPE) {
      *cur_gtype = 0x00;
    }
  }

  if (report) {
    message_out("\nEOF reached\n");
#ifdef HAVE_FLUSHCONSOLE
    FLUSH();
#endif
  }
  *iferror = 0;
  read_ok = 1; /* everything is fine */
  /* we want drop through to tidy up the db's */

  /* Error conditions */
 read_error:
  if (!read_ok) {
    /* This is inserted because it is possible to have a read error on the first line */
    message_out("Read error after reading %d records.\n", *iferror);
  }
  index_destroy(chip_db);
  index_destroy(snp_db); 
  return;
 open_error:
  *iferror = -2;
  return;
}
