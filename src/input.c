/* Read genotype data in long format

filename    Name of input text file (4 fields per row) 
ifconf      If +ve, high score on confidence score is good, if -ve high 
            score is bad. If zero no score is present in input.
threshold   Quality threshold for inclusion of data
nchip       Number of chips
chip_id     Array of chip identifier strings
nsnp        Number of snps
snp_id      Array of snp identifier strings
ifX         Indicator for X chromosome
female      If X chromosome, indicator array for female sex
progress    1 if a progess report to stdout, else 0.
gtypes      Array of char's, length nchip*nsnp, to hold genotype data
counts      Various counts
iferror     Error indicator

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <R.h>
#include "zlib.h"

#define MAX_ID 64
#define MAX_GT 16
#define LINE_BUFFER_SIZE 1024

void insnp_long(char **filename, 
	   int *nchip, char **chip_id, int *nsnps, char **snp_id,
		int *ifX, int *female,
		char **codes, int *ifconf, double *threshold, 
		int *progress,
                unsigned char *gtypes, int *counts,  int *iferror) {
  int conf = *ifconf;
  int report = *progress;
  if (report) 
    Rprintf("\nInput file assumed already sorted\nFilling matrix\n");

  char *fname = *filename;
  gzFile ginfile = NULL;
  ginfile = gzopen(fname, "rb");
  if (!ginfile)
    goto open_error;

  int nread = 0;
  int nskip = 0;
  int no_call = 0;
  int not_called = 0;
  int not_conf = 0;
  int called = 0;
  int not_found = 0;
  int repeats = 0;
  int disagree = 0;

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
  char *snp_target, *chip_target;
  double thr_in;
  int cmp=-1;
  int full = (*nsnps)*(*nchip);
  int pcpnt = full/100; /* One percentage point */
  int ij=0;
  int i = 0, j = 0;
  *snp_in = (char) 0;
  *chip_in = (char) 0;
  int incmplt = 0;
  for (j=0; j<(*nsnps); j++) {
    incmplt = 1;
    snp_target = snp_id[j];
    for (i=0; i<(*nchip); i++, incmplt=0) {
      int gij = -1;
      chip_target = chip_id[i];
      int jcmp = strcoll(snp_in, snp_target);
      if (jcmp>0)
	cmp = 1;
      else if (jcmp==0) {
	int icmp = strcoll(chip_in, chip_target);
	cmp = icmp<0 ? -1 : (icmp>0 ? +1: 0);
      }
      else 
	cmp = -1;

      while (cmp<=0) { /* Read until ahead */
	int get1 = 0;
	/* If At, update current genotype, gij, with read value, gijw */

	int gijw;
	if (!cmp) { /* record found */
	  if (*ifX && !female[i]) { 
	    if (!strcmp(code_ay, gt_in)) 
	      gijw = 1;
	    else if (!strcmp(code_by, gt_in)) 
	      gijw = 3;
	    else
	      gijw = 0;
	  }
	  else {
	    if (!strcmp(code_aa, gt_in))
	      gijw = 1;
	    else if (!strcmp(code_ab, gt_in))
	      gijw = 2;
	    else if (!strcmp(code_bb, gt_in))
	      gijw = 3;
	    else
	      gijw = 0;
	  }
	  if (gijw == 0) 
	    no_call++;
	  else if (((conf>0) && (thr_in < *threshold)) ||
		   ((conf<0) && (thr_in > *threshold))) {
	    not_conf++;
	    gijw = 0;
	  }
	  /* Update gij */
	  if (gij < 0) 
	    gij = gijw;
	  else {
	    repeats ++;
	    if (gij < 4) {
	      if (gijw) {
		if (!gij)
		  gij = gijw;
		else if (gijw != gij)
		  gij = 4; /* Call conflict */
	      }
	    }
	  }
	}
	else if (*chip_in) {
	  nskip++;
	}

	/* Read next record */

	gzgets(ginfile, line_buffer, LINE_BUFFER_SIZE);
	if (gzeof(ginfile)) {
	  gzclose(ginfile); 
	  cmp = 2;
	  break;
	}
	/* gzeof() is sometimes unreliable until the next read */
	if ((get1 = gzgetc(ginfile)) != -1) {
	  /* get1 was successful, put it back */
	  if (gzungetc(get1, ginfile) != get1) {
	    goto read_error;
	  }
	} else {
	  if (gzeof(ginfile)) {
	    gzclose(ginfile); 
	    cmp = 2;
	    break;
	  }
	}
	
	if (conf==0) {
	  if (sscanf(line_buffer, " %s %s %s", 
		     chip_in, snp_in, gt_in)!=3)
	    goto read_error;
	}
	else {
	  if (sscanf(line_buffer, " %s %s %s %lf", 
		     chip_in, snp_in, gt_in, &thr_in)!=4)
	    goto read_error;
	}
	nread++;

	/* Before, At or Past target? */

	jcmp = strcoll(snp_in, snp_target);
	if (jcmp>0)
	  cmp = 1;
	else if (jcmp==0) {
	  int icmp = strcoll(chip_in, chip_target);
	  cmp = icmp<0 ? -1 : (icmp>0 ? +1: 0);
	}
	else 
	  cmp = -1;
      }

      if (gij == 0)
	not_called ++;
      else if (gij == 4) { 
	disagree ++;
	gij = 0;
      }
      else if (gij > 0)
	called ++;
      else {
	not_found ++;
	gij = 0;
      }
      gtypes[ij++] = (unsigned char) gij;
      if (report && pcpnt && !(called % pcpnt)) {
	Rprintf("%d%%\r", called/pcpnt);
	R_FlushConsole();
      }
      if (cmp==2) 
	break;
    }
    if (cmp==2)
      break;
  }
  
  /* Finish normally  */

  if (report) {
    if (cmp==2) 
      Rprintf("\nEOF reached\n");
    else
      Rprintf("\n");
    R_FlushConsole();
  }
  while (ij<full) {
    not_found ++;
    gtypes[ij++] = (unsigned char)0;
  }
  counts[0] = called;
  counts[1] = not_called;
  counts[2] = not_found;
  counts[3] = nread;
  counts[4] = repeats;
  counts[5] = disagree;
  counts[6] = nskip;
  counts[7] = no_call;
  counts[8] = not_conf;
  if (cmp==2 && incmplt) 
    *iferror = -3; /* EOF before matrix filled */
  else
    *iferror = 0;
  return;
  
  /* Error conditions */ 
  
 open_error: 
  *iferror = (-2); return;
 read_error: 
  Rprintf("Read error after reading %d records.\n", nread);
  *iferror = nread; return;
}

/* Find non-empty rows and columns */

void empty(int *N, int *M, unsigned char *matrix, int *row, int *col){
  int i, j, ij;
  for (i=0; i<(*N); i++)
    row[i] = 0;
  for (j=0, ij=0; j<(*M); j++) {
    int colj = 0;
    for (i=0; i<(*N); i++, ij++) {
      if(matrix[ij]) {
	row[i] = 1;
	colj = 1;
      }
    }
    col[j] = colj;
  }
}
