/* 
   Read a plink .bed file as a snp.matrix
 
*/


#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Rmissing.h"


SEXP readbed(SEXP Bed, SEXP Id, SEXP Snps) {
  const unsigned char recode[4] = {'\x01', '\x00', '\x02', '\x03'};
  const unsigned char mask = '\x03';
  int nrow = LENGTH(Id);
  int ncol = LENGTH(Snps);
  const char *file = CHAR(STRING_ELT(Bed, 0));
  FILE *in = fopen(file, "r");
  if (!in)
    error("Couln't open input file: %s", file);
  unsigned char start[3];
  if (fread(start, 1, 3, in)!=3)
    error("Failed to read first 3 bytes");
  if (start[0]!='\x6C' || start[1]!='\x1B')
    error("Input file does not appear to be a .bed file (%X, %X)", 
	  start[0], start[1]);

  /* Create output object */

  SEXP Result, Dimnames, Package, Class;
  PROTECT(Result = allocMatrix(RAWSXP, nrow, ncol));
  PROTECT(Dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(Dimnames, 0, Id);
  SET_VECTOR_ELT(Dimnames, 1, Snps);
  setAttrib(Result, R_DimNamesSymbol, Dimnames);
  PROTECT(Class = allocVector(STRSXP, 1));
  SET_STRING_ELT(Class, 0, mkChar("snp.matrix"));
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("snpMatrix"));
  setAttrib(Class, install("package"), Package);
  classgets(Result, Class);
  SET_S4_OBJECT(Result);
  
  unsigned char *result = RAW(Result); 
  int ncell = nrow*ncol;
  memset(result, 0x00, ncell);

  /* Read in data */

  int snp_major = start[2];
  int part=0, ij=0, i=0, j=0;
  while (1) {
    unsigned char byte;
    if (!part) {
      if (feof(in) || !fread(&byte, 1, 1, in)) 
	error("Unexpected end of file reached");
      part = 4;
    }
    unsigned char code = byte & mask;
    byte = byte >> 2;
    part--;
    result[ij] = recode[code];
    if (snp_major) {
      ij++;
      i++;
      if (i==nrow) {
	i = part = 0;
	j++;
	if (j==ncol)
	  break;
      }
    }	
    else {
      ij += nrow;
      j++;
      if (j==ncol){
	j = part = 0;
	i++;
	if (i==nrow)
	  break;
	ij = i;
      }
    }
  }

  UNPROTECT(4);
  return Result;
}

