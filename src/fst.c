#include <string.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/* Fst calculation -- as used in HapMap */

SEXP Fst(SEXP Snps, SEXP Group) {
  
  /* Process Snps argument */

  const char *classS = NULL;

  classS = CHAR(STRING_ELT(getAttrib(Snps, R_ClassSymbol), 0));
  int ifX = 0; /* default Not X */
  if (!strncmp(classS, "snp", 3))
    ifX = 0;
  else if (!strncmp(classS, "X.snp", 5))
    ifX = 1;
  else {
    ifX = 0; /* to avoid warning message */
    error("Argument error - class(Snps)");
  }
  if(!IS_S4_OBJECT(Snps)) {
    error("Argument error - Snps is not S4 object");
  }

  if (TYPEOF(Snps)!=RAWSXP)
    error("Argument error - Snps");
  const unsigned char *snps = RAW(Snps);
  int N, nsnp;
  if (strlen(classS)>5) {
    int *dim = INTEGER(getAttrib(Snps, R_DimSymbol));
    N = dim[0];
    nsnp = dim[1];
  }
  else {
    N = LENGTH(Snps);
    nsnp = 1;
  }
     
  const int *female = NULL;
  if (ifX) {
    SEXP Female = R_do_slot(Snps, mkString("Female")); 
    female = LOGICAL(Female);
  }

  /* Process Group argument */

  classS = CHAR(STRING_ELT(getAttrib(Group, R_ClassSymbol), 0));
  if (strcmp(classS, "factor"))
    error("Argument error - class(Group)");
  if (LENGTH(Group)!=N)
    error("Non-conformant arguments");
  int ngrp = nlevels(Group);
  int *group = INTEGER(Group);

  /* Setup result object */

  SEXP Fst, Weight;
  PROTECT(Fst = allocVector(REALSXP, nsnp));
  PROTECT(Weight = allocVector(REALSXP, nsnp));
  double *fst = REAL(Fst);
  double *weight = REAL(Weight);

  /* Work arrays */

  int *na2 = (int *) Calloc(ngrp, int);
  int *na = (int *) Calloc(ngrp, int);
  double *gwts = (double *) Calloc(ngrp, double);
  
  /* Calculate group weights */

  memset(na, 0x00, ngrp*sizeof(int));
  for (int i=0; i<N; i++) {
    int gi = group[i];
    if (gi != NA_INTEGER) {
      gi--;
      if (ifX) {
	int fi = female[i];
	  na[gi] += fi? 2: 1;
      }
      else {
	na[gi] += 2;
      }	
    }
  }
  double sgw = 0.0;
  for (int g=0; g<ngrp; g++) {
    double w = (double)na[g]*(double)(na[g]-1);
    sgw += w;
    gwts[g] = w;
  }
  for (int g=0; g<ngrp; g++) 
    gwts[g] /= sgw;


  /* Calculate Fst for each SNP, plus its weight in the overall estimate  */
    
  for (int j=0, ij=0; j<nsnp; j++) {
    memset(na, 0x00, ngrp*sizeof(int));
    memset(na2, 0x00, ngrp*sizeof(int));
    for (int i=0; i<N; i++, ij++) {
      int gi = group[i];
      int sij = (int) snps[ij];
      if (gi != NA_INTEGER && sij) {
	gi--;
	if (ifX) {
	  int fi = female[i];
	  na[gi] += fi? 2: 1;
	  na2[gi] += fi? (sij-1): (sij==3);
	}
	else {
	  na[gi] += 2;
	  na2[gi] += (sij-1);
	}
      }
    }
    int nt = 0, nt2 = 0;
    double pq = 0.0;
    for (int g=0; g<ngrp; g++) {
      int nag = na[g];
      if (nag>1) {
	int na2g = na2[g];
	double pg = (double)na2g/ (double)nag;
	nt += nag;
	nt2 += na2g;
	pq += gwts[g]*pg*(1.0-pg)*(double)nag/(double)(nag-1);
      }
    }
    if (nt>1) {
      double p = (double)nt2 / (double)nt;
      double pqm = p*(1.0-p)*(double)nt/(double)(nt-1);
      fst[j] = 1.0-pq/pqm;
      weight[j] = pqm;
    }
    else {
      fst[j] = NA_REAL;
      weight[j] = NA_REAL;
    }
  }

  /* Tidy up and return */

  Free(na);
  Free(na2);
  Free(gwts);

  SEXP Result, Names;
  PROTECT(Result = allocVector(VECSXP, 2));
  PROTECT(Names = allocVector(STRSXP, 2));
  SET_STRING_ELT(Names, 0, mkChar("Fst"));
  SET_STRING_ELT(Names, 1, mkChar("weight"));
  setAttrib(Result, R_NamesSymbol, Names);
  SET_VECTOR_ELT(Result, 0, Fst);
  SET_VECTOR_ELT(Result, 1, Weight);
  
  UNPROTECT(4);
  return Result;
}

      