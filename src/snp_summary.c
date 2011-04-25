#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdio.h>
#include <string.h>


SEXP X_snp_summary(const SEXP Snps, const SEXP Uncertain) {

  /* SNPs ---- an X.snp.matrix */

  int *ifFemale;
  SEXP Female = R_do_slot(Snps, mkString("Female"));
  if (TYPEOF(Female)!=LGLSXP)
    error("Argument error -  Female slot wrong type");
  ifFemale = LOGICAL(Female);
  if (TYPEOF(Snps)!=RAWSXP)
    error("Argument error - Snps wrong type");
  if (Snps == R_NilValue) {
    error("Argument error - Snps = NULL");
  }
  if(!IS_S4_OBJECT(Snps)) {
    error("Argument error - Snps is not S4 object");
  }
  const unsigned char *snps = RAW(Snps);
  int N, M;
  int *dim = INTEGER(getAttrib(Snps, R_DimSymbol));
  N = dim[0];
  M = dim[1];
  SEXP names = getAttrib(Snps, R_DimNamesSymbol);
  if (names == R_NilValue) {
    error("Argument error - Snps object has no names");
  }
  SEXP snpNames = VECTOR_ELT(names, 1);
  if (snpNames == R_NilValue) {
    error("Argument error - Snps object has no snp names");
  }


  /* Handling of uncertain genotypes */

  if (TYPEOF(Uncertain) != LGLSXP)
    error("Argument error: Uncertain is wrong type");
  /* int uncert = *LOGICAL(Uncertain); */

  /* Output object */

  SEXP Result, Calls, Call_rate, MAF, P_AA, P_AB, P_BB, P_AY, P_BY, Z_HWE, Calls_female;
  PROTECT(Result = allocVector(VECSXP, 10));
  PROTECT(Calls = allocVector(INTSXP, M));
  SET_VECTOR_ELT(Result, 0, Calls);
  PROTECT(Call_rate = allocVector(REALSXP, M));
  SET_VECTOR_ELT(Result, 1, Call_rate);
  PROTECT(MAF = allocVector(REALSXP, M));
  SET_VECTOR_ELT(Result, 2, MAF);
  PROTECT(P_AA = allocVector(REALSXP, M));
  SET_VECTOR_ELT(Result, 3, P_AA);
  PROTECT(P_AB = allocVector(REALSXP, M));
  SET_VECTOR_ELT(Result, 4, P_AB);
  PROTECT(P_BB = allocVector(REALSXP, M));
  SET_VECTOR_ELT(Result, 5, P_BB);
  PROTECT(P_AY = allocVector(REALSXP, M));
  SET_VECTOR_ELT(Result, 6, P_AY);
  PROTECT(P_BY = allocVector(REALSXP, M));
  SET_VECTOR_ELT(Result, 7, P_BY);
  PROTECT(Z_HWE = allocVector(REALSXP, M));
  SET_VECTOR_ELT(Result, 8, Z_HWE);
  PROTECT(Calls_female = allocVector(INTSXP, M));
  SET_VECTOR_ELT(Result, 9, Calls_female);

  SEXP Names;
  PROTECT(Names = allocVector(STRSXP, 10));
  SET_STRING_ELT(Names, 0, mkChar("Calls"));
  SET_STRING_ELT(Names, 1, mkChar("Call.rate"));
  SET_STRING_ELT(Names, 2, mkChar("MAF"));
  SET_STRING_ELT(Names, 3, mkChar("P.AA"));
  SET_STRING_ELT(Names, 4, mkChar("P.AB"));
  SET_STRING_ELT(Names, 5, mkChar("P.BB"));
  SET_STRING_ELT(Names, 6, mkChar("P.AY"));
  SET_STRING_ELT(Names, 7, mkChar("P.BY"));
  SET_STRING_ELT(Names, 8, mkChar("z.HWE"));
  SET_STRING_ELT(Names, 9, mkChar("Calls.female"));

  int *calls = INTEGER(Calls);
  double *call_rate = REAL(Call_rate);
  double *maf = REAL(MAF);
  double *p_aa = REAL(P_AA);
  double *p_ab = REAL(P_AB);
  double *p_bb = REAL(P_BB);
  double *p_ay = REAL(P_AY);
  double *p_by = REAL(P_BY);
  double *z_hwe = REAL(Z_HWE);
  int *calls_female = INTEGER(Calls_female);

  setAttrib(Result, R_NamesSymbol, Names);
  setAttrib(Result, R_RowNamesSymbol, duplicate(snpNames));
  SEXP dfClass;
  PROTECT(dfClass = allocVector(STRSXP, 1));
  SET_STRING_ELT(dfClass, 0, mkChar("data.frame"));
  setAttrib(Result, R_ClassSymbol, dfClass);

  /* Calculations */

  int *obs = (int *) R_alloc(N, sizeof(int));
  int i, j, ij;
  for (i=0; i<N; i++)
    obs[i] = 0;

  for (j=0, ij=0; j<M; j++) {
    int aa = 0, ab = 0, bb = 0, ay=0, by=0;
    for (i=0; i<N; i++) {
      int g = (int) snps[ij++];
      if (g) {
	obs[i] = 1;
	if (ifFemale[i])
	  switch (g) {
	  case 1: aa++; break;
	  case 2: ab++; break;
	  case 3: bb++;
	  }
	else
	  switch (g) {
	  case 1: ay++; break;
	  case 3: by++;
	  }
      }
    }
    double nv = aa + ab + bb;
    double na = 2*aa + ab;
    double p = na/(2.0*nv);
    double q = 1.0 - p;
    double den = 2*p*q*sqrt(nv);
    double z = den>0.0? (ab - 2*p*q*nv)/den: NA_REAL;
    double ny = ay + by;
    double nc = 2*nv + ny;
    p = (na + ay)/nc;
    if (p>0.5)
      p = 1.0 - p;
    calls[j] = nv+ny;
    calls_female[j] = nv;
    call_rate[j] = (nv+ny)/(double) N;
    maf[j] = nc>0? p: NA_REAL;
    p_aa[j] = nv>0? ((double) aa)/nv: NA_REAL;
    p_ab[j] = nv>0? ((double) ab)/nv: NA_REAL;
    p_bb[j] = nv>0? ((double) bb)/nv: NA_REAL;
    p_ay[j] = ny>0? ((double) ay)/ny: NA_REAL;
    p_by[j] = ny>0? ((double) by)/ny: NA_REAL;
    z_hwe[j] = z;
  }
  int Nobs = 0;
  for (i=0; i<N; i++)
    Nobs += obs[i];
  if (Nobs < N) {
    warning("%d rows were empty - ignored when calculating call rates",
	    N - Nobs);
    double infl = (double) N / (double) Nobs;
    if (Nobs) {
      for (j=0; j<M; j++)
	call_rate[j] *= infl;
    }
    else {
      error("Empty matrix");
    }
  }
  UNPROTECT(13);
  return Result;
}


SEXP snp_summary(const SEXP Snps, const SEXP Uncertain) {

  /* SNPs ---- a snp.matrix */

  if (TYPEOF(Snps)!=RAWSXP)
    error("Argument error - Snps wrong type");
  if (Snps == R_NilValue) {
    error("Argument error - Snps = NULL");
  }
  if(!IS_S4_OBJECT(Snps)) {
    error("Argument error - Snps is not S4 object");
  }
  const unsigned char *snps = RAW(Snps);
  int N, M;
  int *dim = INTEGER(getAttrib(Snps, R_DimSymbol));
  N = dim[0];
  M = dim[1];
  SEXP names = getAttrib(Snps, R_DimNamesSymbol);
  if (names == R_NilValue) {
    error("Argument error - Snps object has no names");
  }
  SEXP snpNames = VECTOR_ELT(names, 1);
  if (snpNames == R_NilValue) {
    error("Argument error - Snps object has no snp names");
  }

  /* Handling of uncertain genotypes */

  if (TYPEOF(Uncertain) != LGLSXP)
    error("Argument error: Uncertain is wrong type");
  /* int uncert = *LOGICAL(Uncertain); */

  /* Output object */

  SEXP Result, Calls, Call_rate, MAF, P_AA, P_AB, P_BB, Z_HWE;
  PROTECT(Result = allocVector(VECSXP, 7));
  PROTECT(Calls = allocVector(INTSXP, M));
  SET_VECTOR_ELT(Result, 0, Calls);
  PROTECT(Call_rate = allocVector(REALSXP, M));
  SET_VECTOR_ELT(Result, 1, Call_rate);
  PROTECT(MAF = allocVector(REALSXP, M));
  SET_VECTOR_ELT(Result, 2, MAF);
  PROTECT(P_AA = allocVector(REALSXP, M));
  SET_VECTOR_ELT(Result, 3, P_AA);
  PROTECT(P_AB = allocVector(REALSXP, M));
  SET_VECTOR_ELT(Result, 4, P_AB);
  PROTECT(P_BB = allocVector(REALSXP, M));
  SET_VECTOR_ELT(Result, 5, P_BB);
  PROTECT(Z_HWE = allocVector(REALSXP, M));
  SET_VECTOR_ELT(Result, 6, Z_HWE);

  SEXP Names;
  PROTECT(Names = allocVector(STRSXP, 7));
  SET_STRING_ELT(Names, 0, mkChar("Calls"));
  SET_STRING_ELT(Names, 1, mkChar("Call.rate"));
  SET_STRING_ELT(Names, 2, mkChar("MAF"));
  SET_STRING_ELT(Names, 3, mkChar("P.AA"));
  SET_STRING_ELT(Names, 4, mkChar("P.AB"));
  SET_STRING_ELT(Names, 5, mkChar("P.BB"));
  SET_STRING_ELT(Names, 6, mkChar("z.HWE"));

  int *calls = INTEGER(Calls);
  double *call_rate = REAL(Call_rate);
  double *maf = REAL(MAF);
  double *p_aa = REAL(P_AA);
  double *p_ab = REAL(P_AB);
  double *p_bb = REAL(P_BB);
  double *z_hwe = REAL(Z_HWE);

  setAttrib(Result, R_NamesSymbol, Names);
  setAttrib(Result, R_RowNamesSymbol, duplicate(snpNames));
  SEXP dfClass;
  PROTECT(dfClass = allocVector(STRSXP, 1));
  SET_STRING_ELT(dfClass, 0, mkChar("data.frame"));
  setAttrib(Result, R_ClassSymbol, dfClass);

  /* Calculations */

  int *obs = (int *) R_alloc(N, sizeof(int));
  int i, j, ij;
  for (i=0; i<N; i++)
    obs[i] = 0;

  for (j=0, ij=0; j<M; j++) {
    int aa = 0, ab = 0, bb = 0;
    for (i=0; i<N; i++) {
      int g = (int) snps[ij++];
      if (g) {
	obs[i] = 1;
	switch (g) {
	case 1: aa++; break;
	case 2: ab++; break;
	case 3: bb++;
	}
      }
    }
    double nv = aa + ab + bb;
    double na = 2*aa + ab;
    double p =  na/(2*nv);
    double q = 1.0 - p;
    double den = 2*p*q*sqrt(nv);
    double z = den>0.0? (ab - 2*p*q*nv)/den: NA_REAL;
    if (p>0.5)
      p = 1.0 - p;
    calls[j] = nv;
    call_rate[j] = nv/(double) N;
    maf[j] = nv>0? p: NA_REAL;
    p_aa[j] = nv>0? ((double) aa)/nv: NA_REAL;
    p_ab[j] = nv>0? ((double) ab)/nv: NA_REAL;
    p_bb[j] = nv>0? ((double) bb)/nv: NA_REAL;
    z_hwe[j] = z;
  }

  int Nobs = 0;
  for (i=0; i<N; i++)
    Nobs += obs[i];
  if (Nobs < N) {
    warning("%d rows were empty - ignored when calculating call rates",
	    N - Nobs);
    double infl = (double) N / (double) Nobs;
    if (Nobs) {
      for (j=0; j<M; j++)
	call_rate[j] *= infl;
    }
    else {
      error("Empty matrix");
    }
  }
  UNPROTECT(10);
  return Result;
}


SEXP row_summary(const SEXP Snps) {

  if (TYPEOF(Snps)!=RAWSXP)
    error("Argument error - Snps wrong type");
  if (Snps == R_NilValue) {
    error("Argument error - Snps = NULL");
  }
  if(!IS_S4_OBJECT(Snps)) {
    error("Argument error - Snps is not S4 object");
  }
  const unsigned char *snps = RAW(Snps);
  int N, M;
  int *dim = INTEGER(getAttrib(Snps, R_DimSymbol));
  N = dim[0];
  M = dim[1];
  SEXP names = getAttrib(Snps, R_DimNamesSymbol);
  if (names == R_NilValue) {
    error("Argument error - Snps object has no names");
  }
  SEXP rowNames = VECTOR_ELT(names, 0);
  if (rowNames == R_NilValue) {
    error("Argument error - Snps object has no row names");
  }
  /* Output object */

  SEXP Result, Call_rate, Het;
  PROTECT(Result = allocVector(VECSXP, 2));
  PROTECT(Call_rate = allocVector(REALSXP, N));
  SET_VECTOR_ELT(Result, 0, Call_rate);
  PROTECT(Het = allocVector(REALSXP, N));
  SET_VECTOR_ELT(Result, 1, Het);

  SEXP Names;
  PROTECT(Names = allocVector(STRSXP, 2));
  SET_STRING_ELT(Names, 0, mkChar("Call.rate"));
  SET_STRING_ELT(Names, 1, mkChar("Heterozygosity"));

  double *call_rate = REAL(Call_rate);
  double *het = REAL(Het);

  setAttrib(Result, R_NamesSymbol, Names);
  setAttrib(Result, R_RowNamesSymbol, duplicate(rowNames));
  SEXP dfClass;
  PROTECT(dfClass = allocVector(STRSXP, 1));
  SET_STRING_ELT(dfClass, 0, mkChar("data.frame"));
  setAttrib(Result, R_ClassSymbol, dfClass);

  /* Calculations */
  int i, j, ij;
  for (i=0; i<N; i++) {
    int ncall = 0, nhet=0;
    for (j=0, ij=i; j<M; j++, ij+=N) {
      unsigned char g = snps[ij];
      if (g) {
	ncall++;
	if (g == 0x02) /* Het */
	  nhet++;
      }
    }
    call_rate[i] = (double) ncall/ (double) M;
    het[i] = (double) nhet/ (double) ncall;
  }
  UNPROTECT(5);
  return Result;
}
