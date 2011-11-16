/* 
   Single SNP tests - 1df and 2df, controlling for a stratification 
   
   Calculates score test and permutation variance (variance under 
   random permutations of phenotype vector within strata) and then the 
   conventional U^2/V chi-squared test
 
*/

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Rmissing.h"

SEXP one_at_a_time(const SEXP Phenotype, const SEXP Stratum, const SEXP Snps, 
		   const SEXP Subset, const SEXP Snp_subset) {
  const double tol=1.e-2;
  double vec[4] = {1.0, 0.0, 0.0, 0.0};
  int i, j, k, km, m, t, su;

  /* Phenotype */

  if (TYPEOF(Phenotype)!=REALSXP)
    error("Argument error - Phenotype");
  const double *phenotype = REAL(Phenotype);
  int n = LENGTH(Phenotype);
  
  /* Stratum */

  const int *stratum = NULL;
  int nstrata = 1;
  SEXPTYPE stype = TYPEOF(Stratum);
  if (stype==INTSXP) {
    if (LENGTH(Stratum)!=n)
      error("Dimension error - Stratum");
    stratum = INTEGER(Stratum);
    for (i=0; i<n; i++) {
      int si = stratum[i];
      if (si>nstrata) nstrata = si;
    }
  }
  else if (stype!=NILSXP)
    error("Argument error - Stratum");

  /* SNPs ---- should be a snp.matrix or an X.snp.matrix */

  const char *classS = NULL;
  if (TYPEOF(R_data_class(Snps, FALSE)) == STRSXP) {
    classS = CHAR(STRING_ELT(R_data_class(Snps, FALSE), 0));
  } else {
    classS = CHAR(STRING_ELT(getAttrib(Snps, R_ClassSymbol), 0));
  }
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
  if (N!=n)
    error("Dimension error - Snps");

  /* Subset */

  int *subset = NULL;
  int nsubj = n;
  SEXPTYPE sutype = TYPEOF(Subset);
  if (sutype==INTSXP) {
    if (LENGTH(Subset)>n)
      error("Dimenion error - Subset");
    subset = INTEGER(Subset);
    nsubj = LENGTH(Subset);
  }
  else if (sutype!=NILSXP)
    error("Argument error - Subset");

  
  /* SNP subset */

  int *snp_subset = NULL;
  int ntest = nsnp;
  SEXPTYPE sstype = TYPEOF(Snp_subset);
  if (sstype==INTSXP) {
    if (LENGTH(Snp_subset)>nsnp)
      error("Dimenion error - Snp_subset");
    snp_subset = INTEGER(Snp_subset);
    ntest = LENGTH(Snp_subset);
  }
  else if (sstype!=NILSXP)
    error("Argument error - Snp.subset");

  /* Sex */

  const int *female = NULL;
  if (ifX) {
    SEXP Female = R_do_slot(Snps, mkString("Female")); 
    female = LOGICAL(Female);
  }
   

  /* Output objects */
 
  SEXP Result, Chi_1df, Chi_2df, p_1df, p_2df, Used;
  PROTECT(Result = allocVector(VECSXP, 5));
  PROTECT(Chi_1df =  allocVector(REALSXP, ntest) ); /* Chi-squared(1) */
  SET_VECTOR_ELT(Result, 0, Chi_1df);
  double *chi_1df = REAL(Chi_1df); 
  PROTECT(Chi_2df =  allocVector(REALSXP, ntest) ); /* Chi-squared(2) */
  SET_VECTOR_ELT(Result, 1, Chi_2df);
  double *chi_2df = REAL(Chi_2df);
  PROTECT(p_1df =  allocVector(REALSXP, ntest) ); /* p.1df */
  PROTECT(p_2df =  allocVector(REALSXP, ntest) ); /* p.2df */
  SET_VECTOR_ELT(Result, 2, p_1df);
  SET_VECTOR_ELT(Result, 3, p_2df);
  PROTECT(Used =  allocVector(INTSXP, ntest) ); /* N used */
  SET_VECTOR_ELT(Result, 4, Used);
  int *Nused = INTEGER(Used);


  /* Do calculations */

  
  if (!ifX) {

    /* Work arrays */
    
    double **UV = (double **) calloc(nstrata, sizeof(double *));
    for (i=0; i<nstrata; i++) 
      UV[i] = (double *) calloc(10, sizeof(double));
    
    for (t=0; t<ntest; t++) {
      int i = snp_subset? snp_subset[t] - 1: t; 
      const unsigned char *snpsi = snps + n*i;
      /* Initialise score and score variance array */
      for (j=0; j<nstrata; j++) {
	double *uvj = UV[j];
	for (k=0; k<10; k++)
	  uvj[k] = 0.0;
      }
      int nu = 0;
      for (su=0; su<nsubj; su++) {
	int j = subset? subset[su] - 1: su;
	unsigned char zij = snpsi[j];
	if (zij) {
	  nu++;
	  vec[1] = phenotype[j];
	  vec[2] = zij - 1;
	  vec[3] = (zij==2);
	  int strat = (nstrata==1)? 0: (stratum[j]-1);
	  double *uv = UV[strat];
	  for (k=0, km=0; k<4; k++) {
	    double vk = vec[k];
	    for (m=0; m<=k; m++)
	      uv[km++] += vk*vec[m];
	  }
	}
      }
      double u1=0.0, u2=0.0, v11=0.0, v12=0.0, v22=0.0; 
      for (j=0; j<nstrata; j++) {
	double *uvj = UV[j];
	double N = uvj[0];
	if (N>1.0) {
	  u1 += (uvj[4] - uvj[3]*uvj[1]/N);
	  u2 += (uvj[7] - uvj[6]*uvj[1]/N);
	  double vy = (uvj[2] - uvj[1]*uvj[1]/N)/(N-1.0);
	  v11 += vy*(uvj[5] - uvj[3]*uvj[3]/N);
	  v12 += vy*(uvj[8] - uvj[3]*uvj[6]/N);
	  v22 += vy*(uvj[9] - uvj[6]*uvj[6]/N);
	}
      }
      /* 1 df test */
      double u = u1;
      double v = v11;
      chi_1df[t] =  v>0.0? u*u/v: NA_REAL;
      /* 2 df test */
      double r2 = v12*v12/(v11*v22);
      if (v11 <= 0.0 || v22 <= 0.0 || (1.0-r2) < tol) 
	chi_2df[t] = NA_REAL; /* Not positive definite --- enough! */
      else 
	chi_2df[t] = (u1*u1/v11 + u2*u2/v22 - 2.0*r2*u1*u2/v12)/(1.0-r2);
      Nused[t] = nu;
    }

    /* Return work arrays */

    for (i=0; i<nstrata; i++) 
      free(UV[i]);
    free(UV);
  }
  else {

    /* Work arrays */
    
    double **UVM = (double **) calloc(nstrata, sizeof(double *));
    double **UVF = (double **) calloc(nstrata, sizeof(double *));
    for (i=0; i<nstrata; i++) {
      UVM[i] = (double *) calloc(10, sizeof(double));
      UVF[i] = (double *) calloc(10, sizeof(double));
    }
    
    for (t=0; t<ntest; t++) {
      int i = snp_subset? snp_subset[t] - 1: t; 
      const unsigned char *snpsi = snps + n*i;
      /* Initialise score and score variance arrays */
      for (j=0; j<nstrata; j++) {
	double *uvmj = UVM[j];
	double *uvfj = UVF[j];
	for (k=0; k<10; k++)
	  uvmj[k] = uvfj[k] = 0.0;
      }
      int nu = 0;
      for (su=0; su<nsubj; su++) {
	int j = subset? subset[su] - 1: su;
	unsigned char zij = snpsi[j];
	if (zij) {
	  nu++;
	  vec[1] = phenotype[j];
	  vec[2] = zij - 1;
	  vec[3] = (zij==2);
	  int strat = (nstrata==1)? 0: (stratum[j]-1);
	  double *uv = female[j] ? UVF[strat]: UVM[strat];
	  for (k=0, km=0; k<4; k++) {
	    double vk = vec[k];
	    for (m=0; m<=k; m++)
	      uv[km++] += vk*vec[m];
	  }
	}
      }
      /* rewrite this loop */
      double u=0, v=0, u1=0.0, u2=0.0, v11=0.0, v12=0.0, v22=0.0;
      for (j=0; j<nstrata; j++) {
	double *uvmj = UVM[j];
	double Nm = uvmj[0];
	double *uvfj = UVF[j];
	double Nf = uvfj[0];
	double Nt = Nf + Nm;
	if (Nt>0) {
	  double ybar = (uvmj[1] + uvfj[1])/Nt;
	  double yb2 = ybar*ybar;
	  double af = (uvfj[3] + uvmj[3]/2)/(2*Nf+Nm);
	  u += uvmj[4] + uvfj[4] - ybar*(uvmj[3]+uvfj[3]);
	  /* Variance of male contribution (Bernoulli) */
	  if (Nm>0) {
	    double ssy = uvmj[2] - 2*ybar*uvmj[1] + Nm*yb2;
	    v += 4*af*(1.0-af)*ssy;
	  }
	  /* Female contribution */
	  if (Nf>1) {
	    u1 += uvfj[4] - uvfj[3]*uvfj[1]/Nf;
	    u2 += uvfj[7] - uvfj[6]*uvfj[1]/Nf;
	    double ssy = (uvfj[2] - 2*ybar*uvfj[1] + Nf*yb2)/(Nf-1.0);
	    double w = (uvfj[5] - uvfj[3]*uvfj[3]/Nf);
	    v += ssy*w;
	    ssy = (uvfj[2] - uvfj[1]*uvfj[1]/Nf)/(Nf-1.0);
	    v11 += ssy*w;
	    v12 += ssy*(uvfj[8] - uvfj[3]*uvfj[6]/Nf);
	    v22 += ssy*(uvfj[9] - uvfj[6]*uvfj[6]/Nf);
	  }
	}
      }
      /* 1 df test */
      chi_1df[t] =  v>0.0? u*u/v: NA_REAL;
      /* 2 df test */
      double r2 = v12*v12/(v11*v22);
      if (v11 <= 0.0 || v22 <= 0.0 || (1.0-r2) < tol) 
	chi_2df[t] = NA_REAL; /* Not positive definite --- enough! */
      else 
	chi_2df[t] = chi_1df[t] + 
	  (r2*u1*u1/v11 + u2*u2/v22 - 2.0*r2*u1*u2/v12)/(1.0-r2);
      Nused[t] = nu;
    }

    /* Return work arrays */

    for (i=0; i<nstrata; i++) {
      free(UVM[i]);
      free(UVF[i]);
    }
    free(UVM);
    free(UVF);
  }

  /* calculating the p-values */
  for (t=0; t<ntest; t++) {
    REAL(p_1df)[t] = pchisq(chi_1df[t], 1.0, 0,0);
    REAL(p_2df)[t] = pchisq(chi_2df[t], 2.0, 0,0);
  }

  /* Attributes of output object */

  SEXP cNames, rnames, dfClass;
  SEXP snpNames = VECTOR_ELT(getAttrib(Snps, R_DimNamesSymbol), 1);
  SEXPTYPE sntype = TYPEOF(snpNames);
  PROTECT(rnames = allocVector(sntype, ntest));
  for (t=0; t<ntest; t++) {
    int i = snp_subset? snp_subset[t] - 1: t; 
    if (sntype==INTSXP) 
      INTEGER(rnames)[t] = INTEGER(snpNames)[i];
    else 
      SET_STRING_ELT(rnames, t, mkChar(CHAR(STRING_ELT(snpNames,i))));
  } 
  setAttrib(Result, R_RowNamesSymbol,  rnames);

  PROTECT(cNames = allocVector(STRSXP, 5));
  SET_STRING_ELT(cNames, 0, mkChar("chi2.1df"));
  SET_STRING_ELT(cNames, 1, mkChar("chi2.2df"));
  SET_STRING_ELT(cNames, 2, mkChar("p.1df"));
  SET_STRING_ELT(cNames, 3, mkChar("p.2df"));
  SET_STRING_ELT(cNames, 4, mkChar("N"));
  setAttrib(Result, R_NamesSymbol, cNames);
  PROTECT(dfClass = allocVector(STRSXP, 1));
  SET_STRING_ELT(dfClass, 0, mkChar("data.frame"));
  setAttrib(Result, R_ClassSymbol, dfClass);
  UNPROTECT(9);

  return(Result);
}
