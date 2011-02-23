/* 
   Single trio-based SNP tests - 1df and 2df
   
   Calculates score test and variance 
 
*/

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hash_index.h"
#include "imputation.h"
#include "Rmissing.h"

SEXP score_tdt(const SEXP Proband, const SEXP Father, const SEXP Mother, 
	       const SEXP Cluster, const SEXP Snps,  const SEXP Rules, 
               const SEXP Snp_subset, const SEXP Check, const SEXP Robust){
 
  int mendelian[36] = {1, 0, 0, 1, 1, 0, 0, 1, 0,
		       1, 1, 0, 1, 1, 1, 0, 1, 1,
		       0, 1, 0, 0, 1, 1, 0, 0, 1,
                       1, 0, 0, 1, 0, 1, 0, 0, 1};
   
  /* Trios */

  if (TYPEOF(Proband)!=INTSXP || TYPEOF(Father)!=INTSXP || 
      TYPEOF(Mother)!=INTSXP || TYPEOF(Cluster)!=INTSXP )
    error("Argument error - Proband|Father|Mother|Cluster");
  const int *proband = INTEGER(Proband);
  const int *father = INTEGER(Father);
  const int *mother = INTEGER(Mother);
  const int *cluster = INTEGER(Cluster);

  int ntrio = LENGTH(Proband);
  
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
  int nsubj, nsnp;
  if (strlen(classS)>5) {
    int *dim = INTEGER(getAttrib(Snps, R_DimSymbol));
    nsubj = dim[0];
    nsnp = dim[1];
  }
  else {
    nsubj = LENGTH(Snps);
    nsnp = 1;
  }
  SEXP snp_names = VECTOR_ELT(getAttrib(Snps, R_DimNamesSymbol), 1);
  index_db name_index = create_name_index(snp_names);

  /* Rules */

  int nrules = 0;
  if (!isNull(Rules)) {
    const char *classR = NULL;
    if (TYPEOF(R_data_class(Rules, FALSE)) == STRSXP) {
      classR = CHAR(STRING_ELT(R_data_class(Rules, FALSE), 0));
    } else {
      classR = CHAR(STRING_ELT(getAttrib(Rules, R_ClassSymbol), 0));
    }
    if (strcmp(classR, "snp.reg.imputation")!=0) 
      error("Argument error - Rules");
    nrules = LENGTH(Rules);
  }

  /* SNP subset */

  int *snp_subset = NULL;
  int ntest = nrules? nrules: nsnp;
  SEXPTYPE sstype = TYPEOF(Snp_subset);
  if (sstype==INTSXP) {
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

  /* Check inheritance */

  if (TYPEOF(Check)!=LGLSXP)
    error("Argument error `Check'");
  int check = *LOGICAL(Check);
 
  /* Robust option */

  if (TYPEOF(Robust)!=LGLSXP)    
    error("Argument error `Robust'");
  int robust = *LOGICAL(Robust);

  /* Output objects */
 
  SEXP Result, Used, U, V;
  PROTECT(Result = allocVector(VECSXP, 4));
  
  if (female) {
    PROTECT(U =  allocMatrix(REALSXP, ntest, 3) ); /* Scores */
    PROTECT(V =  allocMatrix(REALSXP, ntest, 4) ); /* Score variances */
  }
  else {
    PROTECT(U =  allocMatrix(REALSXP, ntest, 2) ); /* Scores */
    PROTECT(V =  allocMatrix(REALSXP, ntest, 3) ); /* Score variances */
  }
  SET_VECTOR_ELT(Result, 0, U);
  double *umat = REAL(U); 
  SET_VECTOR_ELT(Result, 1, V);
  double *vmat = REAL(V);

  PROTECT(Used =  allocVector(INTSXP, ntest) ); /* N used */
  int *Nused = INTEGER(Used);
  SET_VECTOR_ELT(Result, 2, Used);

  SEXP Nr2;
  double *nr2;
  if (isNull(Rules)) {
    PROTECT(Nr2 = allocVector(REALSXP, 0));
    nr2 = NULL;
  }
  else {
    PROTECT(Nr2 = allocVector(REALSXP, ntest)); /* N used x R-squared */
    nr2 = REAL(Nr2);
  }
  SET_VECTOR_ELT(Result, 3, Nr2);
  
  /* Space to hold imputed values */

  double *ximp = NULL;
  if (nrules) {
    ximp = (double *) Calloc(nsubj, double);
  }

  /* Do calculations */
  
  int nonmend=0, Xerrors=0;
  
  const double thrsix = 3.0/16.0;
  const double quart = 1.0/4.0;

  if (!ifX) {

    for (int t=0; t<ntest; t++) {
      int i = snp_subset? snp_subset[t] - 1: t; 
      const unsigned char *snpsi = NULL;
      double r2 = 0.0;
      if (nrules) {
	if (i >= nrules)
	  error("snp_subset out of range");
	SEXP Rule =  VECTOR_ELT(Rules, i);
	if (isNull(Rule)){ /* Monomorphic */
	  for (int j=0; j<nsubj; j++){
	    ximp[j] = 0.0;
	  }
	}
	else {
	  do_impute(snps, nsubj, NULL, nsubj, name_index, Rule, ximp, NULL);
	  r2 = *REAL(VECTOR_ELT(Rule, 0));
	}
      }
      else {
	if (i >= nsnp)
	  error("snp_subset out of range");
	snpsi = snps + nsubj*i;
      }

      /* Initialise score and score variance array */

      double u1=0.0, u2=0.0, v11=0.0, v12=0.0, v22=0.0; 
      double up1=0.0, up2=0.0; 

      int nu = 0;
      for (int j=0, jn=1; j<ntrio; j=jn, jn++) {
	int pj = proband[j] - 1;
	int fj = father[j] - 1;
	int mj = mother[j] - 1;
	if (nrules) {
	  double xp = ximp[pj];
	  double xf = ximp[fj];
	  double xm = ximp[mj];
	  if (!(ISNA(xp) || ISNA(xf) || ISNA(xm))) {
	    nu++;
	    up1 += xp - (xf+xm)/2.0;
	    up2 += (xp*xp - xf*xm)/4.0;
	  }
	}
	else {
	  unsigned char sp = snpsi[pj];
	  unsigned char sf = snpsi[fj];
	  unsigned char sm = snpsi[mj];
	  if (sp && sf && sm) {
	    sp--;
	    sf--;
	    sm--;
	    int jm = mendelian[sp + 3*sm + 9*sf];
	    if (!check || jm)  {
	      nu++;
	      int hetf = (sf==1);
	      int hetm = (sm==1);
	      if (hetf || hetm) { 
		int homp = (sp==2);
		int ss = sf + sm;
		if (!robust) {
		  u1 += (double) sp - (double) ss /2.0;
		  v11 += (double) (hetm + hetf)/4.0;
		  u2 += (double) homp - (double)(sf*sm)/4.0;
		  if (hetm) {
		    if (hetf) {
		      v22 += thrsix;
		      v12 += quart;
		    }
		    else if (sf==2) {
		      v22 += quart;
		      v12 += quart;
		    }
		  }
		  else if (sm==2) {
		    v22 += quart;
		    v12 += quart;
		  }
		}
		else {
		  up1 += (double) sp - (double) ss/2.0;
		  up2 += (double) homp - (double)(sf*sm)/4.0;
		}
	      }
	    }
	    else {
	      nonmend += !jm;
	    }
	  }
	}
	if (robust) {
	  if ((jn==ntrio) || (cluster[jn]!=cluster[j])) {
	    v11 += up1*up1;
	    v12 += up1*up2;
	    v22 += up2*up2;
	    u1 += up1;
	    u2 += up2;
	    up1 = up2 = 0.0;
	  }
	}
      }
      /* Store results in output object */
    
      Nused[t] = nu;
      if (nrules)
	nr2[t] = nu*r2;
      umat[t] = u1;
      umat[ntest+t] = u2;
      vmat[t] = v11;
      vmat[ntest+t] = v12;
      vmat[2*ntest+t] = v22;
    }
  }

  else {

    for (int t=0; t<ntest; t++) {
      int i = snp_subset? snp_subset[t] - 1: t; 
      const unsigned char *snpsi = NULL;
      double r2 = 0.0;
      if (nrules) {
	if (i >= nrules)
	  error("snp_subset out of range");
	SEXP Rule = VECTOR_ELT(Rules, i);
	if (isNull(Rule)) {
	  for (int j=0; j<nsubj; j++)
	    ximp[j] = 0.0;
	}
	else {
	  do_impute(snps, nsubj, NULL, nsubj, name_index, Rule, ximp, NULL);
	  r2 = *REAL(VECTOR_ELT(Rule, 0));
	}
      }
      else {
	if (i >= nsnp)
	  error("snp_subset out of range");
	snpsi = snps + nsubj*i;
      }

      double u = 0, v=0, u1=0.0, u2=0.0, v11=0.0, v12=0.0, v22=0.0; 
      double up=0.0, up1=0.0, up2=0.0; 
      
      int nu = 0;

      for (int j=0, jn=1; j<ntrio; j=jn, jn++) {
	int pj = proband[j] - 1;
	int fj = father[j] - 1;
	int mj = mother[j] - 1;
	int fp = female[pj];
	if (ISNA(fp)) 
	  continue;
	if (nrules) {
	  double xp = ximp[pj];
	  double xf = ximp[fj];
	  double xm = ximp[mj];
	  if (!(ISNA(xp) || ISNA(xf) || ISNA(xm))) {
	    nu++;
	    u +=  xp - (xf+xm)/2.0;
	    if (fp) {
	      up1 += xp - (xf+xm)/2.0;
	      up2 += (xp*xp - xf*xm)/4.0;
	    }
	  }
	}
	else {
	  int sp = (int) snpsi[pj];
	  int sf = (int) snpsi[fj];
	  int sm = (int) snpsi[mj];
	  if (sp && sf && sm) {
	    sp--;
	    sf--;
	    sm--;
	    int malehet = (sf==1) || (!fp && (sp==1));
            int jm = fp? mendelian[sp + 3*sm + 9*sf]: 
                         mendelian[sp + 3*sm + 27];
	    if (!malehet && (!check || jm))  {
	      nu++;
	      int hetm = (sm==1);
	      if (hetm) {
		int homp = (sp==2);
		int homf = (sf==2);
		double esp = fp? (double) 0.5 + homf: 1.0;
		double uw = (double) sp - esp;
		if (!robust) {
		  u += uw;
		  if (fp) {
		    u1 += uw;
		    v11 += quart; 
		    if (homf) {
		      u2 += (double) homp - 0.5;
		      v22 += quart;
		      v12 += quart;
		    }
		    v += quart;
		  }
		  else {
		    v += 1.0;
		  }
		}
		else {
		  up += uw;
		  if (fp) {
		    up1 += uw;
		    if (homf)
		      up2 += (double) homp - 0.5;
		  }
		}
	      }
	    }
	    else {
	      nonmend += !jm;
	      Xerrors += malehet;
	    }
	  }
	}
	if (robust) {
	  if ((jn==ntrio) || (cluster[jn]!=cluster[j])) {
	    v += up*up;
	    v11 += up1*up1;
	    v12 += up1*up2;
	    v22 += up2*up2;
	    u += up;
	    u1 += up1;
	    u2 += up2;
	    up = up1 = up2 = 0.0;
	  }
	}
      }
      Nused[t] = nu;
      if (nrules)
	nr2[t] = nu*r2;
      umat[t] = u;
      umat[ntest+t] = u1;
      umat[2*ntest+t] = u2;
      vmat[t] = v;
      vmat[ntest+t] = v11;
      vmat[2*ntest+t] = v12;
      vmat[3*ntest+t] = v22;
    }
  }
 
   /* Warnings */

  if (Xerrors) 
    warning("%d instances of a male coded as heterozygous at an X locus", 
	    Xerrors);
  if (nonmend)
    warning("%d misinheritances were detected", nonmend);

  /* Tidy up */

  index_destroy(name_index);
  if (nrules) {
    Free(ximp);
  }

  /* Attributes of output object */

  SEXP Names;
  PROTECT(Names = allocVector(STRSXP, 4));
  SET_STRING_ELT(Names, 0, mkChar("U"));
  SET_STRING_ELT(Names, 1, mkChar("V"));
  SET_STRING_ELT(Names, 2, mkChar("N"));
  SET_STRING_ELT(Names, 3, mkChar("N.r2"));

  setAttrib(Result, R_NamesSymbol, Names);
  UNPROTECT(6);

  return(Result);
}

