#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <string.h>

#include "Rmissing.h"

#include "glm_test.h"

#define MAX_NAME_LENGTH 81

SEXP snp_lhs_score(const SEXP Y, const SEXP X, const SEXP Stratum, 
		   const SEXP Z, const SEXP Snp_subset, 
		   const SEXP Robust, const SEXP Cluster, const SEXP Control) {

  /* Y should be a snp.matrix or an X.snp.matrix */
  const char *classY = NULL;
  if (TYPEOF(R_data_class(Y, FALSE)) == STRSXP) {
    classY = CHAR(STRING_ELT(R_data_class(Y, FALSE), 0));
  } else {
    classY = CHAR(STRING_ELT(getAttrib(Y, R_ClassSymbol), 0));
  }
  if(!IS_S4_OBJECT(Y)) {
    error("Y in snp_lhs_score is missing S4 Object bit");
  }
  int ifX = 0;
  if (!strncmp(classY, "snp", 3))
    ifX = 0;
  else if (!strncmp(classY, "X.snp", 5))
    ifX = 1;
  else 
    error("Argument error - class(Y)");

  /* Y and its dimensions */

  if (TYPEOF(Y)!=RAWSXP)
    error("Argument error - Y");
  const unsigned char *y = RAW(Y);
  int *dim;
  int N, nsnp;
  if (strlen(classY)>5) {
    dim = INTEGER(getAttrib(Y, R_DimSymbol));
    N = dim[0];
    nsnp = dim[1];
  }
  else {
    N = LENGTH(Y);
    nsnp = 1;
  }

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


  /* If X-chromosome, indicators of female sex */

  const int *female;
  if (ifX) {
    SEXP Female = R_do_slot(Y, mkString("Female")); 
    female = LOGICAL(Female);
  }
  else
    female = NULL;
  
  /* X and its dimensions */

  double *x = NULL;
  int M = 0;
  if (TYPEOF(X)==REALSXP) {
    x = REAL(X);
    const int *dim = INTEGER(getAttrib(X, R_DimSymbol));
    M = dim[1];
    if (dim[0]!=N)
      error("Dimension error - X"); 
  }
  else if (TYPEOF(X)!=NILSXP)
    error("Argument error - X");

  /* Stratum */

  int i=0, t=0;
  int S = 1;
  const int *stratum = NULL;
  if (TYPEOF(Stratum)==INTSXP) {
    if (LENGTH(Stratum)!=N) 
      error("Dimension error - Stratum"); 
    stratum = INTEGER(Stratum);
    for (i=0; i<N; i++) 
      if (stratum[i]>S) S = stratum[i];
  }
  else if (TYPEOF(Stratum)!=NILSXP)
    error("Argument error - Stratum");

  /* Z and its dimensions */

  if (TYPEOF(Z)!=REALSXP)
    error("Argument error - Z");
  double *z = REAL(Z);
  dim = INTEGER(getAttrib(Z, R_DimSymbol));
  int P = dim[1];
  if (dim[0]!=N)
    error("Dimension error - Z"); 


  /* Robust */

  if (TYPEOF(Robust)!=LGLSXP)
    error("Argument error - Robust");
  int robust = LOGICAL(Robust)[0];

  /* Cluster */

  int C = robust? 1: 0;
  int *cluster = NULL;
  if (TYPEOF(Cluster)==INTSXP) {
    if (LENGTH(Cluster)!=N)
      error("Dimension error - Cluster"); 
    cluster = INTEGER(Cluster);
    for (i=0; i<N; i++)
      if (cluster[i]>C) C = cluster[i];
  }
  else if (TYPEOF(Cluster)!=NILSXP)
    error("Argument error - Cluster");
  
  /* Control object */

  if (TYPEOF(Control)!=VECSXP || LENGTH(Control)!=3) 
    error("Argument error - Control");
  SEXP Maxit = VECTOR_ELT(Control, 0);
  if (TYPEOF(Maxit)!=INTSXP || length(Maxit)!=1)
    error("Argument error - Maxit");
  int maxit = INTEGER(Maxit)[0];
  SEXP Epsilon = VECTOR_ELT(Control, 1);
  if (TYPEOF(Epsilon)!=REALSXP || LENGTH(Epsilon)!=1)
    error("Argument error - Epsilon");
  double epsilon = REAL(Epsilon)[0];
  SEXP R2Max = VECTOR_ELT(Control, 2);
  if (TYPEOF(R2Max)!=REALSXP || LENGTH(R2Max)!=1)
    error("Argument error - R2Max");
  double r2Max = REAL(R2Max)[0];
  
  /* Work arrays */

  double *yd = (double *) R_alloc(N, sizeof(double));
  double *fitted = (double *) R_alloc(N, sizeof(double));
  double *prior = (double *) R_alloc(N, sizeof(double));
  double *resid = (double *) R_alloc(N, sizeof(double));
  double *weights = (double *) R_alloc(N, sizeof(double));
  double *xb;
  if (x)
    xb = (double *) R_alloc(N*M, sizeof(double));
  else 
    xb = NULL;

  /* Output arrays */

  SEXP Result, Chi2, Df, Df_r;
  PROTECT(Result = allocVector(VECSXP, 3));
  PROTECT(Chi2 = allocVector(REALSXP, ntest));
  SET_VECTOR_ELT(Result, 0, Chi2); /* Chi-squared */
  double *chi2 = REAL(Chi2); 
  PROTECT(Df = allocVector(INTSXP, ntest));
  SET_VECTOR_ELT(Result, 1, Df); /* Df */
  int *df = INTEGER(Df);
  PROTECT(Df_r = allocVector(INTSXP, ntest));
  SET_VECTOR_ELT(Result, 2, Df_r); /* Df.resid */
  int *df_resid = INTEGER(Df_r);


  /* Do calculations */

  for (t=0; t<ntest; t++) {
    int j = snp_subset? snp_subset[t] - 1: t; 
    const unsigned char *yj = y + N*j;
    int mono = 1, yv = 0;
 

    /* Load SNP as Binomial y-variate, with prior weights */

    if (ifX) {
      for (i=0; i<N; i++) {
	int yij = (int) yj[i];
	if (yij) {
	  if (!yv)
	    yv = yij;
	  else if (mono) 
	    mono = (yv == yij);
	  prior[i] = female[i]? 2.0: 1.0;
	  yd[i] = ((double) (yij - 1))/2.0;
	}
	else {
	  prior[i] = yd[i] = 0.0;
	}
      }
    }
    else {
      for (i=0; i<N; i++) {
	int yij = (int) yj[i];
	if (yij) {
	  if (!yv)
	    yv = yij;
	  else if (mono) 
	    mono = (yv == yij);
	  prior[i] = 2.0;
	  yd[i] = ((double) (yij - 1))/2.0;
	}
	else {
	  prior[i] = yd[i] = 0.0;
	}	
      }
    }
    if (mono) { /* Monomorphic SNP */
      chi2[t] = NA_REAL;
      df[t] = NA_INTEGER;
      df_resid[t] = NA_INTEGER;
    }
    else {
      
      /* Fit base model */
    
      int rank, dfr;
      double scale;
      double rc = glm_fit(BINOMIAL, LOGIT, N, M, S, yd, prior, x, stratum, 
			  maxit, epsilon, 0, 
			  &rank, xb, fitted, resid, weights, &scale, &dfr);
      if (rc) 
	warning("Failure to converge while fitting base model for SNP %d",j+1);
    
      /* Score test */
      
      df_resid[t] = dfr;
      if (dfr) {
	int dfj;
	double chi2j;
	glm_score_test(N, rank, S, stratum, P, z, C, cluster,
		       resid, weights, xb, scale,
		       r2Max, &chi2j, &dfj);
	if (dfj) {
	  chi2[t] = chi2j;
	  df[t] = dfj;
	}
	else {
	  chi2[t] = NA_REAL;
	  df[t] = NA_INTEGER;
	}
      }
      else {
	chi2[t] = NA_REAL;
	df[t] = NA_INTEGER;
      }
    }
  }
  
/* Attributes of output object */

  SEXP cNames, rnames, dfClass;
  SEXP snpNames = VECTOR_ELT(getAttrib(Y, R_DimNamesSymbol), 1);
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

  PROTECT(cNames = allocVector(STRSXP, 3));
  SET_STRING_ELT(cNames, 0, mkChar("Chi.squared"));
  SET_STRING_ELT(cNames, 1, mkChar("Df"));
  SET_STRING_ELT(cNames, 2, mkChar("Df.residual"));
  setAttrib(Result, R_NamesSymbol, cNames);
  PROTECT(dfClass = allocVector(STRSXP, 1));
  SET_STRING_ELT(dfClass, 0, mkChar("data.frame"));
  setAttrib(Result, R_ClassSymbol, dfClass);
  UNPROTECT(7);

  return(Result);
}


SEXP snp_rhs_score(SEXP Y, SEXP family, SEXP link, 
		   SEXP X, SEXP Stratum, SEXP Z, 
		   SEXP Prior, SEXP Tests, SEXP Robust, SEXP Cluster, 
		   SEXP Control, SEXP MissAllow) {

  char testname[MAX_NAME_LENGTH];
  int max_name_length =  MAX_NAME_LENGTH -1;
  testname[max_name_length] = (char) 0; /* Make sure its null terminated */

  /* Y and its dimensions */

  if (TYPEOF(Y)!=REALSXP)
    error("Argument error - Y");
  double *y = REAL(Y);
  int N = LENGTH(Y);

  /* Family */

  if (TYPEOF(family)!=INTSXP || LENGTH(family)!=1)
    error("Argument error - family");
  int fam = INTEGER(family)[0];
  if (fam<1 || fam>4)
    error("Illegal family argument");
  
  /* Link */

  if (TYPEOF(link)!=INTSXP || LENGTH(link)!=1)
    error("Argument error - link");
  int lnk = INTEGER(link)[0];
  if (lnk<1 || lnk>4)
    error("Illegal link argument");

  /* X and its dimensions */

  double *x = NULL;
  int M = 0;
  int *dim = NULL;
  if (TYPEOF(X)==REALSXP) {
    x = REAL(X);
    dim = INTEGER(getAttrib(X, R_DimSymbol));
    M = dim[1];
    if (dim[0]!=N)
      error("Dimension error - X"); 
  }
  else if (TYPEOF(X)==NILSXP) {
    M = 0; 
    x = NULL;
  }
  else
    error("Argument error - X");

  /* Stratum */

  int i=0, j=0, ij=0, test=0;
  int S = 1;
  int *stratum = NULL;
  if (TYPEOF(Stratum)==INTSXP) {
    if (LENGTH(Stratum)!=N) 
      error("Dimension error - Stratum"); 
    stratum = INTEGER(Stratum);
    for (i=0; i<N; i++) 
      if (stratum[i]>S) S = stratum[i];
  }
  else if (TYPEOF(Stratum)!=NILSXP)
    error("Argument error - Stratum");

  /* Z should be a snp.matrix or an X.snp.matrix */

  const char *classZ = NULL;
  if (TYPEOF(R_data_class(Z, FALSE)) == STRSXP) {
    classZ = CHAR(STRING_ELT(R_data_class(Z, FALSE), 0));
  } else {
    classZ = CHAR(STRING_ELT(getAttrib(Z, R_ClassSymbol), 0));
  }
  if(!IS_S4_OBJECT(Z)) {
    error("Z in snp_rhs_score is missing S4 Object bit");
  }
  /*
  int ifX = 0;
  if (!strncmp(classZ, "snp", 3))
    ifX = 0;
  else if (!strncmp(classZ, "X.snp", 5))
    ifX = 1;
  else 
    error("Argument error - class(Z)");
  */

  /* Z and its dimensions */

  int Nz;
  if (TYPEOF(Z)!=RAWSXP)
    error("Argument error - Z");
  const unsigned char *z = RAW(Z);
  int nsnp;
  if (strlen(classZ)>5) {
    dim = INTEGER(getAttrib(Z, R_DimSymbol));
    Nz = dim[0];
    nsnp = dim[1];
  }
  else {
    Nz = LENGTH(Z);
    nsnp = 1;
  }
  if (Nz!=N)
    error("Dimension error - Z");

  /* If X-chromosome, indicators of female sex */

  /*
  int *female;
  if (ifX) {
    SEXP Female = R_do_slot(Z, mkString("Female")); 
    female = LOGICAL(Female);
  }
  else
    female = NULL;
  */
  
  /* Prior weights */

  double *prior_all = NULL;
  if (TYPEOF(Prior)==REALSXP) {
    if (LENGTH(Prior)!=N)
      error("Dimension error - Prior");
    prior_all = REAL(Prior);
  }
  else if (TYPEOF(Prior)!=NILSXP)
    error("Argument error - Prior");

  /* Tests. For now snps referred to by col number not name */

  int test_size = 1;
  int *test_int = NULL;
  int ntest = nsnp;
  SEXPTYPE tests_type = TYPEOF(Tests);
  if (tests_type==VECSXP) {
    ntest = LENGTH(Tests);
    for (i=0; i<ntest; i++) {
      SEXP testi =  VECTOR_ELT(Tests, i);
      if (TYPEOF(testi)!=INTSXP)
	error("Non-integer list - Tests[i]");
      int leni = LENGTH(testi);
      if (leni>test_size)
	test_size = leni;
    }
    if (test_size==1)
      error("List of multi-SNP tests doesn't contain a multi-SNP spec");
  }
  else if (tests_type==INTSXP) {
    test_int = INTEGER(Tests);
    ntest = LENGTH(Tests);
  }
  else if (tests_type!=NILSXP)
    error("Argument error - tests");

  /* Robust */

  if (TYPEOF(Robust)!=LGLSXP)
    error("Argument error - Robust");
  int robust = LOGICAL(Robust)[0];

  /* Cluster */

  int C = robust? 1: 0;
  int *cluster = NULL;
  if (TYPEOF(Cluster)==INTSXP) {
    if (LENGTH(Cluster)!=N)
      error("Dimension error - Cluster"); 
    cluster = INTEGER(Cluster);
    for (i=0; i<N; i++)
      if (cluster[i]>C) C = cluster[i];
  }
  else if (TYPEOF(Cluster)!=NILSXP)
    error("Argument error - Cluster");
  
  /* Control object */

  if (TYPEOF(Control)!=VECSXP || LENGTH(Control)!=3) 
    error("Argument error - Control");
  SEXP Maxit = VECTOR_ELT(Control, 0);
  if (TYPEOF(Maxit)!=INTSXP || length(Maxit)!=1)
    error("Argument error - Maxit");
  int maxit = INTEGER(Maxit)[0];
  SEXP Epsilon = VECTOR_ELT(Control, 1);
  if (TYPEOF(Epsilon)!=REALSXP || LENGTH(Epsilon)!=1)
    error("Argument error - Epsilon");
  double epsilon = REAL(Epsilon)[0];
  SEXP R2Max = VECTOR_ELT(Control, 2);
  if (TYPEOF(R2Max)!=REALSXP || LENGTH(R2Max)!=1)
    error("Argument error - R2Max");
  double r2Max = REAL(R2Max)[0];
  
  /* Allowable missingness */

  if (TYPEOF(MissAllow)!=REALSXP || LENGTH(MissAllow)!=1)
    error("Argument error - MissAllow");
  int mallowed = N * REAL(MissAllow)[0];

  /* Work arrays */

  double *zw = (double *) R_alloc(N*test_size, sizeof(double));
  double *prior = (double *) R_alloc(N, sizeof(double));
  double *resid = (double *) R_alloc(N, sizeof(double));
  double *resid_all = (double *) R_alloc(N, sizeof(double));
  double *weights = (double *) R_alloc(N, sizeof(double));
  double *weights_all = (double *) R_alloc(N, sizeof(double));
  double *fitted = (double *) R_alloc(N, sizeof(double));
  double *fitted_all = (double *) R_alloc(N, sizeof(double));
  double *xb, *xb_all;
  if (x){
    xb = (double *) R_alloc(N*M, sizeof(double));
    xb_all = (double *) R_alloc(N*M, sizeof(double));
  }
  else 
    xb =xb_all =  NULL;
  
  /* Fit base model to full data */

  int rank, rank_all, df_r, df_r_all;
  double scale, scale_all;
  int err = glm_fit(fam, lnk, N, M, S, y, prior_all, x, stratum, 
		    maxit, epsilon, 0, 
		    &rank_all, xb_all, 
		    fitted_all, resid_all, weights_all, 
		    &scale_all, &df_r_all);
  if (err) 
    error("failure to converge at initial fitting of base model");


  /* Output arrays */

  SEXP Result, Chi2, Df, Df_r, TestNames;
  PROTECT(Result = allocVector(VECSXP, 3));
  PROTECT(Chi2 = allocVector(REALSXP, ntest));
  SET_VECTOR_ELT(Result, 0, Chi2); /* Chi-squared */
  double *chi2 = REAL(Chi2); 
  PROTECT(Df = allocVector(INTSXP, ntest));
  SET_VECTOR_ELT(Result, 1, Df); /* Df */
  int *df = INTEGER(Df);
  PROTECT(Df_r = allocVector(INTSXP, ntest));
  SET_VECTOR_ELT(Result, 2, Df_r); /* Df.resid */
  int *df_resid = INTEGER(Df_r);
  
  SEXP snpNames = VECTOR_ELT(getAttrib(Z, R_DimNamesSymbol), 1);
  PROTECT(TestNames = allocVector(STRSXP, ntest));
    
	  
  /* Do tests */

  int snp;
  int *snps = &snp;
  for (test=0; test<ntest; test++) {
    int nsnpt = 1;
    if (test_size==1) { /* Single SNP tests */
      if (test_int)  /* Selected list */
	snps = test_int+test;
      else
	snp = test + 1;
    }
    else  { /* List of multi-SNP tests */
      SEXP Snps = VECTOR_ELT(Tests, test);
      nsnpt = LENGTH(Snps);
      snps = INTEGER(Snps);
    }

    if (!prior_all) 
      for (i=0; i<N; i++) prior[i] = 1.0;
    else 
      for (i=0; i<N; i++) prior[i] = prior_all[i];
    
    /* Set up Z matrix, tracking incomplete cases  */

    int missing = 0, /* err=0, */ nsnptu = 0;
    for (j=0, ij=0; j<nsnpt; j++) {
      int mono = 1;
      int snpsj = snps[j] - 1;
      unsigned char zv = 0;
      const unsigned char *zj = z + N*snpsj;
      for (i=0; i<N; i++, ij++) {
	unsigned char zij = zj[i];
	if (zij) {
	  zw[ij] = (double) (zij - 1);
	  if (!zv)
	    zv = zij;
	  else if (mono)
	    mono = (zv==zij);
	}
	else {
	  zw[ij] = 0.0;
	  if (prior[i]) {
	    missing ++;
	    prior[i] = 0.0;
	  }
	}
      }
      if (mono) 
	ij -= N;
      else 
	nsnptu++;
      int len;
      const char *snpnm = CHAR(STRING_ELT(snpNames, snpsj));
      if (j) {
	len = strlen(testname);
	if (len<max_name_length) {
	  testname[len] = '+';
	  len++;
	}
      }
      else
	len = 0;
      int space = MAX_NAME_LENGTH - 1 - len;
      if (space>0)
	strncpy(testname+len, snpnm, space);
    }
    SET_STRING_ELT(TestNames, test, mkChar(testname));
   
    double chi2t;
    int dft;
    if (missing == N) { 

      /* No data */
      
      warning("No valid data for test %d", test+1);
      err = 1;
      dft = df_r = 0;
    }
    else if (!nsnpt) {
      
      /* No polymorphic data */

      warning("No polymorphic markers for test %d", test+1);
      err = 1;
      dft = df_r = 0;
    }
    else if (missing > mallowed) {

      /* Refit the model using current fitted values as start */

      for (i=0; i<N; i++) 
	fitted[i] = fitted_all[i];
      int err = glm_fit(fam, lnk, N, M, S, y, prior, x, stratum, 
			maxit, epsilon, 1,
			&rank, xb, fitted, resid, weights, 
			&scale, &df_r);
      if (err) 
	warning("No convergence while fitting base model for test %d", test+1);
      glm_score_test(N, rank, S, stratum, nsnptu, zw, C, cluster,
		     resid, weights, xb, scale,
		     r2Max, &chi2t, &dft);
 
    }
    else {

      /* Use  GLM fit to total data */

      df_r = df_r_all - missing;
      rank = rank_all;
      for (i=0; i<N; i++) 
	weights[i] = prior[i]? weights_all[i]: 0.0;
      glm_score_test(N, rank, S, stratum, nsnptu, zw, C, cluster,
		     resid_all, weights, xb_all, scale_all,
		     r2Max, &chi2t, &dft);
    }
    df_resid[test] = df_r;
    if (dft) {
      chi2[test] = chi2t;
      df[test] = dft;
    }
    else {
      chi2[test] = NA_REAL;
      df[test] = NA_INTEGER;
    }
  }
  
  /* Attributes of output object */

  SEXP cNames, dfClass;
  setAttrib(Result, R_RowNamesSymbol,  TestNames);
  PROTECT(cNames = allocVector(STRSXP, 3));
  SET_STRING_ELT(cNames, 0, mkChar("Chi.squared"));
  SET_STRING_ELT(cNames, 1, mkChar("Df"));
  SET_STRING_ELT(cNames, 2, mkChar("Df.residual"));
  setAttrib(Result, R_NamesSymbol, cNames);
  PROTECT(dfClass = allocVector(STRSXP, 1));
  SET_STRING_ELT(dfClass, 0, mkChar("data.frame"));
  setAttrib(Result, R_ClassSymbol, dfClass);
  UNPROTECT(7);

  return(Result);
}
