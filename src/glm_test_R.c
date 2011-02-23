#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <string.h>

#include "Rmissing.h"

#include "glm_test.h"
#include "hash_index.h"
#include "imputation.h"
#include "invert.h"

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
      error("Dimension error - Snp_subset");
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

  int S = 1;
  const int *stratum = NULL;
  if (TYPEOF(Stratum)==INTSXP) {
    if (LENGTH(Stratum)!=N) 
      error("Dimension error - Stratum"); 
    stratum = INTEGER(Stratum);
    for (int i=0; i<N; i++) 
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
    for (int i=0; i<N; i++)
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

  SEXP Result, Rnames;
  PROTECT(Result = allocVector(VECSXP, ntest));

  SEXP snpNames = VECTOR_ELT(getAttrib(Y, R_DimNamesSymbol), 1);
  SEXPTYPE sntype = TYPEOF(snpNames);
  PROTECT(Rnames = allocVector(sntype, ntest));

  /* Do calculations */

  for (int t=0; t<ntest; t++) {

    SEXP U, V;
    PROTECT(U = allocVector(REALSXP, P));
    PROTECT(V = allocVector(REALSXP, (P*(P+1))/2));
    double *u = REAL(U);
    double *v = REAL(V);
       
    int j = snp_subset? snp_subset[t] - 1: t; 
    const unsigned char *yj = y + N*j;
    int mono = 1, yv = 0;
 

    /* Load SNP as Binomial y-variate, with prior weights */

    if (ifX) {
      for (int i=0; i<N; i++) {
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
      for (int i=0; i<N; i++) {
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
      memset(u, 0x00, P*sizeof(double));
      memset(v, 0x00, sizeof(double)*(P*(P+1))/2); 
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
      
      if (dfr) {
	glm_score_test(N, rank, S, stratum, P, z, C, cluster,
		       resid, weights, xb, scale,
		       r2Max, u, v);
      }
      else {
	memset(u, 0x00, P*sizeof(double));
	memset(v, 0x00, sizeof(double)*(P*(P+1))/2);
      }
    }
    SEXP Score, Chi2, Df, Nunit;
    PROTECT(Score = allocVector(VECSXP, 6));
    PROTECT(Chi2 = allocVector(REALSXP, 1));
    PROTECT(Df = allocVector(INTSXP, 1));
    PROTECT(Nunit = allocVector(INTSXP, 1));
    int *nunit = INTEGER(Nunit);
    for (int i; i<N; i++) 
      if (weights[i]) (*nunit)++;
	    
    double *chi2 = REAL(Chi2);
    int *df = INTEGER(Df);
    if (P>1) {
      if (qform(P, u, v, NULL, chi2, df)) {
	warning("Matrix not positive semi-definite in test ", t+1);
	*chi2 = NA_REAL;
	*df = NA_INTEGER;
      }
    }
    else {
      if ((*v)>0) {
	*chi2 = (*u)*(*u)/(*v);
	*df = 1;
      }
      else {
	*chi2 = NA_REAL;
	*df = NA_INTEGER;
      }
    }
    SEXP Names;
    PROTECT(Names = allocVector(STRSXP, 6));
    SET_STRING_ELT(Names, 0, mkChar("parameters"));
    SET_STRING_ELT(Names, 1, mkChar("chi.squared"));
    SET_STRING_ELT(Names, 2, mkChar("df"));
    SET_STRING_ELT(Names, 3, mkChar("U"));
    SET_STRING_ELT(Names, 4, mkChar("V"));
    SET_STRING_ELT(Names, 5, mkChar("N"));
    setAttrib(Score, R_NamesSymbol, Names);
    SET_VECTOR_ELT(Score, 0, R_NilValue); /* Placeholder */
    SET_VECTOR_ELT(Score, 1, Chi2);
    SET_VECTOR_ELT(Score, 2, Df);
    SET_VECTOR_ELT(Score, 3, U);
    SET_VECTOR_ELT(Score, 4, V);
    SET_VECTOR_ELT(Score, 5, Nunit);
    SET_VECTOR_ELT(Result, t, Score);
    UNPROTECT(7);
    
    /* Name of list element */
    if (sntype==INTSXP) 
      INTEGER(Rnames)[t] = INTEGER(snpNames)[j];
    else 
      SET_STRING_ELT(Rnames, t, mkChar(CHAR(STRING_ELT(snpNames,j))));

  }
  setAttrib(Result, R_NamesSymbol,  Rnames);
  
  SEXP Class, Package;
  PROTECT(Class = allocVector(STRSXP, 1));
  SET_STRING_ELT(Class, 0, mkChar("snp.tests.glm"));
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("snpMatrix"));
  setAttrib(Class, install("package"), Package);
  classgets(Result, Class);
  SET_S4_OBJECT(Result);

  UNPROTECT(4);

  return(Result);
}


SEXP snp_rhs_score(SEXP Y, SEXP family, SEXP link, 
		   SEXP X, SEXP Stratum, SEXP Z, SEXP Rules, 
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

  int S = 1;
  int *stratum = NULL;
  if (TYPEOF(Stratum)==INTSXP) {
    if (LENGTH(Stratum)!=N) 
      error("Dimension error - Stratum"); 
    stratum = INTEGER(Stratum);
    for (int i=0; i<N; i++) 
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
  int ifX = 0;
  if (!strncmp(classZ, "snp", 3))
    ifX = 0;
  else if (!strncmp(classZ, "X.snp", 5))
    ifX = 1;
  else 
    error("Argument error - class(Z)");

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

  int *female;
  if (ifX) {
    SEXP Female = R_do_slot(Z, mkString("Female")); 
    female = LOGICAL(Female);
  }
  else
    female = NULL;

  /* If imputation involved, calculate snp name index */

  index_db name_index = NULL;
  SEXP Snp_names =  VECTOR_ELT(getAttrib(Z, R_DimNamesSymbol), 1);
  SEXP Rule_names;			    
  if (TYPEOF(Rules)!=NILSXP) {
    name_index = create_name_index(Snp_names);
    Rule_names = getAttrib(Rules, R_NamesSymbol);
  }

  
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
    for (int i=0; i<ntest; i++) {
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
    for (int i=0; i<N; i++)
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


  /* Output list */

  SEXP Result, TestNames;
  PROTECT(Result = allocVector(VECSXP, ntest));
  SEXP NameTests = getAttrib(Tests, R_NamesSymbol);
  int gen_names = TYPEOF(NameTests)==NILSXP;
  if (gen_names)
    PROTECT(TestNames =  allocVector(STRSXP, ntest));
	  
  /* Do tests */

  int snp;
  int *snps = &snp;
  for (int test=0; test<ntest; test++) {
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
      for (int i=0; i<N; i++) prior[i] = 1.0;
    else 
      for (int i=0; i<N; i++) prior[i] = prior_all[i];
    
    /* Names of SNPs tested */

    SEXP Tnames;
    PROTECT(Tnames = allocVector(STRSXP, nsnpt));

    /* Set up Z matrix, tracking incomplete cases  */

    int missing = 0, err=0;
    for (int j=0, ij=0; j<nsnpt; j++) {
      int len;
      if (gen_names) {
	if (j){
	  len = strlen(testname);
	  if (len<max_name_length) {
	    testname[len] = '+';
	    len++;
	  }
	  else
	    len=0;
	}
      }
      int space = max_name_length - len;
      int snpsj = snps[j];
      if (snpsj>0) {
	snpsj--;
	SEXP Snp_namej =  STRING_ELT(Snp_names, snpsj);
	SET_STRING_ELT(Tnames, j, Snp_namej);
	if (gen_names && space>0)
	  strncpy(testname+len, CHAR(Snp_namej), space);
	const unsigned char *zj = z + N*snpsj;
	for (int i=0; i<N; i++, ij++) {
	  unsigned char zij = zj[i];
	  if (zij) {
	    zw[ij] = (double) (zij - 1);
	  }
	  else {
	    zw[ij] = 0.0;
	    if (prior[i]) {
	      missing ++;
	      prior[i] = 0.0;
	    }
	  }
	}
      }
      else {
	snpsj = -(1+snpsj);
	SEXP Rule_namej = STRING_ELT(Rule_names, snpsj);
	SET_STRING_ELT(Tnames, j, Rule_namej);
	if (gen_names && space>0)
	  strncpy(testname+len, CHAR(Rule_namej), space);	
	SEXP Rule =  VECTOR_ELT(Rules, snpsj);
	if (!isNull(Rule)){ /* Not monomorphic */
	  do_impute(z, N, NULL, N, name_index, Rule, zw+ij, NULL);
	  for (int i=0; i<N; i++, ij++) {
	    if (ISNA(zw[ij])) {
	      zw[ij] = 0.0;
	      if (prior[i]) {
		missing ++;
		prior[i] = 0.0;
	      }
	    }
	  }   
	}
	else {
	  for (int i=0; i<N; i++, ij++)
	    zw[ij] = 0.0;
	}
      }
    }
    if (gen_names)
      SET_STRING_ELT(TestNames, test, mkChar(testname));
   	
    SEXP U, V;
    PROTECT(U = allocVector(REALSXP, nsnpt));
    PROTECT(V = allocVector(REALSXP, (nsnpt*(nsnpt+1))/2));
    double *u = REAL(U);
    double *v = REAL(V);

    if (missing == N) { 

      /* No data */
      
      warning("No valid data for test %d", test+1);
      memset(u, 0x00, nsnpt*sizeof(double));
      memset(v, 0x00, sizeof(double)*(nsnpt*(nsnpt+1))/2);

    }
    else if (!nsnpt) {
      
      /* No polymorphic data */

      warning("No polymorphic markers for test %d", test+1);
      memset(u, 0x00, nsnpt*sizeof(double));
      memset(v, 0x00, sizeof(double)*(nsnpt*(nsnpt+1))/2);
 
    }
    else if (missing > mallowed) {

      /* Refit the model using current fitted values as start */

      for (int i=0; i<N; i++) 
	fitted[i] = fitted_all[i];
      int err = glm_fit(fam, lnk, N, M, S, y, prior, x, stratum, 
			maxit, epsilon, 1,
			&rank, xb, fitted, resid, weights, 
			&scale, &df_r);
      if (err) 
	warning("No convergence while fitting base model for test %d", test+1);
      glm_score_test(N, rank, S, stratum, nsnpt, zw, C, cluster,
		     resid, weights, xb, scale,
		     r2Max, u, v);
 
    }
    else {

      /* Use  GLM fit to total data */

      df_r = df_r_all - missing;
      rank = rank_all;
      for (int i=0; i<N; i++) 
	weights[i] = prior[i]? weights_all[i]: 0.0;
      glm_score_test(N, rank, S, stratum, nsnpt, zw, C, cluster,
		     resid_all, weights, xb_all, scale_all,
		     r2Max, u, v);
    }
    SEXP Score, Chi2, Df, Nused;
    PROTECT(Score = allocVector(VECSXP, 6));
    PROTECT(Chi2 = allocVector(REALSXP, 1));
    PROTECT(Df = allocVector(INTSXP, 1));
    PROTECT(Nused = allocVector(INTSXP, 1));
    int *nused = INTEGER(Nused);
    *nused = 0;
    for (int i=0; i<N; i++)
      if(weights[i]) (*nused)++;
    double *chi2 = REAL(Chi2);
    int *df = INTEGER(Df);
    if (nsnpt>1) {
      if (qform(nsnpt, u, v, NULL, chi2, df)) {
	warning("Matrix not positive semi-definite in test ", test+1);
	*chi2 = NA_REAL;
	*df = NA_INTEGER;
      }
    }
    else {
      if ((*v)>0.0) {
	*chi2 = (*u)*(*u)/(*v);
	*df = 1;
      }
      else {
	*chi2 = NA_REAL;
	*df = NA_INTEGER;
      }
    }
    SEXP Names;
    PROTECT(Names = allocVector(STRSXP, 6));
    SET_STRING_ELT(Names, 0, mkChar("parameters"));
    SET_STRING_ELT(Names, 1, mkChar("chi.squared"));
    SET_STRING_ELT(Names, 2, mkChar("df"));
    SET_STRING_ELT(Names, 3, mkChar("U"));
    SET_STRING_ELT(Names, 4, mkChar("V"));
    SET_STRING_ELT(Names, 5, mkChar("N"));
    setAttrib(Score, R_NamesSymbol, Names);
    SET_VECTOR_ELT(Score, 0, Tnames);
    SET_VECTOR_ELT(Score, 1, Chi2);
    SET_VECTOR_ELT(Score, 2, Df);
    SET_VECTOR_ELT(Score, 3, U);
    SET_VECTOR_ELT(Score, 4, V);
    SET_VECTOR_ELT(Score, 5, Nused);
    SET_VECTOR_ELT(Result, test, Score);
    setAttrib(Result, R_NamesSymbol, Names);
    UNPROTECT(8);
  }

  /* Return hash table memory */

  if (name_index)
    index_destroy(name_index);

  SEXP Class, Package;
  PROTECT(Class = allocVector(STRSXP, 1));
  SET_STRING_ELT(Class, 0, mkChar("snp.tests.glm"));
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("snpMatrix"));
  setAttrib(Class, install("package"), Package);
  classgets(Result, Class);
  SET_S4_OBJECT(Result);

  if (gen_names) {
    setAttrib(Result, R_NamesSymbol,  TestNames);
    UNPROTECT(4);
  }
  else {
    setAttrib(Result, R_NamesSymbol, NameTests);
    UNPROTECT(3);
  }

  return(Result);
  
}

/* Pooling two objects */

SEXP pool2_glm(SEXP X, SEXP Y) {
  
  int N = LENGTH(X);
  if (N!=LENGTH(Y))
    error("pool2_glm: unequal length arguments");
  SEXP Result;
  PROTECT(Result = allocVector(VECSXP, N));
  for (int i=0; i<N; i++) {
    SEXP Xi = VECTOR_ELT(X, i);
    SEXP Yi = VECTOR_ELT(Y, i);
    SEXP XiU = VECTOR_ELT(Xi, 3);
    double *xiu = REAL(XiU);
    SEXP XiV = VECTOR_ELT(Xi, 4);
    double *xiv = REAL(XiV);
    SEXP YiU = VECTOR_ELT(Yi, 3);
    double *yiu = REAL(YiU);
    SEXP YiV = VECTOR_ELT(Yi, 4);
    double *yiv = REAL(YiV);
    int nu = LENGTH(XiU);
    int nv = LENGTH(XiV);
    if (LENGTH(YiU)!=nu)
      error("attempt to pool tests on unequal numbers of parameters");
    SEXP RU, RV, Rchi2, Rdf, Rnused;
    PROTECT(RU = allocVector(REALSXP, nu));
    double *ru = REAL(RU);
    PROTECT(RV = allocVector(REALSXP, nv));
    double *rv = REAL(RV);
    memset(ru, 0x00, nu*sizeof(double));
    memset(rv, 0x00, nv*sizeof(double));
    for (int j=0; j<nu; j++) 
      ru[j] = xiu[j] + yiu[j];
    for (int j=0; j<nv; j++) 
      rv[j] = xiv[j] + yiv[j];
    PROTECT(Rchi2 = allocVector(REALSXP, 1));
    PROTECT(Rdf = allocVector(INTSXP, 1));
    double *chi2 = REAL(Rchi2);
    int *df = INTEGER(Rdf);
    if (qform(nu, ru, rv, NULL, chi2, df)) {
      warning("Matrix not positive semi-definite in pooled test ", i+1);
      *chi2 = NA_REAL;
      *df = NA_INTEGER;
    }
    PROTECT(Rnused = allocVector(INTSXP, 1));
    *(INTEGER(Rnused)) = *(INTEGER(VECTOR_ELT(X, 5))) + 
			 *(INTEGER(VECTOR_ELT(Y, 5)));
    SEXP Resi, Names;
    PROTECT(Resi = allocVector(VECSXP, 6));
    PROTECT(Names = allocVector(STRSXP, 6));
    SET_STRING_ELT(Names, 0, mkChar("parameters"));
    SET_STRING_ELT(Names, 1, mkChar("chi.squared"));
    SET_STRING_ELT(Names, 2, mkChar("df"));
    SET_STRING_ELT(Names, 3, mkChar("U"));
    SET_STRING_ELT(Names, 4, mkChar("V"));
    SET_STRING_ELT(Names, 5, mkChar("N"));
    setAttrib(Resi, R_NamesSymbol, Names);
    SET_VECTOR_ELT(Resi, 0, VECTOR_ELT(Xi, 0));
    SET_VECTOR_ELT(Resi, 1, Rchi2);
    SET_VECTOR_ELT(Resi, 2, Rdf);
    SET_VECTOR_ELT(Resi, 3, RU);
    SET_VECTOR_ELT(Resi, 4, RV);
    SET_VECTOR_ELT(Resi, 5, Rnused);
    SET_VECTOR_ELT(Result, i, Resi);
    UNPROTECT(7);
  }

  UNPROTECT(1);

  return(Result);
}

    
