#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "Rmissing.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "covwin.h"
#include "hash_index.h"
#include "imputation.h"

int bin_search(const double *sorted, const int len, const double value);
int nearest_N(const double *sorted, const int len, const double value, 
	      const int N);
double covariances(int i, int j, va_list ap);
double snpcov(const unsigned char *x, const unsigned char *y, 
	      const int *female, const int N);
double snpmean(const unsigned char *x, 
	       const int *female, const int N);
void utinv(double *, const int);


SEXP snp_impute(const SEXP X, const SEXP Y, const SEXP Xord, const SEXP Yord,
		const SEXP Xpos, const SEXP Ypos, const SEXP Phase, 
		const SEXP Try, const SEXP R2stop, const SEXP MaX){

  int pmax = *INTEGER(MaX);   /* Maximum number of predictor variables */
  int try =  *INTEGER(Try);   /* Number to search */
  double r2stop = *REAL(R2stop); /* R^2 value to stop inclusion */
  if (r2stop>1.0)
    r2stop = 1.0;
  int phase = *LOGICAL(Phase);   /* haploid or diploid computation */
  
  const double *xpos = REAL(Xpos);  /* Sorted list of X positions */
  int nx = LENGTH(Xpos);      /* Number of X s */
  int *xord = INTEGER(Xord);  /* Corresponding columns in X */
  const double *ypos = REAL(Ypos);  /* Sorted list of Y positions  */
  int ny = LENGTH(Ypos);      /* Number of Y s */
  int *yord = INTEGER(Yord);  /* Corresponding columns in Y */
  int nsubject = nrows(X);
  unsigned char *x = RAW(X);
  unsigned char *y = RAW(Y);
  SEXP Xsnpnames = VECTOR_ELT(getAttrib(X, R_DimNamesSymbol), 1); 
  int *female = NULL;
  SEXP cl = GET_CLASS(X);
  if (TYPEOF(cl) != STRSXP) {
    cl = R_data_class(X, FALSE); /* S4 way of getting class attribute */
  }
  if (!strcmp(CHAR(STRING_ELT(cl, 0)), "X.snp.matrix")) {
    SEXP Female = R_do_slot(X, mkString("Female"));
    female = LOGICAL(Female);
  }

  /* Work arrays */

  double *xy = (double *) Calloc(try, double);       /* XY covariances */
  double *xxd = (double *) Calloc(try, double);      /* Diagonals of XX */
  double *xxi = (double *) Calloc(try*pmax, double); /* Row of XX covars */
  int *sel = (int *) Calloc(pmax, int);              /* Selected SNPs */
  double *coef =(double *) Calloc((pmax*(pmax+1))/2, double);/* Coefficients*/ 
  double *ycoef = (double *) Calloc(pmax, double);   /* Y coefficients */
  COV_WIN_PTR cache = new_window(try, 0);            /* Covariance cache */

  /* Result */
 
  SEXP Result;
  PROTECT(Result = allocVector(VECSXP, ny));
  setAttrib(Result, R_NamesSymbol, 
	    VECTOR_ELT(getAttrib(Y, R_DimNamesSymbol),1));

  /* Main loop */

  for (int i=0; i<ny; i++) {
    unsigned char *yi = y + nsubject*(yord[i]-1);
    int start = nearest_N(xpos, nx, ypos[i], try);
    double yy = snpcov(yi, yi, female, nsubject);
    if (yy>0) {
      for (int j=0; j<try; j++) { 
	int jx = nsubject*(xord[start+j]-1);
	xy[j] = snpcov(x+jx, yi, female, nsubject);
      }
      move_window(cache, start);
      get_diag(cache, xxd, covariances, x, nsubject, xord, female, phase);
      int nregr = 0;
      double resid = yy;
      double rsq = 0.0;
      int ic = 0;
      while (1) {
	/* Find next snp to be included */
	double max_due = 0.0;
	int best = -1;
	for (int j=0; j<try; j++) {
	  double xxj = xxd[j];
	  if (xxj>0.0) {
	    double xyj = xy[j];
	    double xyj2 = xyj*xyj;
	    if (xyj2>(xxj*resid)) {
	      xy[j] = xyj>0.0? sqrt(xxj*resid): -sqrt(xxj*resid);
	      best = j;
	      max_due = resid;
	    }
	    else {
	      double due = xyj2/xxj;
	      if (due>max_due) {
		max_due = due;
		best = j;
	      }
	    }
	  }
	}
	sel[nregr] = best; /* Save index */
	double bestc = xy[best]/xxd[best]; 
	ycoef[nregr] = bestc; /* Save regression coefficient */
	double dX2 = (double) (nsubject-nregr-1) * max_due/resid;
	resid -= max_due;
	rsq = 1.0 - resid/yy; 
	nregr++;
	double *xxin = xxi + try*(nregr-1);
	get_row(cache, start+best, xxin, 
		covariances, x, nsubject, xord, female, phase);
	double *xxik = xxi;
	for (int k=0; k<(nregr-1); k++, xxik+=try) {
	  int selk = sel[k];
	  double ck = xxin[selk]/xxik[selk];
	  coef[ic++] = ck;
	  for (int j=0; j<try; j++)
	    xxin[j] -= ck*xxik[j];
	}
	int stop = (rsq>=r2stop)||(nregr==pmax)||((r2stop==1)&&(dX2<=2.0)); 
	if (stop) {
	  break;
	}
	else {
	  for (int j=0; j<try; j++) 
	    xy[j] -= bestc*xxin[j];
	  double vn = xxd[best];
	  for (int j=0; j<try; j++) 
	    xxd[j] -= xxin[j]*xxin[j]/vn;
	}
      }
      for (int k=0; k<nregr; k++)
	coef[ic++] = ycoef[k];
      utinv(coef, nregr+1);

 
      SEXP Rule, Rlnames, R2, Pnames, Coefs;
      PROTECT(Rule = allocVector(VECSXP, 3));
      
      PROTECT(Rlnames = allocVector(STRSXP, 3));
      SET_STRING_ELT(Rlnames, 0, mkChar("r.squared"));
      SET_STRING_ELT(Rlnames, 1, mkChar("snps"));
      SET_STRING_ELT(Rlnames, 2, mkChar("coefficients"));
      setAttrib(Rule, R_NamesSymbol, Rlnames);
      
      PROTECT(R2 = allocVector(REALSXP, 1));
      *REAL(R2) = rsq;
      PROTECT(Pnames = allocVector(STRSXP, nregr));
      PROTECT(Coefs = allocVector(REALSXP, nregr+1));
      double intcpt = snpmean(yi, female, nsubject);
      for (int j=0, ic=(nregr*(nregr-1))/2; j<nregr; j++, ic++) {
	int xsnp = xord[start+sel[j]]-1;
	SET_STRING_ELT(Pnames, j, STRING_ELT(Xsnpnames, xsnp));
	double beta =  (-coef[ic]);
	REAL(Coefs)[j+1] = beta;
	intcpt -= beta*snpmean(x+nsubject*xsnp, female, nsubject);
      }
      *REAL(Coefs) = intcpt; 
      SET_VECTOR_ELT(Rule, 0, R2);
      SET_VECTOR_ELT(Rule, 1, Pnames);
      SET_VECTOR_ELT(Rule, 2, Coefs);
      SET_VECTOR_ELT(Result, yord[i]-1, Rule);
      UNPROTECT(5);
    }
    else {
      SET_VECTOR_ELT(Result, yord[i]-1, R_NilValue);
    }
  }
  SEXP IrClass;
  PROTECT(IrClass = allocVector(STRSXP, 1));
  SET_STRING_ELT(IrClass, 0, mkChar("snp.reg.imputation"));
  setAttrib(Result, R_ClassSymbol, IrClass);
  SET_S4_OBJECT(Result);

  /* Tidy up */

  Free(xy);
  Free(xxi);
  Free(xxd);
  Free(coef);
  Free(ycoef);
  free_window(cache);
  UNPROTECT(2);
  return Result;
}

double covariances(int i, int j, va_list ap) {
  unsigned char *snps = va_arg(ap, unsigned char *);
  int N = va_arg(ap, int);
  int *cols = va_arg(ap, int *);
  int *female = va_arg(ap, int *);
  /* int phase = va_arg(ap, int); */
  int ik = N*(cols[i]-1), jk = N*(cols[j]-1);
  return snpcov(snps+ik, snps+jk, female, N);
}

double snpcov(const unsigned char *x, const unsigned char *y, 
	      const int *female, const int N) {
  int sum=0, sumi=0, sumj=0, sumij=0;
  if (female) {
    for (int k=0; k<N; k++) {
      int wt = female[k]? 2: 1;
      int si = (int) *(x++);
      int sj = (int) *(y++);
     if (si && sj) { 
	sum += wt;
	sumi += wt*si;
	sumj += wt*sj;
	sumij += wt*si*sj;
      }
    }
  }
  else {
    for (int k=0; k<N; k++) {
      int si = (int) *(x++);
      int sj = (int) *(y++);
      if (si && sj) {
	sum ++;
	sumi += si;
	sumj += sj;
	sumij += si*sj;
      }
    }
  }
  if (sum>1) 
    return ((double) sumij - ((double) sumi* (double) sumj/(double) sum))/
      (double) (sum - 1);
  else
    return 0.0;
}

double snpmean(const unsigned char *x, const int *female, const int N) {
  int sum=0, sumx=0;
  if (female) {
    for (int i=0; i<N; i++) {
      int wt = female[i]? 2: 1;
      int w = (int) *(x++);
      if (w) {
	sum += wt;
	sumx += wt*w;
      }
    }
  }
  else {
    for (int i=0; i<N; i++) {
      int w = (int) *(x++);
      if (w) {
	sum++;
	sumx += w;
      }
    }
  }
  if (sum)
    return (double) sumx/(double) sum - 1.0;
  else
    return NA_REAL;
}

/* Inverse of unit triangular matrix -- diagonal not stored */
void utinv(double *mat, const int N){
  if (N<2)
    return;
  for (int j=1, ij=0; j<N; j++) {
    for (int i=0, is=0; i<j; i++, ij++) {
      double w = mat[ij];
      for (int k=(i+1), k1=ij+1, k2=is; k<j; k++){ 
	w += mat[k1]*mat[k2];
	k1++;
	k2+= (k+1);
      }
      mat[ij] = (-w);
      is += (i+2);
    }
  }
}
  

int bin_search(const double *sorted, const int len, const double value) {
  int low = 0, high = len - 1;
  int mid = (low + high)/2;
  while (low<mid) {
    if (sorted[mid] > value)
      high = mid;
    else if (sorted[mid] < value)
      low = mid;
    else
      return mid;
    mid = (low + high)/2;
  }
  return low;
}

int nearest_N(const double *sorted, const int len, const double value, 
	      const int N) {
  int last = len - N;
  int res = bin_search(sorted, len, value) - N/2;
  if (res<0)
    res = 0;
  if (res>last)
    res = last;
  if ((value - sorted[res])>(sorted[res+N-1] - value)) {
    while (res<last) {
      res++;
      if ((value - sorted[res])<=(sorted[res+N-1] - value))
	return res;
    }
  }
  else {
    while (res>0) {
      res--;
      if ((value - sorted[res])>=(sorted[res+N-1] - value))
	return res;
    }
  }
  return res;
}

/* Create a hash index from a set of names */

index_db create_name_index(const SEXP names) {
  if (TYPEOF(names)!=STRSXP) 
    error("Names not character variable");
  int N = LENGTH(names);
  index_db res = index_create(N);
  for (int i=0; i<N; i++) {
    if (index_insert(res, CHAR(STRING_ELT(names, i)), i)!=0)
      error("Duplicate names");
  }
  return res;
}



/* Do an imputation (on selected rows) - no class/type checking yet */

void do_impute(const unsigned char *snps, const int nrow, 
	       const int *rows, int nuse, 
	       index_db snp_names,
	       SEXP Rule, 
	       double *value_a, double *value_d) {
  SEXP Snps = VECTOR_ELT(Rule, 1);
  int nsnp = LENGTH(Snps);
  SEXP Coefs = VECTOR_ELT(Rule, 2);
  double *coefs = REAL(Coefs);
  double alpha = *coefs;
  if (!rows)
    nuse = nrow;

  for (int j=0; j<nsnp; j++) {
    int jj = index_lookup(snp_names, CHAR(STRING_ELT(Snps, j)));
    if (jj<0)
      error("Couldn't match snp name: %s", CHAR(STRING_ELT(Snps, j)));
    double beta = coefs[j+1];
    for (int r=0, ist=nrow*jj; r<nuse; r++) {
      int i = rows? rows[r]-1: r;
      unsigned char sij = snps[ist+i];
      double var = j? value_a[r]: alpha;
      if (sij && !ISNA(var))
	value_a[r] = var + beta*((double)(sij - 1));
      else 
	value_a[r] = NA_REAL;
    }
  }
  /* I think one may be able to do better than this */
  if (value_d) {
    for (int r=0; r<nuse; r++) {
      double w = value_a[r];
      value_d[r] = w*w/4.0;
    }
  }
}
 

SEXP impute_snps(const SEXP Rules, const SEXP Snps, const SEXP Subset) { 
  SEXP names = getAttrib(Snps, R_DimNamesSymbol);
  index_db name_index = create_name_index(VECTOR_ELT(names, 1));
  int N = nrows(Snps);
  const unsigned char *snps = RAW(Snps);
  int M = LENGTH(Rules);
  int *subset = NULL;
  int nsubj = N;
  SEXPTYPE sutype = TYPEOF(Subset);
  if (sutype==INTSXP) {
    if (LENGTH(Subset)>N)
      error("Dimension error - Subset");
    subset = INTEGER(Subset);
    nsubj = LENGTH(Subset);
  }
  else if (sutype!=NILSXP)
    error("Argument error - Subset");

  SEXP Result, Dimnames;
  PROTECT(Result = allocMatrix(REALSXP, nsubj, M));
  double *result = REAL(Result);
  PROTECT(Dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(Dimnames, 0, VECTOR_ELT(names, 0));
  SET_VECTOR_ELT(Dimnames, 1, getAttrib(Rules, R_NamesSymbol));
  setAttrib(Result, R_DimNamesSymbol, Dimnames);
  for (int j=0; j<M; j++, result+=nsubj) {
    SEXP Rule = VECTOR_ELT(Rules, j);
    do_impute(snps, N, subset, nsubj, name_index, Rule, result, NULL); 
  }
  index_destroy(name_index);
  UNPROTECT(2);
  return Result;
}

/* Summarize an imputation rule set */

SEXP r2_impute(const SEXP Rules) {
  int M = LENGTH(Rules);
  SEXP Result;
  PROTECT(Result = allocMatrix(REALSXP, M, 2));
  double *result = REAL(Result);
  for (int i=0; i<M; i++) {
    SEXP Rule = VECTOR_ELT(Rules, i);
    if (TYPEOF(Rule)==NILSXP){
      result[i] = 1.0;
      result[i+M] = 0;
    }
    else {
      result[i] = *REAL(VECTOR_ELT(Rule, 0));
      result[i+M] = LENGTH(VECTOR_ELT(Rule, 1));
    }
  }
  UNPROTECT(1);
  return Result;
}

