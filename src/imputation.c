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

#define min2(x, y) (x<y? x: y)

int bin_search(const double *sorted, const int len, const double value);
int nearest_N(const double *sorted, const int len, const double value, 
	      const int N);
double covariances(int i, int j, va_list ap);
double snpcov(const unsigned char *x, const unsigned char *y, 
	      const int *female, const int N, const int phase, 
	      const double minA);
double snpmean(const unsigned char *x, 
	       const int *female, const int N);
void utinv(double *, const int);


SEXP snp_impute(const SEXP X, const SEXP Y, const SEXP Xord, const SEXP Yord,
		const SEXP Xpos, const SEXP Ypos, const SEXP Phase, 
		const SEXP Try, const SEXP Stop,
		const SEXP Hapcontr, const SEXP EMcontr,
		const SEXP MinA){

  int try =  *INTEGER(Try);   /* Number to search */
  if (LENGTH(Stop)!=3)
    error("Stop argument not of length 3");
  double r2stop = REAL(Stop)[0]; /* R^2 value to stop inclusion */
  if (r2stop>1.0)
    r2stop = 1.0;
  int pmax = (int) REAL(Stop)[1];   /* Maximum number of predictor variables */
  if (pmax>try)
    pmax = try;
  double dR2 = REAL(Stop)[2]; /* Minimum increase in R^2 to include */
  if (dR2<0)
    dR2 = NA_REAL;
  int phase = *LOGICAL(Phase);   /* haploid or diploid computation */
  double minA = *REAL(MinA);     /* min paired data test */

  /* Fully phased haplotype-based  imputation control args */

  if (LENGTH(Hapcontr)!=2)
    error("Hapcntr argument not of length 2");
  double hapr2 = REAL(Hapcontr)[0]; /* R^2 value to force try */
  double hapimp = REAL(Hapcontr)[1]; /* Required gain in 1-R^2 to persist */
  if (LENGTH(EMcontr)!=2)
    error("EMcontr argument not of length 2");
  int maxit = REAL(EMcontr)[0];  /* Max EM iterations */
  double emtol = REAL(EMcontr)[1]; /* EM convergence tolerance */

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

  double *xy = (double *)Calloc(try, double);       /* XY covariances */
  double *xxd = (double *)Calloc(try, double);      /* Diagonals of XX */
  double *xxi = (double *)Calloc(try*pmax, double); /* Row of XX covars */
  int *sel = (int *)Calloc(pmax, int);              /* Selected SNPs */
  double *coef =(double *)Calloc((pmax*(pmax+1))/2, double);/* Coefficients*/ 
  double *ycoef = (double *)Calloc(pmax, double);   /* Y coefficients */
  COV_WIN_PTR cache = new_window(try, 0);            /* Covariance cache */

  /* Work arrays for haplotype phasing etc. */

  int *contin=NULL, *hcontin=NULL;
  int tmax=0, hmax=0;
  double *phap=NULL;
  GTYPE **tables=NULL;
  int *tcell=NULL;
  if (hapr2>0.0) {
    tmax = (1 << 2*(pmax+1));  /* Space for 4x4x..x4 table */
    hmax = (1 << (pmax+1));    /* Space for 2x2x..x2 table */
    tcell = (int *)Calloc(nsubject, int); /* addresses in table */
    contin = (int *)Calloc(tmax, int); 
    if (female) 
      hcontin = (int *)Calloc(tmax, int);
    phap = (double *)Calloc(hmax, double);
    /* gtype->htype lookup tables */
    tables = (GTYPE **)Calloc(pmax+1, GTYPE *);
    for (int i=0; i<=pmax; i++)
      tables[i] = create_gtype_table(i+1);
  }

  /* Result */
 
  SEXP Result;
  PROTECT(Result = allocVector(VECSXP, ny));
  setAttrib(Result, R_NamesSymbol, 
	    VECTOR_ELT(getAttrib(Y, R_DimNamesSymbol),1));

  /* Main loop */

  int n_em_fail=0, maxpred=0;
  for (int i=0; i<ny; i++) {
    unsigned char *yi = y + nsubject*(yord[i]-1);
    /* Minor allele frequency */
    int ng=0, na=0;
    for (int j=0; j<nsubject; j++) {
      int yij = (int) yi[j];
      if (yij) {
	ng++;
	na += yij;
      }
    }
    if (ng>0) {
      double maf = (double) (na - ng)/ (double) (2*ng);
      if (maf>0.5)
	maf = 1.0 - maf;
      double yy = snpcov(yi, yi, female, nsubject, phase, minA);
      if (!ISNA(yy)) {
	int start = nearest_N(xpos, nx, ypos[i], try);
	for (int j=0; j<try; j++) { 
	  int jx = nsubject*(xord[start+j]-1);
	  xy[j] = snpcov(x+jx, yi, female, nsubject, phase, minA);
	}
	move_window(cache, start);
	get_diag(cache, xxd, covariances, x, nsubject, xord, female, phase, 
		 minA);
	double resid = yy;
	double rsq = 0.0;
	int nregr=0, ic = 0;
	while (1) {
	  /* Find next snp to be included */
	  double max_due = 0.0;
	  int best = -1;
	  for (int j=0; j<try; j++) {
	    double xxj = xxd[j];
	    double xyj = xy[j];
	    if (xxj==0.0 || ISNA(xxj) || ISNA(xyj))
	      continue;
	    double xyj2 = xyj*xyj;
	    if (xyj2>(xxj*resid)) { /* r^2 > 1 */
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
	  if (best<0)
	    break;
	  double *xxin = xxi + try*nregr;
	  get_row(cache, start+best, xxin, 
		  covariances, x, nsubject, xord, female, phase, minA);
	  double *xxik = xxi;
	  int reject = 0;
	  for (int k=0; k<nregr; k++, xxik+=try) {
	    int selk = sel[k];
	    double xxink = xxin[selk], xxikk = xxik[selk];
	    if (ISNA(xxink) || ISNA(xxikk)) {
	      reject = 1;
	      ic -= k;
	      break;
	    }
	    double ck = xxink/xxikk;
	    coef[ic++] = ck;
	    for (int j=0; j<try; j++) {
	      double w1 = xxin[j], w2 = xxik[j];
	      xxin[j] = ISNA(w1) || ISNA(w2)? NA_REAL: w1 - ck*w2;
	    }
	  }
	  if (reject) {
	    xxd[best] = 0.0;
	    continue;
	  }
	  sel[nregr] = best; /* Save index */
	  double bestc = xy[best]/xxd[best]; 
	  ycoef[nregr] = bestc; /* Save regression coefficient */
	  if (ISNA(dR2))
	    dR2 = 2/ (double) (nsubject-nregr-2);
	  double deltaR2 =  max_due/resid;
	  resid -= max_due;
	  rsq = 1.0 - resid/yy; 
	  nregr++;
	  int stop = (rsq>=r2stop)||(nregr==pmax)||(deltaR2<dR2); 
	  if (stop) 
	    break;

	  double vn = xxd[best];
	  xxd[best] = 0.0;
	  for (int j=0; j<try; j++) {
	    double w = xxin[j];
	    if (!ISNA(w)) {
	      xy[j] -=  bestc*w;
	      w = xxd[j]-w*w/vn;
	      xxd[j] = w>0.0? w: 0.0;
	    }
	  }
	}

	/* Use regression imputation or phased haplotypes? */

	double r2hap = 0.0, r2gain = -1.0;
	if (nregr>1 && rsq<hapr2) {
	  /* Calculate phased haplotypes */  
	  for (int j=0; j<nsubject; j++) {
	    tcell[j] = (int) yi[j];
	  }
	  for (int k=0, sh=2; k<nregr; k++, sh+=2) {
	    unsigned char *xk = x + nsubject*(xord[start+sel[k]]-1);
	    for (int j=0; j<nsubject; j++) {
	      int  xkj = (int) xk[j];
	      tcell[j] = tcell[j] | (xkj << sh);
	    }
	  }
	  /* Calculate contingency table */
	  int dim = nregr+1;
	  memset(contin, 0x00, tmax*sizeof(int));
	  if (hcontin)
	    memset(hcontin, 0x00, tmax*sizeof(int));
	  for (int j=0; j<nsubject; j++) {
	    if (female && !female[j])
	      hcontin[tcell[j]]++;
	    else
	      contin[tcell[j]]++;
	  }
	  /* EM algorithm */
	  int em_fail = emhap(dim, contin, hcontin, tables[nregr], 
			      maxit, emtol, phap);
	  if (em_fail>=0) {
	    n_em_fail += em_fail; 
	    r2hap = gen_r2(nregr, phap, tables[nregr-1]);
	    r2gain = (r2hap - rsq)/(1.0 - rsq);
	  }
	}

	if (nregr>0) {

	  /* save imputation rule */
	
	  if (nregr>maxpred)
	    maxpred = nregr;

	  SEXP Rule, Rlnames, Maf, R2, Pnames, Coefs;
	  PROTECT(Rule = allocVector(VECSXP, 4));
	  
	  PROTECT(Rlnames = allocVector(STRSXP, 4));
	  SET_STRING_ELT(Rlnames, 0, mkChar("maf"));
	  SET_STRING_ELT(Rlnames, 1, mkChar("r.squared"));
	  SET_STRING_ELT(Rlnames, 2, mkChar("snps"));
	  
	  PROTECT(Maf = allocVector(REALSXP, 1));
	  *REAL(Maf) = maf;
	  PROTECT(R2 = allocVector(REALSXP, 1));
	  PROTECT(Pnames = allocVector(STRSXP, nregr));
	  for (int j=0; j<nregr; j++) {
	    int xsnp =  xord[start+sel[j]]-1;
	    SET_STRING_ELT(Pnames, j, STRING_ELT(Xsnpnames, xsnp));
	  }

	  if (r2gain<hapimp) {

	    /* Save regression imputation */
	    
	    SET_STRING_ELT(Rlnames, 3, mkChar("coefficients"));
	    *REAL(R2) = rsq>1.0? 1.0: rsq;

	    for (int k=0; k<nregr; k++)
	      coef[ic++] = ycoef[k];
	    utinv(coef, nregr+1);
	    PROTECT(Coefs = allocVector(REALSXP, nregr+1));
	    double *coefs = REAL(Coefs);
	    double intcpt = snpmean(yi, female, nsubject);
	    for (int j=0, ic=(nregr*(nregr-1))/2; j<nregr; j++, ic++) {
	      int xsnp = xord[start+sel[j]]-1;
	      double beta =  (-coef[ic]);
	      coefs[j+1] = beta;
	      intcpt -= beta*snpmean(x+nsubject*xsnp, female, nsubject);
	    }
	    coefs[0] = intcpt;
	  }
	  else {

	    /* save phased  haplotype imputation */
	    
	    SET_STRING_ELT(Rlnames, 3, mkChar("hap.probs"));
	    *REAL(R2) = r2hap;
	    int lenp = 1 << (nregr+1);
	    PROTECT(Coefs = allocVector(REALSXP, lenp));
	    double *coefs = REAL(Coefs);
	    for (int j=0; j<lenp; j++) 
	      coefs[j] = phap[j];
	  }
	  SET_VECTOR_ELT(Rule, 0, Maf);
	  SET_VECTOR_ELT(Rule, 1, R2);
	  SET_VECTOR_ELT(Rule, 2, Pnames);
	  SET_VECTOR_ELT(Rule, 3, Coefs);
	  setAttrib(Rule, R_NamesSymbol, Rlnames);
	  SET_VECTOR_ELT(Result, yord[i]-1, Rule);
	  UNPROTECT(6);
	}
	else {
	  /* No valid predictors */
	  SET_VECTOR_ELT(Result, yord[i]-1, R_NilValue);
	}
      }
      else {
	/*MAF too low  */
	SET_VECTOR_ELT(Result, yord[i]-1, R_NilValue);
      }
    }
    else {
      /* No data */
      SET_VECTOR_ELT(Result, yord[i]-1, R_NilValue);
    }
  }
  if (n_em_fail)
    warning("Maximum iterations for %d haplotype imputation rules", n_em_fail);

  SEXP IrClass, Package;
  PROTECT(IrClass = allocVector(STRSXP, 1));
  SET_STRING_ELT(IrClass, 0, mkChar("snp.reg.imputation"));
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("snpMatrix"));
  setAttrib(IrClass, install("package"), Package);
  classgets(Result, IrClass);
  SEXP Maxpred;
  PROTECT(Maxpred = allocVector(INTSXP, 1));
  INTEGER(Maxpred)[0] = maxpred;
  setAttrib(Result, install("Max.predictors"), Maxpred);
  SET_S4_OBJECT(Result);

  /* Tidy up */

  Free(xy);
  Free(xxi);
  Free(xxd);
  Free(sel);
  Free(coef);
  Free(ycoef);
  free_window(cache);
  if (hapr2>0.0) {
    Free(contin);
    if (hcontin)
      Free(hcontin);
    Free(phap);
    Free(tcell);
    for (int i=0; i<=pmax; i++)
      destroy_gtype_table(tables[i], i+1);
    Free(tables);
  }
  UNPROTECT(4);
  return Result;
}

double covariances(int i, int j, va_list ap) {
  unsigned char *snps = va_arg(ap, unsigned char *);
  int N = va_arg(ap, int);
  int *cols = va_arg(ap, int *);
  int *female = va_arg(ap, int *);
  int phase = va_arg(ap, int); 
  double minA = va_arg(ap, double);
  int ik = N*(cols[i]-1), jk = N*(cols[j]-1);
  return snpcov(snps+ik, snps+jk, female, N, phase, minA);
}


double snpcov(const unsigned char *x, const unsigned char *y, 
	  const int *female, const int N, const int phase, const double minA) {
  int n1=0, n2=0, nt=0, sx=0, sy=0, sxy=0;
  double cov, n11;
  if (phase) {
    if (female)
      error("phase=TRUE not yet implemented for the X chromosome");
    error("phase=TRUE not yet implemented");
    return NA_REAL;
  }
  else {
    if (female) {
      for (int k=0; k<N; k++) {
	int xk = (int) *(x++);
	int yk = (int) *(y++);
	if (xk && yk) { 
	  xk--;
	  yk--;
	  if (female[k]) {
	    n2++;
	  }
	  else {
	    n1++;
	    xk/=2;
	    yk/=2;
	  }
	  sx += xk;
	  sy += yk;
	  sxy += xk*yk;
	}
      }
      nt = 2*n2 + n1;
      if (nt<2)
	return NA_REAL;
      int nt_1 = nt-1;
      double p2 = (double)(2*n2)/(double)nt;
      double ps = (double)sx * (double)sy;
      cov = ((double)sxy - (1.0+p2)*ps/(double)nt)/
	((double)nt_1 - p2);
      n11 = (double)nt_1*(sxy - p2*ps/(double)nt_1)/((double)nt_1-p2);
    }
    else {
      for (int k=0; k<N; k++) {
	int xk = (int) *(x++);
	int yk = (int) *(y++);
	if (xk && yk) {
	  xk--;
	  yk--;
	  n2++;
	  sx += xk;
	  sy += yk;
	  sxy += xk*yk;
	}
      }
      if (n2 < 2)
	return NA_REAL;
      nt = 2*n2;
      double ps = (double)sx * (double)sy;
      double n_1 = n2 - 1;
      cov = 0.5*((double)sxy - ps/(double)n2)/(double)n_1;
      double twon_1 = nt - 1;
      n11 = (double)twon_1*((double)sxy - ps/(double)twon_1)/
	(2.0*(double)n_1);
    }
    double test = (cov > 0.0)?
      min2(n11, nt - sx - sy + n11):
      min2(sx - n11, sy - n11);
    /*    printf("n11 = %lf, test = %lf, ", n11, test); */
    if (test < minA)
      return NA_REAL;
    return cov;
  }
}

/*  Routine for testing snpcov */

SEXP snpcov_test(const SEXP X, const SEXP i, const SEXP j, const SEXP minA) {
  int ii = *INTEGER(i) - 1;
  int jj = *INTEGER(j) - 1;
  int N = nrows(X);
  double ma = *REAL(minA);
  unsigned char *x = RAW(X);
  double mycov = snpcov(x+N*ii, x+N*jj, 0, N, 0, ma);
  Rprintf("N = %d, cov = ", N);
  if (ISNA(mycov))
    Rprintf("NA_REAL\n");
  else
    Rprintf("%lf\n", mycov);
  SEXP Result = allocVector(REALSXP, 1);
  *REAL(Result) = mycov;
  return Result;
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
      if (ISNA(w))
	warning("Bug: NAs in triangular coefficients matrix");
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
	       SEXP Rule, GTYPE **gt2ht, 
	       double *value_a, double *value_d) {
  SEXP Snps = VECTOR_ELT(Rule, 2);
  int nsnp = LENGTH(Snps);
  SEXP Coefs = VECTOR_ELT(Rule, 3);
  int ncoefs = LENGTH(Coefs);
  double *coefs = REAL(Coefs);
  double alpha = *coefs;
  if (!rows)
    nuse = nrow;

  if (ncoefs==(nsnp+1)) { /* Regression imputation */
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
  else { /* Haplotype imputation */
    int *gt = (int *)Calloc(nuse, int);
    memset(gt, 0x00, nuse*sizeof(int));
    /* Calculate predictor genotypes */
    for (int j=0, sh=0; j<nsnp; j++, sh+=2) {
      int jj = index_lookup(snp_names, CHAR(STRING_ELT(Snps, j)));
      if (jj<0)
	error("Couldn't match snp name: %s", CHAR(STRING_ELT(Snps, j)));
      for (int r=0, ist=nrow*jj; r<nuse; r++) {
	int i = rows? rows[r]-1: r;
	int sij = (int)snps[ist+i];
	gt[r] = gt[r] | (sij << sh);
      }
    }
    /* 
       Score genotypes
       Perhaps more efficient to compute a lookup table at outset
    */

    const GTYPE *gtab = gt2ht[nsnp-1];
    for (int i=0; i<nuse; i++) {
      double score[3];
      int gti = gt[i];
      if (gti) {
	predict_gt(nsnp, gti, coefs, gtab, score);
	value_a[i] = score[1]+2.0*score[2];
	if (value_d)
	  value_d[i] = score[2];
      }
      else {
	value_a[i] = NA_REAL;
	if (value_d)
	  value_d[i] = NA_REAL;
      }
    }
    Free(gt);
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
  int pmax = *INTEGER(getAttrib(Rules, install("Max.predictors")));
  GTYPE **gt2ht = (GTYPE **)Calloc(pmax, GTYPE *); 
  for (int i=0; i<pmax; i++)
    gt2ht[i] = create_gtype_table(i+1);
  for (int j=0; j<M; j++, result+=nsubj) {
    SEXP Rule = VECTOR_ELT(Rules, j);
    if (isNull(Rule))
      for (int i=0; i<nsubj; i++)
	result[i] = NA_REAL;
    else
      do_impute(snps, N, subset, nsubj, name_index, Rule, gt2ht, result, NULL); 
  }
  index_destroy(name_index);
  for (int i=0; i<pmax; i++) 
    destroy_gtype_table(gt2ht[i], i+1);
  Free(gt2ht);
  UNPROTECT(2);
  return Result;
}

/* Extract r-squared and size from an imputation rule set */

SEXP r2_impute(const SEXP Rules) {
  int M = LENGTH(Rules);
  SEXP Result;
  PROTECT(Result = allocMatrix(REALSXP, M, 2));
  double *result = REAL(Result);
  for (int i=0; i<M; i++) {
    SEXP Rule = VECTOR_ELT(Rules, i);
    if (TYPEOF(Rule)==NILSXP){
      result[i] = NA_REAL;
      result[i+M] = NA_REAL;
    }
    else {
      result[i] = *REAL(VECTOR_ELT(Rule, 1));
      int nsnp = LENGTH(VECTOR_ELT(Rule, 2));
      int nco = LENGTH(VECTOR_ELT(Rule, 3));
      result[i+M] = nsnp*2 - 1 + (nco>(1+nsnp));
    }
  }
  UNPROTECT(1);
  return Result;
}

