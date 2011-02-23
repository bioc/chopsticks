#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "mla.h"
#include "glm_test.h"


/* Fit GLM, possibly with a stratification in the RHS

Input:

family       GLM family (see below)
link         Link function (see below)
N            # units
M            # X variables
S            # strata (0 means no intercept)
y            y-variable (N-vector)
prior        prior weights (if present)
X            If M>0, N*M matrix of X variables
stratum      If S>1, stratum assignments coded 1...S (N-vector)
maxit        Maximum number of iterations of IRLS algorithm
conv         Proportional change in weighted sum of squares residuals to
             declare convergence
init         If true (non-zero), the iteration starts from initial estimates 
             of fitted values (see below). This option has no effect if
	     no iteration is required

Output:

rank         rank of X after regression on strata
Xb           orthogonal basis for X space (N*rank matrix)
fitted       fitted values 
resid        working residuals (on linear predictor scale) (N-vector)
weights      weights (N-vector)
scale        scale factor (scalar)
df_resid     residual degrees of freedom

Return

0            convergence
1            no convergence after maxit iterations

*/

int glm_fit(int family, int link, int N, int M, int S,
	    const double *y, const double *prior, const double *X, 
	    const int *stratum, int maxit, double conv, int init, 
	    int *rank, double *Xb, 
	    double *fitted, double *resid, double *weights, 
	    double *scale, int *df_resid) {
  const double eta = 1.e-8;       /* Singularity threshold */
  int i = 0, j=0;
  int Nu, dfr, irls;
  int empty = 0;

  /* Is iteration necessary? */

  irls = (M>0) && !((family==GAUSSIAN) && (link==IDENTITY));

  if (!init || !irls) {
    
    /* Fit intercept and/or strata part of model */

    empty = wcenter(y, N, prior, stratum, S, 0, fitted);
  }

  Nu = 0;
  int invalid = 0;
  for (i=0; i<N; i++) {
    double mu = fitted[i];
    double ri, wi;
    double pi = prior? prior[i] : 1.0;
    if (!muvalid(family, mu)) {
      invalid = 1;
      pi = 0.0;
    }
    if (!(pi)) 
      wi = ri = 0.0;
    else {
      Nu ++;
      double Vmu = varfun(family, mu);
      if (link == family) {
	ri = (y[i] - mu)/Vmu;
	wi = pi*Vmu;
      }
      else {
	double D = dlink(link, mu);
	ri = D*(y[i] - mu);
	wi = pi/(D*D*Vmu);
      }
    }
    weights[i] = wi;
    resid[i] = ri;
  }

  /* If M>0, include covariates */

  int x_rank = 0, convg = 0, iter = 0;
  if (M) {
    convg = 0;
    double wss_last = 0.0;
    if (irls) {

      /* IRLS algorithm */

      double *yw = (double *) calloc(N, sizeof(double));
      while(iter<maxit && !convg) {
	for (i=0; i<N; i++) 
	  yw[i] = resid[i] + linkfun(family, fitted[i]);
	empty = wcenter(yw, N, weights, stratum, S, 1, resid);
	const double *xi = X;
	double *xbi = Xb;
	x_rank = 0;
	for (i=0; i<M; i++, xi+=N) {
	  double ssx = wssq(xi, N, weights);
	  wcenter(xi, N, weights, stratum, S, 1, xbi);
	  double *xbj = Xb;
	  for (j=0; j<x_rank; j++, xbj+=N)
	    wresid(xbi, N, weights, xbj, xbi);
	  double ssr = wssq(xbi, N, weights);
	  if (ssr/ssx > eta) {
	    wresid(resid, N, weights, xbi, resid);
	    x_rank++;
	    xbi+=N;
	  }
	}
	double wss=0.0;
	Nu = 0;
	for (i=0; i<N; i++) {
	  double D, Vmu, ri, wi;
	  double mu = invlink(family, yw[i] - resid[i]);
	  fitted[i] = mu;
	  double pi = prior? prior[i] : 1.0;
	  if (!(pi && muvalid(family, mu) && (weights[i]>0.0))) 
	    wi = ri = 0.0;
	  else {
	    Vmu = varfun(family, mu);
	    Nu ++;
	    if (link == family) {
	      ri = (y[i] - mu)/Vmu;
	      wi = pi*Vmu;
	    }
	    else {
	      D = dlink(link, mu);
	      ri = D*(y[i] - mu);
	      wi = pi+D*D*Vmu;
	    }
	    wss += wi*ri*ri;
	  }
	  weights[i] = wi;
	  resid[i] = ri;
	}
	convg = (family==2) || (Nu<=0) ||
	  (iter && (abs(wss-wss_last)/wss_last < conv));
	wss_last = wss;
	iter ++;
      }
      for (i=0; i<N; i++)
	fitted[i] =  invlink(family, yw[i] - resid[i]);
      free(yw);
    }
    else {  

      /* Simple linear Gaussian case */

      const double *xi = X;
      double *xbi = Xb;
      x_rank = 0;
      for (i=0; i<M; i++, xi+=N) {
	double ssx = wssq(xi, N, weights);
	wcenter(xi, N, weights, stratum, S, 1, xbi);
	double *xbj = Xb;
	for (j=0; j<x_rank; j++, xbj+=N)
	  wresid(xbi, N, weights, xbj, xbi);
	double ssr = wssq(xbi, N, weights);
	if (ssr/ssx > eta) {
	  wresid(resid, N, weights, xbi, resid);
	  x_rank++;
	  xbi+=N;
	}
      }
      wss_last = wssq(resid, N, weights);
      for (i=0; i<N; i++)
	fitted[i] = y[i] - resid[i];
    }
    dfr = Nu  - S + empty - x_rank;
    if (family>2) 
      *scale = wss_last/(dfr);
    else
      *scale = 1.0;
  }
  else {
    if ((S>1) && invalid)  /* Need to recalculate empty stratum count  */
      empty = wcenter(fitted, N, weights, stratum, S, 0, fitted); 
    dfr = Nu - S + empty;
    if (family>2) 
      *scale = wssq(resid, N, weights)/(dfr);
    else
      *scale = 1.0;    
    x_rank = 0;
  }
  *df_resid = dfr>0? dfr : 0;
  *rank = x_rank;
  return(irls && !convg);
}


/* 

Variance function

family:
  1    Binomial
  2    Poisson
  3    Gaussian
  4    Gamma

*/
      
double varfun(int family, double mu){
  switch (family) {
  case 1: return((mu*(1.0-mu)));  /* Binomial */
  case 2: return(mu);             /* Poisson */
  case 3: return(1.0);            /* Gaussian */
  case 4: return(mu*mu);          /* Gamma */
  default: return(0.0);
  }
}

/* Valid values for fitted value, mu. 

If, during iteration, an invalid value is returned, the case is omitted 

*/

int muvalid(int family, double mu) {
  const double minb = 0.0001, maxb = 0.9999, minp = 0.0001;
  switch (family) {
  case 1: return(mu>minb && mu<maxb);    /* Binomial */
  case 2: return(mu>minp);               /* Poisson */
  default: return(1);                     
  }
}

/* Link function

Link
  1    Logit
  2    Log
  3    Identity
  4    Inverse

Note that a canonical link shares the code of the corresponding family
so that the test for canonical link is (link==family)

*/


double linkfun(int link, double mu) {
  switch (link) {
  case 1: return(log(mu/(1.0-mu)));    /* Logit */
  case 2: return(log(mu));             /* Log */
  case 3: return(mu);                  /* Identity */
  case 4: return(1.0/mu);              /* Inverse */
  default: return 0.0;
  }
}

double invlink(int link, double eta) {
  switch (link) {
  case 1: return(exp(eta)/(1.0+exp(eta))); /* Logit */
  case 2: return(exp(eta));                /* Log */
  case 3: return(eta);                     /* Identity */
  case 4: return(1.0/eta);                 /* Inverse */
  default: return(0.0);
  }
}

double dlink(int link, double mu) {
  switch (link) {
  case 1: return(1.0/(mu/(1.0-mu)));     /* Logit */
  case 2: return(1.0/mu);                /* Log */
  case 3: return(1.0);                   /* Identity */
  case 4: return(1.0/(mu*mu));           /* Inverse */
  default: return 0.0;
  }
}


/* 

GLM score test 

Input:

P         Number of new explanatory variables to be added 
Z         N*P matrix containing these (destroyed)
C         If robust variance estimate to be used, number of clusters
          (if C==1, each unit forms a cluster)
cluster   If C>1, cluster assignments code 1...C (N-vector)
max_R2    For P>1, the maximum value of R^2 between each column and revious 
          columns (after regression on X and strata)

For all other input arguments, see glm_fit, but note that M now coincides 
with rank -- the number of columns in Xb

Output:

score       P*1 vector of scores
score-var   Upper triangle of variance-covariance matrix of scores

*/

void glm_score_test(int N, int M, int S, const int *stratum, 
		    int P, const double *Z, int C, const int *cluster,
		    const double *resid, const double *weights, 
		    const double *Xb, double scale,
		    double max_R2, 
		    double *score, double *score_var) {
  const double eta = 1.e-8;   /* First stage singularity test */
  const double *Zi = Z;
  double *Zr, *Zri;
  double *U = NULL, *Ui = NULL;

  /* Work array */

  Zri = Zr = (double *) calloc(N*P, sizeof(double));
  int nc = 0;
  if (C) {
    nc = (C==1)? N: C;
    Ui = U = (double *) calloc(nc*P, sizeof(double));
  }
    

  /* Main algorithm */
 
  for (int i=0, ij=0; i<P; i++, Zi+=N, Zri+=N) {

    /* Regress each column of Z on strata indicators and X basis */

    double ssz = wssq(Zi, N, weights);
    wcenter(Zi, N, weights, stratum, S, 1, Zri);
    const double *Xbj = Xb;
    for (int j=0; j<M; j++, Xbj+=N) 
      wresid(Zri, N, weights, Xbj, Zri);
    double ssr = wssq(Zri, N, weights);

    if (ssr/ssz > eta) {     /* Singularity test */
 
      if (C) {

        /* Add new column to matrix of score contributions */

	if (C==1) {
	  for (int k=0; k<N; k++)
	    Ui[k] = Zri[k]*weights[k]*resid[k];
	}
	else {
	  memset(Ui, 0x00, nc*sizeof(double));
	  for (int k=0; k<N; k++) {
	    int kc = cluster[k] - 1;
	    Ui[kc] += Zri[k]*weights[k]*resid[k];
	  }
	}

	/* Update score and score-variance */

	score[i] = wsum(Ui, nc, NULL);
	double *Uj = U;
	for (int j=0; j<i; j++, Uj+=nc) 
	  score_var[ij++] = wspr(Ui, Uj, nc, NULL);
	score_var[ij++] = wssq(Ui, nc, NULL);
      }
      else {

	/* Update score and score-variance */

	score[i] = wspr(Zri, resid, N, weights);
	double *Zrj = Zr;
	for (int j=0; j<i; j++, Zrj+=N) 
	  score_var[ij++] = scale*wspr(Zri, Zrj, N, weights);
	score_var[ij++] = scale*wssq(Zri, N, weights);
      }
    }
    else {
      memset(Zri, 0x00, N*sizeof(double));
      score[i] = 0.0;
      for (int j=0; j<=i; j++) 
	score_var[ij++] = 0.0;
    }
  }
  free(Zr);
  if (C)
    free(U);
}
