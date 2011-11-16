int wcenter(const double *y, int n, const double *weight, const int *stratum, 
	    int nstrata, int resid, double *ynew);

int wresid(const double *y, int n, const double *weight, const double *x, 
	   double *ynew);

double wssq(const double *y, int n, const double *weight);
