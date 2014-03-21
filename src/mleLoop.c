# include <R.h>
# include <Rmath.h>

void mleLoopC(double *alpha, double *doobs, double *c1, double *c2, double *beta, double *irr, double *wtr, double *kz, double *dosat, int *nobs){
	int i;
	double a1;
	// double alpha=*alpha, doobs=*doobs, c1=*c1, c2=*c2, beta=*beta, irr=*irr, wtr=*wtr, kz=*kz, dosat=*dosat;
	int ni=*nobs;
	for(i=1; i<ni; i++){
		a1 = *c1 * irr[i - 1] + *c2 * log(wtr[i - 1]) + kz[i - 1] * dosat[i - 1];
		alpha[i] = a1 / kz[i - 1] + -exp(-1 * kz[i - 1]) * a1 / kz[i - 1] + beta[i - 1] * alpha[i-1];
	}
}
