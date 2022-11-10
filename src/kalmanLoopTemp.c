# include <R.h>
# include <Rmath.h>

void kalmanLoopTempC(double *nlls, double *alpha, double *c1, double *P, double *Q, double *H,  double *beta, double *watts, double *wtr, int *nobs){
	int i;
	double eta, Eff;
	int ni=*nobs;
	for(i=1; i<ni; i++){

		*alpha = *beta * (*alpha) + *c1 * watts[i-1];

		*P = pow(*beta, 2) * *P + *Q;
		eta = wtr[i] - *alpha;
		Eff = *P + *H;
		*alpha = *alpha + *P/Eff*eta;
		*P = *P - pow(*P,2)/Eff;

		nlls[i] = 0.5*log(2*M_PI) + 0.5*log(Eff) + 0.5*pow(eta,2)/Eff;
	}
}
