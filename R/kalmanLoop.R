# kalmanLoopC(double *alpha, double *doobs, double *c1, double *c2, double *P, double *Q, double *H,  double *beta, double *irr, double *wtr, double *kz, double *dosat, int *nobs)
kalmanLoopR <- function(nlls, alpha, doobs, c1, c2, P, Q, H, beta, irr, wtr, kz, dosat){
	nobs <- length(doobs)
	a.loop <- .C("kalmanLoopC", nlls=as.double(nlls), as.double(alpha), as.double(doobs), as.double(c1), as.double(c2), as.double(P), as.double(Q), as.double(H), as.double(beta), as.double(irr), as.double(wtr), as.double(kz), as.double(dosat), as.integer(nobs), PACKAGE="LakeMetabolizer")
	return(a.loop[["nlls"]])
}
