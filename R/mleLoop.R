# double *alpha, double *doobs, double *c1, double *c2, double *beta, double *irr, double *wtr, double *kz, double *dosat, int *nobs
mleLoopR <- function(alpha, doobs, c1, c2, beta, irr, wtr, kz, dosat){
	# dyn.load("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/src/LakeMetabolizer.so")
	# dyn.unload("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/src/LakeMetabolizer.so")
	nobs <- length(doobs)
	a.loop <- .C("mleLoopC", alpha=as.double(alpha), as.double(doobs), as.double(c1), as.double(c2), as.double(beta), as.double(irr), as.double(wtr), as.double(kz), as.double(dosat), as.integer(nobs), PACKAGE="LakeMetabolizer")
	# a.loop[["alpha"]]
	return(a.loop[["alpha"]])
}
