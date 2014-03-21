
# Main recursion written in C
KFnllDO2 <- function(Params, do.obs, do.sat, k.gas, z.mix, irr, wtr){
	
	# ===========================
	# = Unpack and set initials =
	# ===========================
	#!Pseudocode #1: Initial guesses for B, C, and Q t
	c1 <- Params[1] #PAR coeff
	c2 <- Params[2] #log(Temp) coeff
	Q <- exp(Params[3]) # Variance of the process error
	H <- exp(Params[4]) # Variance of observation error
	
	# See KalmanDO_smooth.R comments for explanation of beta
	kz <- k.gas/z.mix # K and Zmix are both vector of length nobs
	# beta <- 1-kz # beta is a vector of length nobs (this beta is for difference equation form)
	beta <- exp(-kz) # This beta is for using the differential equation form
	
	# Set first true value equal to first observation
	alpha <- do.obs[1]#Let's give this model some starting values
	
	# Set process covariance, P, equal to Q
	P <- Q #starting value
	
	# Empty vector for nll's
	nlls <- rep(0,length(do.obs))
		
	# ==================
	# = Main Recursion =
	# ==================
	# kalmanLoopC(double *alpha, double *doobs, double *c1, double *c2, double *P, double *Q, double *H,  double *beta, double *irr, double *wtr, double *kz, double *dosat, int *nobs)
	nlls <- kalmanLoopR(nlls=nlls, alpha=alpha, doobs=do.obs, c1=c1, c2=c2, P=P, Q=Q, H=H, beta=beta, irr=irr, wtr=wtr, kz=kz, dosat=do.sat)

	return(sum(nlls)) # return the sum of nll's
}#End function
