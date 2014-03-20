

metab.mle <- function(do.obs, do.sat, k.gas, z.mix, irr, wtr, ...){
	Q0 <- ((diff(range(do.obs,na.rm=TRUE))-mean(do.obs,na.rm=TRUE))^2/n.obs)
	guesses <- c(1E-4,1E-4,log(Q0))
	fit <- optim(guesses, fn=mleNLL, do.obs=do.obs, do.sat=do.sat, k.gas=k.gas, z.mix=z.mix, irr=irr, wtr=wtr, ...)
	pars0 <- fit$par
	
	pars <- c("gppCoeff"=pars0[1], "rCoeff"=pars0[2], "Q"=exp(pars0[3]))
	
	# ====================================
	# = Use fits to calculate metabolism =
	# ====================================
	GPP <- sum(pars[1]*irr, na.rm=TRUE)
	R <- sum(pars[2]*log(wtr), na.rm=TRUE)
	
	return(list("params"=pars, "metab"=c("GPP"=GPP,"R"=R)))
}


mleNLL <- function(Params, do.obs, do.sat, k.gas, z.mix, irr, wtr){
	c1 <- Params[1] #PAR coeff
	c2 <- Params[2] #log(Temp) coeff
	Q <- exp(Params[3]) # Variance of the process error
	
	# See KalmanDO_smooth.R comments for explanation of beta
	kz <- k.gas/z.mix # K and Zmix are both vector of length nobs
	beta <- exp(-kz) # This beta is for using the differential equation form
	
	# Set first true value equal to first observation
	alpha <- rep(NA, length(do.obs))
	alpha[1] <- do.obs[1]#Let's give this model some starting values

	for(i in 2:length(do.obs)){
		a1 <- c1*irr[i-1] + c2*log(wtr[i-1]) + kz[i-1]*do.sat[i-1]
		alpha[i] <- a1/kz[i-1] + -exp(-kz[i-1])*a1/kz[i-1] + beta[i-1]*alpha[i-1] # NOTE: beta==exp(-kz); kz=K/Zmix
	}
	return(-sum(dnorm(do.obs, alpha, sd=sqrt(Q), log=TRUE), na.rm=TRUE))
}#End function
