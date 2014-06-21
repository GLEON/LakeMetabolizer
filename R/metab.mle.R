
# ====================
# = Main mle wrapper =
# ====================
#'@export
metab.mle <- function(do.obs, do.sat, k.gas, z.mix, irr, wtr, ...){

	nobs <- length(do.obs)
	
	mm.args <- list(...)
	if("datetime"%in%names(mm.args)){ # check to see if datetime is in the ... args
		datetime <- mm.args$datetime # extract datetime
		freq <- calc.freq(datetime) # calculate sampling frequency from datetime
		if(nobs!=freq){ # nobs and freq should agree, if they don't issue a warning
			bad.date <- format.Date(datetime[1], format="%Y-%m-%d")
			warning("number of observations on ", bad.date, " (", nobs, ") ", "does not equal estimated sampling frequency", " (", freq, ")", sep="")
		}
	}else{ # if datetime is *not* in the ... args
		warning("datetime not found, inferring sampling frequency from # of observations") # issue a warning (note checks in addNAs)
		# NOTE: because of the checks in addNA's, it is unlikely a user would receive this warning via metab()
		# warning will only be seen through direct use of metab.bookkeep when datettime is not supplied
		freq <- nobs
	}
	
	chk.list <- list(do.obs, irr, do.sat, z.mix, k.gas, wtr)
	if(!all(sapply(chk.list, is.numeric)) || !all(sapply(chk.list, is.vector))){
		stop('All metab.mle inputs must be numeric vectors.')
	}
 
	if(!all(nobs==sapply(chk.list, length))){
		stop('All input data to metab.mle must be the same length')
	}
	
	Q0 <- ((diff(range(do.obs,na.rm=TRUE)) - mean(do.obs,na.rm=TRUE))^2 / length(do.obs))
	guesses <- c(1E-4, 1E-4, log(Q0))
	
	fit <- optim(guesses, fn=mleNLL, do.obs=do.obs, do.sat=do.sat, k.gas=(k.gas/freq), z.mix=z.mix, irr=irr, wtr=wtr)
	pars0 <- fit$par
	
	pars <- c("gppCoeff"=pars0[1], "rCoeff"=pars0[2], "Q"=exp(pars0[3]))
	
	# ====================================
	# = Use fits to calculate metabolism =
	# ====================================
	GPP <- mean(pars[1]*irr, na.rm=TRUE) * freq
	R <- mean(pars[2]*log(wtr), na.rm=TRUE) * freq
	
	return(list("params"=pars, "metab"=c("GPP"=GPP,"R"=R,"NEP"=GPP+R)))
}

# ==========================
# = The R loop for mle NLL =
# ==========================
mleLoopR <- function(alpha, doobs, c1, c2, beta, irr, wtr, kz, dosat){
	nobs <- length(doobs)
	a.loop <- .C("mleLoopC", alpha=as.double(alpha), as.double(doobs), as.double(c1), as.double(c2), as.double(beta), as.double(irr), as.double(wtr), as.double(kz), as.double(dosat), as.integer(nobs), PACKAGE="LakeMetabolizer")
	return(a.loop[["alpha"]])
}

# ====================
# = mle NLL function =
# ====================
mleNLL <- function(Params, do.obs, do.sat, k.gas, z.mix, irr, wtr){
	c1 <- Params[1] #PAR coeff
	c2 <- Params[2] #log(Temp) coeff
	Q <- exp(Params[3]) # Variance of the process error
	
	# See KalmanDO_smooth.R comments for explanation of beta
	kz <- k.gas/z.mix # K and Zmix are both vector of length nobs
	beta <- exp(-kz) # This beta is for using the differential equation form
	
	# Set first true value equal to first observation
	alpha <- rep(0, length(do.obs))
	alpha[1] <- do.obs[1]#Let's give this model some starting values

	#R version of C loop
	#for(i in 2:length(do.obs)){
	#	a1 <- c1*irr[i-1] + c2*log(wtr[i-1]) + kz[i-1]*do.sat[i-1]
	#	alpha[i] <- a1/kz[i-1] + -exp(-kz[i-1])*a1/kz[i-1] + beta[i-1]*alpha[i-1] # NOTE: beta==exp(-kz); kz=K/Zmix
	#}
	
	alpha <- mleLoopR(alpha=alpha, doobs=do.obs, c1=c1, c2=c2, beta=beta, irr=irr, wtr=wtr, kz=kz, dosat=do.sat)
	return(-sum(dnorm(do.obs, alpha, sd=sqrt(Q), log=TRUE), na.rm=TRUE))
}#End function
