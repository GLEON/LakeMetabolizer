metab.ols <- function(do.obs, do.sat, k.gas, z.mix, irr, wtr, ...){

	nobs <- length(do.obs)
	
	mo.args <- list(...)
	if("datetime"%in%names(mo.args)){ # check to see if datetime is in the ... args
		datetime <- mo.args$datetime # extract datetime
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

	do.diff <- diff(do.obs)


	inst_flux <- (k.gas/freq) * (do.sat - do.obs)  # positive is into the lake
	
	# flux <- (inst_flux[1:(n.obs-1)] + inst_flux[-1])/2
	flux <- inst_flux[-1]

	
	noflux.do.diff <- do.diff - flux/z.mix[-1]
	
	lntemp <- log(wtr)
	mod <- lm(noflux.do.diff ~ irr[-1] + lntemp[-1] -1) # note that if we decied to use resp=rho*log(Temp), you would do lm(do~irr+lntemp-1) (mind the -1)

	rho <- mod[[1]][2] 
	iota <- mod[[1]][1]
	mod.matrix <- model.matrix(mod)
	gpp <- mean(iota*mod.matrix[,1], na.rm=TRUE) * freq
	resp <- mean(rho*mod.matrix[,2], na.rm=TRUE) * freq
	nep <- gpp + resp
	#other ways to get nep:
	# nep = rho + gpp
	# nep = sum(fitted(mod), na.rm=TRUE) # even if there are NA's in the response variable, they shouldn't be included in fitted() ....
	# sum(noflux.do.diff) # i think this should be the same as sum(fitted(mod)) b/c model residuals sum to 0 ... right?
	# also note that NEP is gpp+rho (rho is negative by this convention, which is consistent w/ Kalman, Bayes, mle – unsure of BK)

	results <- list("mod"=mod, "metab"=data.frame("GPP"=gpp, "R"=resp, "NEP"=nep))
	# attr(results, "lm.mod") <- mod
	return(results)
	
}
