metab.ols <- function(do.obs, do.sat, irr, k.gas, z.mix){

	n.obs <- length(do.obs)

	do.diff <- diff(do.obs)

	#basically the average of the flux at time T_0 and T_1
	#flux <- -0.5 * k * (D1 + D2 - 2*S)
	#
	# can be re-arranged
	# flux <- 0.5 * [k * (D1-S) + k * (d2-s)]
	# Average of each inst flux
	inst_flux <- k.gas * (do.sat - do.obs)  # positive is into the lake
	
	flux <- (inst_flux[1:(n.obs-1)] + inst_flux[-1])/2

	
	noflux.do.diff <- do.diff - flux/z.mix

	mod <- lm(noflux.do.diff ~ irr) # note that if we decied to use resp=rho*log(Temp), you would do lm(do~irr+lntemp-1) (mind the -1)

	# also note that this model has different structure than Bayes, mle, Kalman 
	# these have resp as rho*log(Temp), rather than just an intercept – no prob, just need to pick one, easy to change once we get there
	rho <- mod[[1]][1] * length(irr) 
	iota <- mod[[1]][2]
	gpp <- iota + sum(irr)
	nep <- rho+gpp
	#other ways to get nep:
	# nep = rho + gpp
	# nep = sum(fitted(mod), na.rm=TRUE) # even if there are NA's in the response variable, they shouldn't be included in fitted() ....
	# sum(noflux.do.diff) # i think this should be the same as sum(fitted(mod)) b/c model residuals sum to 0 ... right?
	# also note that NEP is gpp+rho (rho is negative by this convention, which is consistent w/ Kalman, Bayes, mle – unsure of BK)

	results <- data.frame("GPP"=gpp, "R"=rho, "NEP"=nep) 
	attr(results, "lm.mod") <- mod
	return(results)
	
}