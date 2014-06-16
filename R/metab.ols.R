metab.ols <- function(do.obs, do.sat, k.gas, z.mix, irr, wtr, ...){

	n.obs <- length(do.obs)

	do.diff <- diff(do.obs)

	#basically the average of the flux at time T_0 and T_1
	#flux <- -0.5 * k * (D1 + D2 - 2*S)
	#
	# can be re-arranged
	# flux <- 0.5 * [k * (D1-S) + k * (d2-s)]
	# Average of each inst flux
	inst_flux <- k.gas * (do.sat - do.obs)  # positive is into the lake
	
	# flux <- (inst_flux[1:(n.obs-1)] + inst_flux[-1])/2
	flux <- inst_flux[-1]

	
	noflux.do.diff <- do.diff - flux/z.mix[-1]
	
	lntemp <- log(wtr)
	mod <- lm(noflux.do.diff ~ irr[-1] + lntemp[-1] -1) # note that if we decied to use resp=rho*log(Temp), you would do lm(do~irr+lntemp-1) (mind the -1)

	rho <- mod[[1]][2] 
	iota <- mod[[1]][1]
	mod.matrix <- model.matrix(mod)
	gpp <- mean(iota*mod.matrix[,1], na.rm=TRUE) * n.obs
	resp <- mean(rho*mod.matrix[,2], na.rm=TRUE) * n.obs
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