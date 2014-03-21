
# Kalman filter metabolism w/ main recursion in NLL function written in C
metab.kalman2 <- function(do.obs, do.sat, k.gas, z.mix, irr, wtr, ...){
	# ==================
	# = Filter and Fit =
	# ==================
	guesses <- c(1E-4,1E-4,log(5),log(5))
	fit <- optim(guesses, fn=KFnllDO2, do.obs=do.obs, do.sat=do.sat, k.gas=k.gas, z.mix=z.mix, irr=irr, wtr=wtr, ...)
	pars0 <- fit$par
	pars <- c("gppCoeff"=pars0[1], "rCoeff"=pars0[2], "Q"=exp(pars0[3]), "H"=exp(pars0[4]))
	
	# ==========
	# = Smooth =
	# ==========
	smoothDO <- KFsmoothDO(pars, do.obs=do.obs, do.sat=do.sat, k.gas=k.gas, z.mix=z.mix, irr=irr, wtr=wtr)
	
	# ====================================
	# = Use fits to calculate metabolism =
	# ====================================
	GPP <- sum(pars[1]*irr, na.rm=TRUE)
	R <- sum(pars[2]*log(wtr), na.rm=TRUE)
	
	return(list("smoothDO"=smoothDO,"params"=pars, "metab"=c("GPP"=GPP,"R"=R)))
}
