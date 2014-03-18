

metab.kalman <- function(do.obs, do.sat, k.gas, z.mix, date.times, irr, wtr, ...){
	# ==================
	# = Filter and Fit =
	# ==================
	# ==================
	# = Fit 1: optim() =
	# ==================
	guesses <- c(1E-4,1E-4,log(5),log(5))
	# KFnllDO(Params=c(1E-4,-1E-4,5,5), do.obs=data[,"do.obs"], do.sat=data[,"do.sat"], K=data[,"K"], Zmix=data[,"Zmix"], irr=data[,"irr"], wtr=data[,"wtr"])
	fit <- optim(guesses, fn=KFnllDO, do.obs=do.obs, do.sat=do.sat, k.gas=k.gas, z.mix=z.mix, irr=irr, wtr=wtr, ...)
	pars0 <- fit$par
	pars <- c(pars0[1], pars0[2], exp(pars0[3:4]))
	
	# ====================
	# = Fit 2: DEoptim() =
	# ====================
	# deLow <- c(0, -10, log(0.1), log(0.1))
	# deHigh <- c(10, 0, log(1E3), log(1E3))
	# fitDE <- DEoptim(fn=KFnllDO, lower=deLow, upper=deHigh, control=list(trace=FALSE), do.obs=do.obs, do.sat=do.sat, k.gas=k.gas, z.mix=z.mix, irr=irr, wtr=wtr)$optim
	# pars0DE <- fit$par
	# pars <- c(pars0[1], pars0[2], exp(pars0[3:4]))
	
	# ==========
	# = Smooth =
	# ==========
	smoothDO <- KFsmoothDO(pars, do.obs=do.obs, do.sat=do.sat, k.gas=k.gas, z.mix=z.mix, irr=irr, wtr=wtr)
	
	# ====================================
	# = Use fits to calculate metabolism =
	# ====================================
	GPP <- sum(pars[1]*data[,"irr"])
	R <- sum(pars[2]*log(data[,"wtr"]))
	
	return(list("smoothDO"=smoothDO,"params"=pars, "metab"=c("GPP"=GPP,"R"=R)))
}



# BayesAns <- metab.bayesian(do.obs=data[,"do.obs"], do.sat=data[,"do.sat"], k.gas=data[,"k.gas"], z.mix=data[,"z.mix"], irr=data[,"irr"], wtr=data[,"wtr"])
