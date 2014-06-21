#'@title Metabolism model based on a maximum likelihood parameter estimation framework.
#'@description This function runs the maximum likelihood metabolism model on the supplied gas concentration and other supporting data. This is a 
#'common approach that allows for the concurrent estimation of metabolism paramters from a timeseries.
#'@param do.obs Vector of dissovled oxygen concentration observations, mg L^-1
#'@param do.sat Vector of dissolved oxygen saturation values based on water temperature. Calculate using \link{o2.at.sat}
#'@param k.gas Vector of kGAS values calculated from any of the gas flux models 
#'(e.g., \link{k.cole}) and converted to kGAS using \link{k600.2.kGAS}
#'@param z.mix Vector of mixed-layer depths in meters. To calculate, see \link{ts.meta.depths}
#'@param irr Vector of photosynthetically active radiation in umoles/m2/s
#'@param wtr Vector of water temperatures in deg C. Used in scaling respiration with temperature
#'@param priors Parameter priors supplied as a named list (example: c("gppMu"=0, "gppSig2"=1E5, "rMu"=0, "rSig2"=1E5, "kSig2"=NA))
#'@param ... additional arguments to be passed
#'@return
#'A named list of parameter estimates.
#'\item{GPP}{Estimated Gross Primary Productivity}
#'\item{R}{Estimated ecosystem respiration}
#'@author Luke A Winslow, Ryan Batt, GLEON Fellows
#'@references
#'Solomon CT, DA Bruesewitz, DC Richardson, KC Rose, MC Van de Bogert, PC Hanson, TK Kratz, B Larget, 
#'R Adrian, B Leroux Babin, CY Chiu, DP Hamilton, EE Gaiser, S Hendricks, V Istvanovics, A Laas, DM O'Donnell, 
#'ML Pace, E Ryder, PA Staehr, T Torgersen, MJ Vanni, KC Weathers, G Zhuw 2013. 
#'\emph{Ecosystem Respiration: Drivers of Daily Variability and Background Respiration in Lakes around the Globe}. 
#'Limnology and Oceanograph 58 (3): 849:866. doi:10.4319/lo.2013.58.3.0849.
#'@seealso
#'\link{metab.bayesian}, \link{metab.bookeep}, \link{metab.ols}
#'@examples
#'\dontrun{
#'library(rLakeAnalyzer)
#'doobs = load.ts(system.file('extdata', 
#'                            'Sparkling.doobs', package="LakeMetabolizer"))
#'wtr = load.ts(system.file('extdata', 
#'                          'Sparkling.wtr', package="LakeMetabolizer"))
#'wnd = load.ts(system.file('extdata', 
#'                          'Sparkling.wnd', package="LakeMetabolizer"))
#'irr = load.ts(system.file('extdata', 
#'                          'Sparkling.par', package="LakeMetabolizer"))
#'
#'#Subset a day
#'mod.date = as.POSIXct('2009-08-12')
#'doobs = doobs[trunc(doobs$datetime, 'day') == mod.date, ]
#'wtr = wtr[trunc(wtr$datetime, 'day') == mod.date, ]
#'wnd = wnd[trunc(wnd$datetime, 'day') == mod.date, ]
#'irr = irr[trunc(irr$datetime, 'day') == mod.date, ]
#'z.mix = ts.thermo.depth(wtr)
#'
#'k600 = k.cole.base(wnd[,2])
#'k.gas = k600.2.kGAS.base(k600, wtr[,3], 'O2')
#'do.sat = o2.at.sat.base(wtr[,3], altitude=300)
#'
#'metab.mle(doobs[,2], do.sat, k.gas, z.mix[,2], irr[,2], wtr[,3])
#'}
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
