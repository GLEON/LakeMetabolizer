

# ====================================
# = Function to write the jags model =
# ====================================
bayes.makeModel <- function(k.gas){
  
	if(!require("R2jags") | !require("R2WinBUGS")){
    stop('metab.bayesian requires R2jags and R2WinBUGS')
	}
	
	finite.1oK <- is.finite(1/k.gas)
	
	choose.allK <- all(finite.1oK)
	choose.noK <- all(!finite.1oK)
	choose.bothK <- any(finite.1oK) & any(!finite.1oK)
	
	choice.mod <- c("allK", "noK", "bothK")[c(choose.allK, choose.noK, choose.bothK)]

	# Write the appropriate bayesian model into a temporary file
	modfile <- tempfile('jags.metab.bayes')
	switch(choice.mod,
		allK = write.model(bayes.mod.allK, modfile),
		noK = write.model(bayes.mod.noK, modfile),
		bothK = write.model(bayes.mod.bothK, modfile)
	)
	return(modfile)
}


# ==================================
# = Bayes model for all non-zero K =
# ==================================
bayes.mod.allK <- function(){
	
	# model process
	for(i in 2:N){	
		Y[i] ~ dnorm(a[i], tauV) # observations (Y) are distributed with mean equivalent to true values, and precision tauV (1/tauV is variance of observation error)
		K[i-1] ~ dnorm(kP[i-1, 1], 1/kP[i-1, 2]) #distributin on K
		
		kz[i-1] <- K[i-1]/Zmix[i-1]
		
		a1[i] <- U[i-1,]%*%C + kz[i-1]*satO[i-1]
		aHat[i] <- a1[i]/kz[i-1] + -exp(-1*kz[i-1])*a1[i]/kz[i-1] + exp(-kz[i-1])*a[i-1]
		
		a[i] ~ dnorm(aHat[i], tauW) # true values have a mean equivalent to estimated values, but accompanied by process error (process precision is tauW)
	}
	
	# Starting values
	a[1] <- a0
	
	#Priors on regression coefficients
	C[1] ~ dnorm(cP[1,1], 1/cP[1,2])
	C[2] ~ dnorm(cP[2,1], 1/cP[2,2])
	
	#Prior on errors
	tauV ~ dgamma(1.0E-3, 1.0E-3)
	tauW ~ dgamma(1.0E-3, 1.0E-3)
	sigmaV <- 1/sqrt(tauV)
	sigmaW <- 1/sqrt(tauW)	
}

# ===============================
# = Bayes model w/ K turned off =
# ===============================
bayes.mod.noK <- function(){

	# model process
	for(i in 2:N){	
		Y[i] ~ dnorm(a[i], tauV) # observations (Y) are distributed with mean equivalent to true values, and precision tauV (1/tauV is variance of observation error)
		# NOT USED
		K[i-1] ~ dnorm(kP[i-1, 1], 1/kP[i-1, 2]) #distribution on K
		
		aHat[i] <- a[i-1] + U[i-1,]%*%C # the process

		a[i] ~ dnorm(aHat[i], tauW) # true values have a mean equivalent to estimated values, but accompanied by process error (process precision is tauW)
	}

	# Starting values
	a[1] <- a0

	#Priors on regression coefficients
	C[1] ~ dnorm(cP[1,1], 1/cP[1,2])
	C[2] ~ dnorm(cP[2,1], 1/cP[2,2])


	#Prior on errors
	tauV ~ dgamma(1.0E-3, 1.0E-3)
	tauW ~ dgamma(1.0E-3, 1.0E-3)
	sigmaV <- 1/sqrt(tauV)
	sigmaW <- 1/sqrt(tauW)	
}

# ===============================================
# = Bayes model that can handle both 0 and !0 K =
# ===============================================
bayes.mod.bothK <- function(){
	
	# model process
	for(i in 2:N){	
		Y[i] ~ dnorm(a[i], tauV) # observations (Y) are distributed with mean equivalent to true values, and precision tauV (1/tauV is variance of observation error)
		K[i-1] ~ dnorm(kP[i-1, 1], 1/kP[i-1, 2]) #distributin on K
		
		# jags cannot handle an expression that includes division by 0 (so can't do blah <- ifelse(kz==0, 1, 1/kz), b/c it'll do 1/kz even when kz==0)
		# so to avoid division by 0 when k.gas==0, have to make a "safe" kz that is 1 if kz is 0
		# when kzSafe is 1 (to avoid division by 0), that means that we just want aHat to be the bio process
		# so we have to cancel out all math that is done by a kzSafe==1 (b/c it is bogus) by multiplying by 0, and add bio process
		# but if kzSafe!=0, we should multiply the math involving kzSafe by 1 (to keep it), and multiply the added bio process by 0 (to remove it)
		# this is a pretty hacked solution, but I don't see another choice (the only control flow in jags is ifelse, no if(){}else{})
			
		kz[i-1] <- K[i-1]/Zmix[i-1]
		kzSafe[i-1] <- ifelse(kz[i-1]==0, 1, kz[i-1]) # if kz is 0, change to 1
		kzCancel[i-1] <- ifelse(kzSafe[i-1]==1, 0, 1) # if kzSafe is 1 (meaning kz is 0), kzCancel needs to be 0
		pCancel[i-1] <- ifelse(kzSafe[i-1]==1, 1, 0) # if kzSafe is just kz, then math involving kzSafe is not bogus, and need to cancel the added process
		
		a1[i] <- U[i-1,]%*%C + kz[i-1]*satO[i-1]
		aHat[i] <- (a1[i]/kzSafe[i-1] + -exp(-1*kzSafe[i-1])*a1[i]/kzSafe[i-1] + exp(-kzSafe[i-1])*a[i-1])*kzCancel[i-1] + (a[i-1] + a1[i])*pCancel[i-1]

		a[i] ~ dnorm(aHat[i], tauW) # true values have a mean equivalent to estimated values, but accompanied by process error (process precision is tauW)
	}
	

	# Starting values
	a[1] <- a0
	
	#Priors on regression coefficients
	C[1] ~ dnorm(cP[1,1], 1/cP[1,2])
	C[2] ~ dnorm(cP[2,1], 1/cP[2,2])

	#Prior on errors
	tauV ~ dgamma(1.0E-3, 1.0E-3)
	tauW ~ dgamma(1.0E-3, 1.0E-3)
	sigmaV <- 1/sqrt(tauV)
	sigmaW <- 1/sqrt(tauW)	
}


# ================================
# = Supply Data and run bayesFit =
# ================================
bayesFit <- function(data, params, mf, tend="median", ...){ #function that writes jags model, traces params, supplies data, etc
	
	bf.args <- list(...)
	
	jags.m <- jags(data, NULL, parameters.to.save=params, mf)

	tF <- function(x, tend){ # tendency function
		switch(tend,
			median=median(x), #median
			mean=mean(x), #mean
			mode1 = unique(x)[which.max(tabulate(match(x, unique(x))))], # mode --- most frequently observed value
			mode2 = { # mode --- based on the highest peak of the posterior probability (from density plot)
				xd <- density(x)
				xd$x[which.max(xd$y)]
			}
			)
	}
	# medSim <- matrix(apply(jags2.m$BUGSoutput$sims.matrix, 2, median)[-(1)], nrow=115, ncol=3)
	simOut <- jags.m$BUGSoutput$sims.matrix
	ctSim <- apply(simOut, 2, tF, tend) # the central tendency metric
	sdSim <- apply(simOut, 2, sd)

	#Figure out the order of the sims.matrix columns ...
	n.obs <- length(data$U[,1])
	GPP <- mean(ctSim[1]*data$U[,1], na.rm=TRUE) * n.obs # gpp coef * par, then sum
	R <- mean(ctSim[2]*data$U[,2], na.rm=TRUE) * n.obs # r coef * log(temp), then sum

	GPPsd <- sqrt(sum(sdSim[1]^2*data$U[,1]^2))
	Rsd <- sqrt(sum(sdSim[2]^2*data$U[,2]^2))
	NEPsd <- GPPsd^2 + Rsd^2

	return(list(
		"model" = jags.m, 
		"params" = ctSim[1:2], 
		"metab.sd" = matrix(c(GPPsd, Rsd, NEPsd), nrow=1, dimnames=list(NULL, c("GPPsd", "Rsd", "NEPsd"))),
		"metab" = matrix(c(GPP, R, GPP+R), nrow=1, dimnames=list(NULL, c("GPP", "R", "NEP")))
	)) # need to clean up format, and maybe include a return of the sd's of the estimates
}


#'@title
#'Metabolism model based on a bayesian parameter estimation framework
#'@description
#'This function runs the bayesian metabolism model on the supplied gas concentration and 
#'other supporting data. This allows for both estimates of metabolism along with uncertainty around the parameters.
#'@usage
#'metab.bayesian(do.obs, do.sat, k.gas, z.mix, irr, wtr, ...)
#'@param do.obs Vector of dissovled oxygen concentration observations, mg L^-1
#'@param do.sat Vector of dissolved oxygen saturation values based on water temperature. Calculate using \link{o2.at.sat}
#'@param k.gas Vector of kGAS values calculated from any of the gas flux models (e.g., \link{k.cole}) and converted to kGAS using \link{k600.2.kGAS}
#'@param z.mix Vector of mixed-layer depths in meters. To calculate, see ts.meta.depths
#'@param irr Vector of photosynthetically active radiation in umoles/m2/s
#'@param wtr Vector of water temperatures in deg C. Used in scaling respiration with temperature
#'@param ... Parameter priors supplied as a named list
#'@return A named list of parameter estimates. 
#'\item{GPP}{Estimated Gross Primary Productivity and R	Estimated ecosystem respiration}
#'\item{R}{Estimated ecosystem respiration}
#'
#'@author Luke Winslow, Ryan Batt
#'@references Holtgrieve, Gordon W., Daniel E. Schindler, Trevor a. Branch, and Z. 
#'Teresa A'mar. 2010. \emph{Simultaneous Quantification of Aquatic Ecosystem Metabolism and 
#'Reaeration Using a Bayesian Statistical Model of Oxygen Dynamics}. 
#'Limnology and Oceanography 55 (3): 1047-1062. doi:10.4319/lo.2010.55.3.1047.
#'@seealso 
#'\link{metab.mle}
#'\link{metab.bookkeep}
#'@examples
#'library(rLakeAnalyzer)
#'\dontrun{
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
#'
#'k600 = k.cole.base(wnd[,2])
#'k.gas = k600.2.kGAS.base(k600, wtr[,3], 'O2')
#'do.sat = o2.at.sat(wtr[,3], altitude=300)
#'priors <- list("gppMu"=0, "gppSig2"=1E5, "rMu"=0, "rSig2"=1E5, "kSig2"=NA)
#'
#'metab.bayesian(irr=irr[,2], z.mix=rep(1, length(k.gas)), 
#'               do.sat=do.sat, wtr=wtr[,2],
#'               k.gas=k.gas, do.obs=doobs[,2],priors)
#'
#'}
#'@export
metab.bayesian = function(do.obs, do.sat, k.gas, z.mix, irr, wtr, ...){
	
	mb.args <- list(...)
	
	# Check for priors supplied in ...
	if("priors" %in% names(mb.args)){
		t.priors <- mb.args$priors
		p.name.logic <- all(c(c("gppMu", "gppSig2", "rMu", "rSig2", "kSig2"))%in%names(t.priors))
		p.class.logic <- is.integer(t.priors) | is.numeric(t.priors)
		if(p.name.logic & p.class.logic){
			priors <- t.priors
		}else{
			stop("supplied priors was not a (properly) named numeric/integer vector")
		}
	}else{
		priors <- c("gppMu"=0, "gppSig2"=1E5, "rMu"=0, "rSig2"=1E5, "kSig2"=NA)
	}
		
	 if(!all(c(is.numeric(do.obs), is.numeric(do.sat), is.numeric(k.gas), 
	           is.numeric(z.mix), is.numeric(irr), is.numeric(wtr)))){
   
	   stop('All inputs to metab.bayes must be numeric vectors')
	 }
	
	require("R2jags")
	require("R2WinBUGS")
	
	# Define model and write to file
	# Model choice depends on k values (all 0, all non-0, mixture)
	modfile <- bayes.makeModel(k.gas=k.gas) # k.gas would be available at higher scope, but passing to be explicit. 
	
	# ===========================================
	# = Define objects to be used in jags model =
	# ===========================================
	#Supply elements of U (PAR, log(temp))
	U <- matrix(NA, nrow=length(irr), ncol=2)
	U[,1] <- irr # PAR Values
	U[,2] <- log(wtr) # log(temp) values

	#Supply kP (for K), and cP (for C)
	kP <- matrix(NA, nrow=length(k.gas), ncol=2)
	# variances for K
	kP[,1] <- k.gas # means for K
	if(is.na(priors["kSig2"])){
		k0.logic <- !is.finite(1/kP[,1]) # test for when k is 0
		kP[,2] <- sum(k.gas)/sum(!k0.logic)*0.1 # k variance = mean of the non-zero K, times 0.1
		kP[k0.logic,2] <- 1E-9
	}else{
		kP[,2] <- priors["kSig2"]
	}
	
	cP <- matrix(NA, nrow=2, ncol=2)
	cP[1,1] <- priors["gppMu"] # prior mean of GPP coefficient (C[1,1]*PAR=GPP)
	cP[1,2] <- priors["gppSig2"] # prior variance of GPP coefficient
	cP[2,1] <- priors["rMu"] # prior mean of R coefficient (C[2,1]*log(Temp)=R)
	cP[2,2] <- priors["rSig2"] # prior variance of R coefficient


	# Put in final format supplied to jags
	data <- list(Y=do.obs, N=length(do.obs), U=U, kP=kP, cP=cP, satO=do.sat, a0=do.obs[1], Zmix=z.mix)
	params <- c("C", "sigmaV", "sigmaW")

	output <- bayesFit(data, params, mf=modfile)
	return(output)
	
}





