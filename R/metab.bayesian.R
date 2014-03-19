# library("rjags")
require("R2jags")
require("R2WinBUGS")


# =========================
# = Define the Jags Model =
# =========================
j2Mod <- function(){
	# model process
	for(i in 2:N){	
		Y[i] ~ dnorm(a[i], tauV) # observations (Y) are distributed with mean equivalent to true values, and precision tauV (1/tauV is variance of observation error)
		K[i] ~ dnorm(kP[i-1, 1], 1/kP[i-1, 2]) #distributin on K
		Uk[i] <- (K[i]*(satO[i-1] - a[i-1]))/Zmix[i-1] # exchange
		aHat[i] <- a[i-1] + U[i,]%*%C + Uk[i] #U[i,]%*%C #+ Uk # the process
		a[i] ~ dnorm(aHat[i], tauW) # true values have a mean equivalent to estimated values, but accompanied by process error (process precision is tauW)
	}

	# K <- kP[1:N,1] #distributin on K

	# Starting values
	# a[1] <- Y[1] # set the first true value equal to first observed value
	a[1] <- a0
	C[1] ~ dnorm(cP[1,1], 1/cP[1,2])
	C[2] ~ dnorm(cP[2,1], 1/cP[2,2])
	#Priors on regression coefficients
	# C[1,1] ~ dnorm(cP[1,1], 1/cP[1,2]) # GPP coefficient
	# C[2,1] ~ dnorm(cP[2,1], 1/cP[2,2]) # R coefficient
	# C[3,1] <- 1 # dummy 1 for exchange

	#Prior on errors
	tauV ~ dgamma(1.0E-3, 1.0E-3)
	tauW ~ dgamma(1.0E-3, 1.0E-3)
	sigmaV <- 1/sqrt(tauV)
	sigmaW <- 1/sqrt(tauW)	
	}

	modfile = tempfile('j2mod')
	write.model(j2Mod, modfile)





	# ================================
	# = Supply Data and run bayesFit =
	# ================================
	bayesFit <- function(data, params, tend="median", ...){ #function that writes jags model, traces params, supplies data, etc

	jags.m <- jags(data, NULL, parameters.to.save=params, modfile, n.chains=3, n.iter=5E3, n.burnin=5E2)

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

	# apply(jags.m$BUGSoutput$sims.matrix, 2, tF, tend="median")
	# apply(jags.m$BUGSoutput$sims.matrix, 2, tF, tend="mean")
	# apply(jags.m$BUGSoutput$sims.matrix, 2, tF, tend="mode1")
	# apply(jags.m$BUGSoutput$sims.matrix, 2, tF, tend="mode2")

	#Figure out the order of the sims.matrix columns ...
	GPP <- sum(ctSim[1]*data$U[,1]) # gpp coef * par, then sum
	R <- sum(ctSim[2]*data$U[,2]) # r coef * log(temp), then sum

	GPPsd <- sqrt(sum(sdSim[1]^2*data$U[,1]^2))
	Rsd <- sqrt(sum(sdSim[2]^2*data$U[,2]^2))

	return(list(
		"model"=jags.m, 
		"params"=ctSim[1:2], 
		"metab"=matrix(c(GPP,GPPsd,R,Rsd), nrow=2, dimnames=list(c("mu", "sd"), c("GPP", "R")))
	)) # need to clean up format, and maybe include a return of the sd's of the estimates
	}



	metab.bayesian = function(do.obs, do.sat, k.gas, z.mix, date.times, irr, wtr, priors=c("gppMu"=0, "gppSig2"=1E5, "rMu"=0, "rSig2"=1E5, "kSig2"=NA)){
 
 
	 if(!all(c(is.numeric(do.obs), is.numeric(do.sat), is.numeric(k.gas), 
	           is.numeric(z.mix), is.numeric(irr), is.numeric(wtr)))){
   
	   error('All inputs to metab.bayes must be numeric vectors')
	 }

	# ===========================================
	# = Define objects to be used in jags model =
	# ===========================================
	#Supply elements of U (PAR, log(temp), 0's)
	U <- matrix(NA, nrow=length(irr), ncol=2)
	U[,1] <- irr #data0[,"PAR"] # filler PAR Values
	U[,2] <- log(wtr) # log(data0[,"Temp"]) # filler log(temp) values

	#Supply kP (for K), and cP (for C)
	# ** I am making up these values for now
	kP <- matrix(NA, nrow=length(k.gas), ncol=2)
	# variances for K
	kP[,1] <- k.gas # means for K
	if(is.na(priors["kSig2"])){
		kP[,2] <- sum(k.gas)/length(k.gas)*0.1
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

	output <- bayesFit(data, params)
	return(output)
  
}





