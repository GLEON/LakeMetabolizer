#library("rjags")


# =========================
# = Define the Jags Model =
# =========================
bayes.makeModel <- function(){
	require("R2jags")
	require("R2WinBUGS")
	#Step 1: use a function to define the model
	bayes.mod <- function(){
		
		# model process
		for(i in 2:N){	
			Y[i] ~ dnorm(a[i], tauV) # observations (Y) are distributed with mean equivalent to true values, and precision tauV (1/tauV is variance of observation error)
			K[i-1] ~ dnorm(kP[i-1, 1], 1/kP[i-1, 2]) #distributin on K
			
			# old lines before gigantic workaround to incorporate the between-time-step calculus for K, 
			# which involves dividing by 0 when K is 0
			
			# Uk[i-1] <- (K[i-1]*(satO[i-1] - a[i-1]))/Zmix[i-1] # exchange
			# aHat[i] <- a[i-1] + U[i-1,]%*%C + Uk[i-1] #U[i,]%*%C #+ Uk # the process
			
			
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
		# a[1] ~ dnorm(Y[1], tauV) # the first true value is the firt observation + OE
		# aHat[1] ~ dnorm(a[1], tauW) # trying to say that the *first* estimate is the first observation + OE + PE, not sure if this is right
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
	
	# Step 2: write the bayesian model into a temporary file
	modfile <- tempfile('jags.metab.bayes')
	write.model(bayes.mod, modfile)
	return(modfile)
}




# ================================
# = Supply Data and run bayesFit =
# ================================
bayesFit <- function(data, params, mf, tend="median", ...){ #function that writes jags model, traces params, supplies data, etc
	
	jags.m <- jags(data, NULL, parameters.to.save=params, mf, n.chains=3, n.iter=2E3)

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



metab.bayesian = function(do.obs, do.sat, k.gas, z.mix, irr, wtr, priors=c("gppMu"=0, "gppSig2"=1E5, "rMu"=0, "rSig2"=1E5, "kSig2"=NA), ...){
	
	 if(!all(c(is.numeric(do.obs), is.numeric(do.sat), is.numeric(k.gas), 
	           is.numeric(z.mix), is.numeric(irr), is.numeric(wtr)))){
   
	   error('All inputs to metab.bayes must be numeric vectors')
	 }
	
	require("R2jags")
	require("R2WinBUGS")
	
	# Define model and write to file
	modfile <- bayes.makeModel() # might need to put this in bayesFit() if jags() can't find the model file. Trying here b/c it will be called less.
	
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

	output <- bayesFit(data, params, mf = modfile)
	return(output)
	
}





