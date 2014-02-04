# library("rjags")
library("R2jags")
library("R2WinBUGS")	


# =========================
# = Define the Jags Model =
# =========================
j2Mod <- function(){
	# model process
	for(i in 2:N){	
		Y[i] ~ dnorm(a[i], tauV) # observations (Y) are distributed with mean equivalent to true values, and precision tauV (1/tauV is variance of observation error)
		K ~ dnorm(kP[i-1, 1], 1/kP[i-1, 2]) #distributin on K
		U[i,3] <- (K*(satO[i-1] - a[i-1]))/Zmix[i-1] # exchange
		aHat[i] <- a[i-1] + U[i,]%*%C # the process
		a[i] ~ dnorm(aHat[i], tauW) # true values have a mean equivalent to estimated values, but accompanied by process error (process precision is tauW)
	}
	
	# Starting values
	a[1] <- Y[1] # set the first true value equal to first observed value
	
	#Priors on regression coefficients
	C[1,1] ~ dnorm(cP[1,1], 1/cP[1,2]) # GPP coefficient
	C[2,1] ~ dnorm(cP[2,1], 1/cP[2,2]) # R coefficient
	C[3,1] <- 1 # dummy 1 for exchange
	
	#Prior on errors
	tauV ~ dgamma(1.0E-3, 1.0E-3)
	tauW ~ dgmma(1.0E-3, 1.0E-3)
	sigmaV <- 1/sqrt(tauV)
	sigmaW <- 1/sqrt(tauW)	
}
write.model(j2Mod, "j2Mod")

# ================
# = Read in data =
# ================
#Supply elements of U (PAR, log(temp), 0's)
U[,1] <- rnorm(10) # filler PAR Values
U[,2] <- rnorm(10) # filler log(temp) values

#Supply kP (for K), and cP (for C)
# ** I am making up these values for now
kP[,1] <- rep(1E-3, 10) # means for K
kP[,2] <- rep(1,10) # variances for K
cP[1,1] <- 1E-5 # prior mean of GPP coefficient (C[1,1]*PAR=GPP)
cP[1,2] <- 1E-5 # prior variance of GPP coefficient
cP[2,1] <- -1E2 # prior mean of R coefficient (C[2,1]*log(Temp)=R)
cP[2,2] <- 1E2 # prior variance of R coefficient

#Supply observed DO
Y <- rnorm(10) # filler values

#Supply saturated DO
satO <- rnorm(10) # filler values

#Supply Zmix
Zmix <- rep(1,10) # filler values

data <- list(Y=Y, N=length(Y), U=U, kP=kP, cP=cP, satO=satO)
params <- c("C", "sigmaV", "sigmaW")


# ================================
# = Supply Data and run bayesFit =
# ================================
bayesFit <- function(data, params, ...){ #function that writes jags model, traces params, supplies data, etc

	jags.m <- jags(data, NULL, parameters.to.save=params, "j2Mod", n.chains=3, n.iter=2E3, n.burnin=5E2)
	
	# medSim <- matrix(apply(jags2.m$BUGSoutput$sims.matrix, 2, median)[-(1)], nrow=115, ncol=3)
	medSim <- apply(jags.m$BUGSoutput$sims.matrix, 2, median)
	sdSim <- apply(jags.m$BUGSoutput$sims.matrix, 2, sd)
	
	#Figure out the order of the sims.matrix columns ...
	GPP <- medSim[,1]%*%data$U[,1] # gpp coef * par, then sum
	R <- medSim[,2]%*%data$U[,2] # gpp coef * par, then sum
	
	return(c(GPP,R)) # need to clean up format, and maybe include a return of the sd's of the estimates
}




