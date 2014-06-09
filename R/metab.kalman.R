



# ======================
# = Kalman filter/ nll =
# ======================
# Main recursion written in C
KFnllDO <- function(Params, do.obs, do.sat, k.gas, z.mix, irr, wtr){
	
	# ===========================
	# = Unpack and set initials =
	# ===========================
	#!Pseudocode #1: Initial guesses for B, C, and Q t
	c1 <- Params[1] #PAR coeff
	c2 <- Params[2] #log(Temp) coeff
	Q <- exp(Params[3]) # Variance of the process error
	H <- exp(Params[4]) # Variance of observation error
	
	# See KalmanDO_smooth.R comments for explanation of beta
	kz <- k.gas/z.mix # K and Zmix are both vector of length nobs
	# beta <- 1-kz # beta is a vector of length nobs (this beta is for difference equation form)
	beta <- exp(-kz) # This beta is for using the differential equation form
	
	# Set first true value equal to first observation
	alpha <- do.obs[1]#Let's give this model some starting values
	
	# Set process covariance, P, equal to Q
	P <- Q #starting value
	
	# Empty vector for nll's
	nlls <- rep(0,length(do.obs))
		
	# ==================
	# = Main Recursion =
	# ==================
	# ===============
	# = Predictions =
	# ===============
	# Equations for Predictions from Harvey
	# a[t|t-1] = T[t]*a[t-1] + c[t] Harvey pg 105 eq. 3.2.2a
	# P[t|t-1] = T[t]*P[t-1]*T'[t] + R[t]*Q[t]*R'[t] Harvey pg 106 eq. 3.2.2b
	
	# Predictions where gas flux not split into beta etc. (I'm pretty sure this is wrong):
	# Uk <- K[i-1]*(do.sat[i-1] - alpha)/Zmix[i-1]
	# alpha <- alpha + c1*irr[i-1] + c2*log(wtr[i-1]) + Uk
	# P <- (Uk*P*Uk) + Q
	
	# Predictions where gas flux is split into beta:
	# Difference Equation Version (original):
	# alpha <- beta[i-1]*alpha + c1*irr[i-1] + c2*log(wtr[i-1]) + kz[i-1]*do.sat[i-1]
	
	# Differential Equation Version:
	# Gordon's code (from BayesMetabSS_indivProcessErr_072012.txt):
	# alpha[f,j] <- P[f,j] - rho[f,j] + KO2[f,j]* DOSat[f,j];
	# DOHat[f+1,j] <- (alpha[f,j]-(alpha[f,j]-KO2[f,j]*DOTrue[f,j])*exp(-KO2[f,j]))/KO2[f,j];
	
	# Define Gordon's "alpha" as "a1":
	# a1 = c1*irr + c2*log(wtr) + kz*do.sat # my a1 is his alpha
	
	# Coefficients in front of alpha are defined above as "beta", and play a particular role in propagating uncertainty
	# Need to algebraically regarrange Gordon's equation for DOHat (my "alpha" is his "DOHat")
	# Note that alpha and a1 are rewritten each iteration of the loop...
	# Therefore, alpha on the left side is "alpha[t]", and alpha on the right side is "alpha[t-1]"
	
	# alpha = (a1 - (a1 - kz*alpha)*exp(-kz))/kz # my alpha is his DOHat
	# alpha = (a1 + (-a1 + kz*alpha)*exp(-kz))/kz # redistribute -1
	# alpha = (a1 + exp(-kz)*(-a1 + kz*alpha))/kz # regarrange
	# alpha = (a1 + -exp(-kz)*a1 + exp(-kz)*kz*alpha)/kz # multiply "exp(kz)" through
	# alpha = a1/kz + -exp(-kz)*a1/kz + exp(-kz)*alpha # multiply "/kz" through
	
	# ======================
	# = Updating Equations =
	# ======================
	# Updating Equations from Harvey
	# a[t] = a[t|t-1] + P[t|t-1]*Z'[t]*F[t]^-1(y[t] - Z[t]*a[t|t-1] - d[t]). Harvey, page 106, 3.2.3a
	# P[t] = P[t|t-1] - P[t|t-1]*Z'[t]*F[t]^-1*Z[t]*P[t|t-1] Harvey, page 106, eq. 3.2.3b
	# F[t] = Z[t]*P[t|t-1]*Z'[t] + H[t] Harvey, page 106, eq. 3.2.3c
	# kalmanLoopC(double *alpha, double *doobs, double *c1, double *c2, double *P, double *Q, double *H,  double *beta, double *irr, double *wtr, double *kz, double *dosat, int *nobs)
	nlls <- kalmanLoopR(nlls=nlls, alpha=alpha, doobs=do.obs, c1=c1, c2=c2, P=P, Q=Q, H=H, beta=beta, irr=irr, wtr=wtr, kz=kz, dosat=do.sat)

	return(sum(nlls)) # return the sum of nll's
}#End function



# ===================
# = Kalman Smoother =
# ===================
KFsmoothDO <- function(Params, do.obs, do.sat, k.gas, z.mix, irr, wtr, Hfac=NULL){
	nobs <- length(do.obs)
	d0 <- double(nobs-1)
	# beta <- 1-KO2zmix #do.obs_t = 1*do.obs_t-1 + -KO2zmix*do.obs_t-1 + Sea%*%Ewe + eta === (1-KO2zmix)*do.obs_t-1.....
	
	# Unpack parameters (these were previously fitted)
	c1 <- Params[1] # irr Coeff
	c2 <- Params[2] # log(wtr) Coeff
	Q <- Params[3] # Variance of the Process Error
	if(is.null(Hfac)){
		H <- Params[4]
	}else{
		H <- Params[4]*Hfac
	}
	 # Variance of Observation Error
	
	# Need to define portion of K multiplied by state variable (DO)
	# Gas flux = K[t-1](do.sat[t-1] - alpha[t-1])/Zmix[t-1]
	# Gas flux = (K[t-1]/Zmix[t-1])*(do.sat[t-1]) - (K[t-1]/Zmix[t-1])*(alpha[t-1])
	# Note that K/Z is essentially a coefficient sitting in front of alpha, the estimate of DO
	# Therefore,
	# alpha[t] = 1*alpha[t-1] + c1*irr[t-1] + c2*log(wtr[t-1]) + K[t-1](do.sat[t-1] - alpha[t-1])/Zmix[t-1]
	# Becomes:
	# kz[t] = k[t]/Zmix[t]
	# alpha[t] = 1*alpha[t-1] + -kz[t-1]*alpha[t-1] + c1*irr[t-1] + c2*log(wtr[t-1]) + kz[t-1]*do.sat[t-1]
	# Or,
	# alpha[t] = (1-kz[t-1])*alpha[t-1] + c1*irr[t-1] + c2*log(wtr[t-1]) + kz[t-1]*do.sat[t-1]
	# Defining kz and redefining (1-kz[t]) as beta[t]:
	kz <- k.gas/z.mix # K and Zmix are both vector of length nobs
	# beta <- 1-kz # beta is a vector of length nobs (this beta is for difference equation form)
	beta <- exp(-kz) # This beta is for using the differential equation form
	
	# Set first true value equal to first observation
	alpha <- do.obs[1]
	
	# Set process covariance, P, equal to Q
	P <- Q # starting value
	
	# Initial values
	aHat <- c(alpha, d0) # aHat[t] == "a[t|t-1]" (estimate of a before updating)
	pHat <- c(P, d0) # pHat[t] == "p[t|t-1]" (estimate of a before updating)
	aVec <- aHat # aVec[t] == "a[t|t]" or "a[t]" (aVec is the "updated" version of aHat)
	pVec <- pHat # pVec[t] == "P[t|t]" or "P[t]" (pVec is the "updated" version of pHat)
	etaVec <- double(nobs)
	
	for(i in 2:nobs){
		# ===============
		# = Predictions =
		# ===============
		# Equations for Predictions from Harvey
		# a[t|t-1] = T[t]*a[t-1] + c[t] Harvey pg 105 eq. 3.2.2a
		# P[t|t-1] = T[t]*P[t-1]*T'[t] + R[t]*Q[t]*R'[t] Harvey pg 106 eq. 3.2.2b
		
		# Predictions where gas flux not split into beta etc.:
		# Uk <- K[i-1]*(do.sat[i-1] - alpha)/Zmix[i-1]
		# alpha <- alpha + c1*irr[i-1] + c2*log(wtr[i-1]) + Uk
		# aHat[i] <- alpha
		# P <- (Uk*P*Uk) + Q
		# pHat[i] <- P
		
		# Predictions where gas flux is split into beta (see explanation above):
		
		# Difference Equation Version:
		# alpha <- beta[i-1]*alpha + c1*irr[i-1] + c2*log(wtr[i-1]) + kz[i-1]*do.sat[i-1]
		
		# Differential Equation Version (see kalmanDO_nll.R for explanation):
		a1 <- c1*irr[i-1] + c2*log(wtr[i-1]) + kz[i-1]*do.sat[i-1]	
		alpha <- a1/kz[i-1] + -exp(-kz[i-1])*a1/kz[i-1] + beta[i-1]*alpha # NOTE: beta==exp(-kz); kz=K/Zmix
		
		aHat[i] <- alpha
		P <- (beta[i-1]*P*beta[i-1]) + Q
		pHat[i] <- P
	
		# ======================
		# = Updating Equations =
		# ======================
		# Updating Equations from Harvey
		# a[t] = a[t|t-1] + P[t|t-1]*Z'[t]*F[t]^-1(y[t] - Z[t]*a[t|t-1] - d[t]). Harvey, page 106, 3.2.3a
		# P[t] = P[t|t-1] - P[t|t-1]*Z'[t]*F[t]^-1*Z[t]*P[t|t-1] Harvey, page 106, eq. 3.2.3b
		# F[t] = Z[t]*P[t|t-1]*Z'[t] + H[t] Harvey, page 106, eq. 3.2.3c
				
		eta <- do.obs[i] - alpha
		Eff <- P + H
		alpha <- alpha + P/Eff*eta
		P <- P - P*P/Eff

		aVec[i] <- alpha
		pVec[i] <- P
		etaVec[i] <- eta
	}
	
	#Kalman Smoother
	aSmooth <- rep(NA,nobs)
	Psmooth <- rep(NA,nobs)
	aSmooth[nobs] <- aVec[nobs] # "starting" value for smoother (smoother starts at end and works backwards)
	# pSmooth[nobs] <- pVec[nobs]
	
	# Filtering is informed by past information
	# Smoothing includes the information from filtering (estimates of parameters), but also future information.
	# "The aim of filtering is to find the expected value of the state vector, alpha[t], conditional on the information available at time t, that is E(alpha[t]|Y[t]). The aim of smoothing is to take account of the information made available after time t. The mean of the distribution of alpha[t], conditional on all the sample, may be written as E(alpha[t]|Y[T]) and is known as the smoothed estimate. THe corresponding estimator is called the SMOOTHER. Since the smoother is based on more information than the filtered estimator, it will have a MSE which, in general, is smaller than that of the filtered estimator; it cannot be greater." ~ Harvey 1989, pgs 149-150.
	 #a[t|T] = a[t] + Pstar[t]*(a[t+1|T] - T[t+1]*a[t])
	# P[t|T] = P[t] + Pstar[t]*(P[t+1|T] - P[t+1|t])*Pstar[t]
	# Pstar[t] = P[t]*T[t+1]/P[t+1|t]
	# t is current time step, T is last time step (when in []), T contains AR parameters (when NOT in [])
	
	for(i in length(d0):1){
		pStar <- pVec[i]*beta[i+1]/pHat[i+1]
		aSmooth[i] <- aVec[i] + pStar*(aSmooth[i+1] - aHat[i+1])
		
		# CAN ALSO SMOOTH P, WHICH GIVES THE SMOOTHED COVARIANCE MATRIX (not a matrix for univariate; gives estimate of accuracy of state estimate)
		# pSmooth[i] <- pVec[i] + pStar*(pSmooth[i+1] - pHat[i+1])*pStar
		}
	return(aSmooth)	# return smoothed DO time series
}

# ===========================================
# = R function that calls the C loop for KF =
# ===========================================
# kalmanLoopC(double *alpha, double *doobs, double *c1, double *c2, double *P, double *Q, double *H,  double *beta, double *irr, double *wtr, double *kz, double *dosat, int *nobs)
kalmanLoopR <- function(nlls, alpha, doobs, c1, c2, P, Q, H, beta, irr, wtr, kz, dosat){
	nobs <- length(doobs)
	a.loop <- .C("kalmanLoopC", nlls=as.double(nlls), as.double(alpha), as.double(doobs), as.double(c1), as.double(c2), as.double(P), as.double(Q), as.double(H), as.double(beta), as.double(irr), as.double(wtr), as.double(kz), as.double(dosat), as.integer(nobs), PACKAGE="LakeMetabolizer")
	return(a.loop[["nlls"]])
}

# ====================
# = Wrapper function =
# ====================
# Kalman filter metabolism w/ main recursion in NLL function written in C
metab.kalman <- function(do.obs, do.sat, k.gas, z.mix, irr, wtr, ...){
	# Filter and fit
	guesses <- c(1E-4,1E-4,log(5),log(5))
	fit <- optim(guesses, fn=KFnllDO, do.obs=do.obs, do.sat=do.sat, k.gas=k.gas, z.mix=z.mix, irr=irr, wtr=wtr, ...)
	pars0 <- fit$par
	pars <- c("gppCoeff"=pars0[1], "rCoeff"=pars0[2], "Q"=exp(pars0[3]), "H"=exp(pars0[4]))
	
	# Smooth
	smoothDO <- KFsmoothDO(pars, do.obs=do.obs, do.sat=do.sat, k.gas=k.gas, z.mix=z.mix, irr=irr, wtr=wtr)
	
	# Use fits to calculate metabolism
	n.obs <- length(do.obs)
	GPP <- mean(pars[1]*irr, na.rm=TRUE) * n.obs
	R <- mean(pars[2]*log(wtr), na.rm=TRUE) * n.obs
	
	return(list("smoothDO"=smoothDO,"params"=pars, "metab"=c("GPP"=GPP,"R"=R, "NEP"=GPP+R)))
}
