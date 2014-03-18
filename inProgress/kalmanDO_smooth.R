KFsmoothDO <- function(Params, do.obs, do.sat, K, Zmix, irr, wtr, Hfac=NULL){
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
	kz <- K/Zmix # K and Zmix are both vector of length nobs
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
	
	# ==========
	# = System =
	# ==========
	# y[t] = Z[t]*alphaStar[t] + 
	
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

	