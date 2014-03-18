
KFnllDO <- function(Params, do.obs, do.sat, K, Zmix, irr, wtr){
	
	# ===========================
	# = Unpack and set initials =
	# ===========================
	#!Pseudocode #1: Initial guesses for B, C, and Q t
	c1 <- Params[1] #PAR coeff
	c2 <- Params[2] #log(Temp) coeff
	Q <- exp(Params[3]) # Variance of the process error
	H <- exp(Params[4]) # Variance of observation error
	
	# See KalmanDO_smooth.R comments for explanation of beta
	kz <- K/Zmix # K and Zmix are both vector of length nobs
	beta <- 1-kz # beta is a vector of length nobs
	
	# Set first true value equal to first observation
	alpha <- do.obs[1]#Let's give this model some starting values
	
	# Set process covariance, P, equal to Q
	P <- Q #starting value
	
	# Empty vector for nll's
	nlls <- double(length(do.obs))
		
	# ==================
	# = Main Recursion =
	# ==================
	for(i in 2:length(do.obs)){
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
		alpha <- beta[i-1]*alpha + c1*irr[i-1] + c2*log(wtr[i-1]) + kz[i-1]*do.sat[i-1]
		P <- (beta[i-1]*P*beta[i-1]) + Q
		
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
		
		# =================
		# = Calculate NLL =
		# =================
		nlls[i] <- 0.5*log(2*pi) + 0.5*log(Eff)+ 0.5*eta*eta/Eff
		} # End recursion
		
		
	return(sum(nlls)) # return the sum of nll's
	}#End function
	