
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
	# beta <- 1-kz # beta is a vector of length nobs (this beta is for difference equation form)
	beta <- exp(-kz) # This beta is for using the differential equation form
	
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
		# Difference Equation Version (original):
		# alpha <- beta[i-1]*alpha + c1*irr[i-1] + c2*log(wtr[i-1]) + kz[i-1]*do.sat[i-1]
		
		# Differential Equation Version:
		# Gordon's code (from BayesMetabSS_indivProcessErr_072012.txt):
		# alpha[f,j] <- P[f,j] - rho[f,j] + KO2[f,j]* DOSat[f,j];
		# DOHat[f+1,j] <- (alpha[f,j]-(alpha[f,j]-KO2[f,j]*DOTrue[f,j])*exp(-KO2[f,j]))/KO2[f,j];
		
		# Define Gordon's "alpha" as "a1":
		# a1 = c1*irr + c2*log(wtr) + kz*do.sat # my a1 is his alpha
		a1 <- c1*irr[i-1] + c2*log(wtr[i-1]) + kz[i-1]*do.sat[i-1]
		
		# Coefficients in front of alpha are defined above as "beta", and play a particular role in propagating uncertainty
		# Need to algebraically regarrange Gordon's equation for DOHat (my "alpha" is his "DOHat")
		# Note that alpha and a1 are rewritten each iteration of the loop...
		# Therefore, alpha on the left side is "alpha[t]", and alpha on the right side is "alpha[t-1]"
		
		# alpha = (a1 - (a1 - kz*alpha)*exp(-kz))/kz # my alpha is his DOHat
		# alpha = (a1 + (-a1 + kz*alpha)*exp(-kz))/kz # redistribute -1
		# alpha = (a1 + exp(-kz)*(-a1 + kz*alpha))/kz # regarrange
		# alpha = (a1 + -exp(-kz)*a1 + exp(-kz)*kz*alpha)/kz # multiply "exp(kz)" through
		# alpha = a1/kz + -exp(-kz)*a1/kz + exp(-kz)*alpha # multiply "/kz" through		
		alpha <- a1/kz[i-1] + -exp(-kz[i-1])*a1/kz[i-1] + beta[i-1]*alpha # NOTE: beta==exp(-kz); kz=K/Zmix
		
		# ==========================================
		# = Test speed & accuracy of alpha algebra =
		# ==========================================
		# f1 <- function(a1=2, kz=4, alpha=10){(a1 - (a1 - kz*alpha)*exp(-kz))/kz}
		# f2 <- function(a1=2, kz=4, alpha=10){(a1 + (-a1 + kz*alpha)*exp(-kz))/kz}
		# f3 <- function(a1=2, kz=4, alpha=10){(a1 + exp(-kz)*(-a1 + kz*alpha))/kz}
		# f4 <- function(a1=2, kz=4, alpha=10){(a1 + -exp(-kz)*a1 + exp(-kz)*kz*alpha)/kz}
		# f5 <- function(a1=2, kz=4, alpha=10){a1/kz + -exp(-kz)*a1/kz + exp(-kz)*alpha}
		# c(f1(), f2(), f3(), f4(), f5())
		# all(f1() == c(f2(), f3(), f4(), f5()))
		# library(microbenchmark)
		# microbenchmark(f1(), f2(), f3(), f4(), f5()) #f4() is slowest, f5() second slowest
		
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
	