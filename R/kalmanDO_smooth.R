KFsmoothDO <- function(Params, DO, Aitch, PAR, DOsatd, KO2zmix, lnTemp){
	Nobs <- length(DO)
	Beta <- 1-KO2zmix #DO_t = 1*DO_t-1 + -KO2zmix*DO_t-1 + Sea%*%Ewe + Eta === (1-KO2zmix)*DO_t-1.....
	
	#Pseudocode #1: Initial guesses for C and Q (B isn't being fitted)
	Sea <- matrix(data=c(Params[1],Params[2], NA), nrow=3, byrow=TRUE) #Collecting coefficients to be fit into a matrix.  Elements that are a part of Params will be fit with nlm()  
	Queue <- Params[3]#Variance of the Process Error
	#!Pseudocode #2: Set first true value equal to first observation
	Alpha <- DO[1]#Let's give this model some starting values
	#!Pseudocode #3: Set process covariance, P, equal to Q
	Pea <- Queue #starting value
	
	AHatvec <- c(Alpha,rep(NA,(Nobs-1)))
	PHatvec <- c(Pea,rep(NA,(Nobs-1)))
	Avec <- c(Alpha,rep(NA,(Nobs-1)))
	Pvec <- c(Pea,rep(NA,(Nobs-1)))
	Etavec <- c(0,rep(NA,(Nobs-1)))

	Ewe <- matrix(data=c(PAR, lnTemp, DOsatd), ncol=3, byrow=FALSE)
	
	#!Pseudocode #4: Starting with 2nd time step, build a Time Series of Alpha and P
	for(i in 2:Nobs){ 
		Sea[3,1] <- KO2zmix[i-1]
		
		#Pseudocode #4a: Predictions
		Alpha <- Beta[i-1]*Alpha + Ewe[i-1,]%*%Sea
		AHatvec[i] <- Alpha
		Pea <- (Beta[i-1]*Pea*Beta[i-1]) + Queue
		PHatvec[i] <- Pea
	
		#Pseudocode #4b: Update Predictions
		Eta <- DO[i] - Alpha
		Eff <- Pea + Aitch
		Alpha <- Alpha + Pea/Eff*Eta #Alpha <- Alpha +P/F(DO[i]-Alpha) == a[t] <- a[t|t-1] + P[t|t-1]%*%t(Z[t])%*%F[t]^-1 *(y[t]-Z[t]*a[t|t-1]-d[t]), where a=Alpha, P=Pea, Z=autoregressive coefficient in the process equation which I don't have, F=Eff, y=DO, and d=deterministic control vector in the measurement equation, which I don't have (Harvey 1989, pg. 106)
		Pea <- Pea - Pea*Pea/Eff #Pea- (Pea%*%Finv_Pea)

		Avec[i] <- Alpha
		Pvec[i] <- Pea
		Etavec[i] <- Eta
	}
	
	#Kalman Smoother
	Asmooth <- rep(NA,Nobs)
	Psmooth <- rep(NA,Nobs)
	Asmooth[Nobs] <- Avec[Nobs]
	Psmooth[Nobs] <- Pvec[Nobs]
	
	for(i in (Nobs-1):1){
		Pstar <- Pvec[i]*Beta[i+1]/PHatvec[i+1]
		Asmooth[i] <- Avec[i] + Pstar*(Asmooth[i+1]-AHatvec[i+1])
		}
	return(list(Asmooth, Etavec))	
		
	}