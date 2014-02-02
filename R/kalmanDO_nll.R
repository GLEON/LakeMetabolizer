#This version not suited for forcing parameters to positive
#Updated to print "non-finite" diagnostics
#v8.0 has the artificial inflation of the likelihood when the optimizer guesses at negative values for Q (log(-x) = Inf)
KFnllDO <- function(Params, DO, Aitch, PAR, Chla, Temp, Atm=FALSE, Wind=NA, Freq=360){
	
	#Pseudocode #1: Initial guesses for B, C, and Q t
	Beta <- 1
	LightGPPCoef <- Params[1] #PAR coeff; in Sea, for GPP
	TempRCoef <- Params[2] #log(Temp) coeff; in Sea, for R

	Sea <- matrix(c(LightGPPCoef,TempRCoef),nrow=1,ncol=2,byrow=TRUE) #Collecting coefficients to be fit into a matrix.  Elements that are a part of Params will be fit with nlm()  
	Queue <- Params[3]
		
	
	#Pseudocode #2: Set first true value equal to first observation
	Alpha <- DO[1]#Let's give this model some starting values
	
	#Pseudocode #3: Set process covariance, P, equal to Q
	Pea <- Queue #starting value
	
	NegLogLikes <- c(0,rep(NA,(length(DO)-1)))
	
	Ewe <- matrix(nrow=length(DO),ncol=2)
	Ewe[,1] <- PAR
	Ewe[,2] <- Temp
		
	#Pseudocode #4: Starting with 2nd time step, build a Time Series of α and P
	for(i in 2:length(DO)){#Could just as well have been DO		
		#Pseudocode #4a: Predictions
		#Alpha <- Beta*Alpha + Sea[1,1]*Ewe[1] + Sea[1,2]*Ewe[2]
		Alpha <- Beta*Alpha + Sea[1,1]*Ewe[i-1,1] + Sea[1,2]*log(Ewe[i-1,2])
		#(Beta%*%Alpha) + (Sea%*%Ewe) #
		Pea <- (Beta^2)*Pea + Queue #(Beta%*%Pea%*%t(Beta)) + Queue #
		
		#Pseudocode #4b: Update Predictions
		Eta <- DO[i] - Alpha
		Eff <- Pea + Aitch
		if(Eff<0){
			return(1E10)
		}

		Alpha <- Alpha + Pea/Eff*Eta #Alpha+ (Pea%*%Finv_Eta)  # #Alpha <- Alpha +P/F(DO[i]-Alpha) == a[t] <- a[t|t-1] + P[t|t-1]%*%t(Z[t])%*%F[t]^-1 *(y[t]-Z[t]*a[t|t-1]-d[t]), where a=Alpha, P=Pea, Z=autoregressive coefficient in the process equation which I don't have, F=Eff, y=DO, and d=deterministic control vector in the measurement equation, which I don't have (Harvey 1989, pg. 106)
		#Pea <- Pea - (Pea*(Pea/Eff)) #Pea- (Pea%*%Finv_Pea)        # 
		Pea <- Pea*(1-(Pea/Eff))
		
		#Pseudocode #5: Calculate the NLL
		NegLogLikes[i] <- .5*log(2*pi) + log(Eff) + .5*Eta*Eta/Eff #0.5*log(2*pi) + 0.5*log(det(Eff)) + (0.5*t(Eta)%*%Finv_Eta)  # 
		}
	NLL <- sum(NegLogLikes)
	return(NLL)
	}