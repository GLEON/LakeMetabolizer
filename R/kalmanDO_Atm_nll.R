#6Feb2011, added warning message
#18Jan2011, changed it back to Alpha.
#17Jan2011, changed Alpha to DO[i] when calculating diffusion... trying to figure out why diffusion is so high, and respiration is so low.  
#In simulations, the DO concentration just shoots up over observed values when I use the model presented here and the parameter estimates determined from the KF.  
#Thus, it seems the KF is overestimating GPP relative to R, and is attributing too much of the reduction in DO to diffusion and not respiration.
KFnllDO <- function(Params, DO, Aitch, PAR, Chla, Temp, Atm=FALSE, Wind=NA, Freq=360, Zmix, Press){
	
	#!Pseudocode #1: Initial guesses for B, C, and Q t
	Beta <- 1
	LightGPPCoef <- Params[1] #PAR coeff; in Sea, for GPP
	TempRCoef <- Params[2] #log(Temp) coeff; in Sea, for R
	

  Sea <- matrix(c(LightGPPCoef,TempRCoef,1),nrow=1,ncol=3,byrow=TRUE) #Collecting coefficients to be fit into a matrix.  
    #Elements that are a part of Params will be fit with nlm()

	Queue <- Params[3]#Variance of the process error
	
	Kt <- KO2(Temp,Freq,Wind)
    RefT <- SatdConc(Temp,Press) #I removed multiplying by .942 because it seems that this has already been done in the sonde's calculations (when the sondes are calibrated in the bubble bath, they go to about 94.2%)
	
	Ewe <- matrix(nrow=length(DO), ncol=3)
	Ewe[,1] <- PAR
	Ewe[,2] <- Temp
	#Ewe[,3] <- (Kt * (RefT-DO))/Zmix #maybe not best to use DO here
	
		
	
	#!Pseudocode #2: Set first true value equal to first observation
	Alpha <- DO[1]#Let's give this model some starting values
	
	#!Pseudocode #3: Set process covariance, P, equal to Q
	Pea <- Queue #starting value
	
	NegLogLikes <- c(0,rep(NA,(length(DO)-1)))#Fill in all likelihood values with 0's
		
	#!Pseudocode #4: Starting with 2nd time step, build a Time Series of Alpha and P
	for(i in 2:length(DO)){#Go through entire data set, save for the first point
	Ewe[i,3] <- ((Kt[i-1] * (RefT[i-1]-Alpha))/Zmix[i-1]) #      #Kt <- KO2(Temp[i],Freq,Wind[i])
      #RefT <- SatdConc(Temp[i],(0.942*760))
			#Ewe <- c(PAR[i], log(Temp[i]), (Kt * (RefT-Alpha))/Zmix[i])#Divide by Zmix here to change diffusion from an areal flux into a change in concentration, which is how everything else is calculated.
			#Light, Temp, Diffusion; Diffusion is just multiplied by 1 (via Sea)
		
		#!Pseudocode #4a: Predictions
		Alpha <- Beta*Alpha + Sea[1,1]*Ewe[i,1] + Sea[1,2]*log(Ewe[i,2]) + Sea[1,3]*Ewe[i,3] #(Beta%*%Alpha) + (Sea%*%Ewe)
		#!Added Kt to variance propagation 14Jul10
		Pea <- (Beta^2)*Pea*(Kt[i]^2) + Queue #(Beta*Kt%*%Pea%*%t(Beta*Kt)) + Queue
		
		#!Pseudocode #4b: Update Predictions
		Eta <- DO[i] - Alpha
		Eff <- Pea + Aitch
		
		#Finv_Pea <- solve(Eff,Pea)
		#Finv_Eta <- solve(Eff,Eta)
		
		Alpha <- Alpha + Pea/Eff*Eta #Alpha+ (Pea%*%Finv_Eta)
		Pea <- Pea - Pea*Pea/Eff #Pea- (Pea%*%Finv_Pea)
		
		#!Pseudocode #5: Calculate the NLL
		NegLogLikes[i] <- .5*log(2*pi) + .5*log(Eff)+ .5*Eta*Eta/Eff #0.5*log(2*pi) + 0.5*log(Eff) + .5*Eta*Finv_Eta #(0.5*t(Eta)%*%Finv_Eta)#Calculate the likelihood for this step
		}#End cycling through data// series
		
		
	NLL <- sum(NegLogLikes)#Sum up all likelihoods for the total for this model/ set of parameters
	return(NLL)#NLL is the item which nlm() will attemp to minimize, thus yielding the 'most likely' combinations of parameter values
	}#End function
	