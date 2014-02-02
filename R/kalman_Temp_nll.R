KFnllTemp <- function(Params,Temp, Watts, Aitch){
	NLL_Term1 <- 0.5*log(2*pi)
	Beta <- Params[1]
	Sea <- Params[2]
	Queue <- exp(Params[3]) #exp to force pos. check _chla too
	Ewe <- Watts
	
	Alpha <- Temp[1]
	Pea <- Queue
	
	NegLogLikes <- c(0,rep(NA,(length(Temp)-1)))
	
	for(i in 2:length(Temp)){		
		#Pseudocode #4a: Predictions
		Alpha <- Beta*Alpha + Sea*Ewe[i-1]
		Pea <- Beta*Pea*Beta + Queue
		
		#Pseudocode #4b: Update Predictions
		Eta <- Temp[i] - Alpha
		Eff <- Pea + Aitch
		
		Alpha <- Alpha + Pea/Eff*Eta
		Pea <- Pea - Pea*Pea/Eff
		
		
		
		#Pseudocode #5: Calculate the NLL
		NegLogLikes[i] <- NLL_Term1 + 0.5*Eta*Eta/Eff  + 0.5*log(Eff)    
		}
	NLL <- sum(NegLogLikes)
	return(NLL)

	}