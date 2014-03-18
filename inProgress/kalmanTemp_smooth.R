#             Function to Calculate a Smoothed Series t
# ****************************************************************
KFsmoothTemp <- function(Params, Temp, Watts, Aitch){
	N<- length(Temp)
	
	#Pseudocode #1: Inital guesses for B, C, and Q
	Beta <- Params[1]
	Sea <- Params[2]
	Queue <- exp(Params[3]) #forced Q to be pos to get rid of warnings; check kfmetab_chla and the NLL function.
	Ewe <- Watts
	#Pseudocode #3: Set process covariance, P, equal to Q
	Pea <- Queue
	
	#Pseudocode #2: Set first true value equal to first observation
	Alpha <- Temp[1]
	
	AHatvec <- c(Alpha,rep(NA,(N-1)))
	PHatvec <- c(Pea,rep(NA,(N-1)))
	Avec <- c(Alpha,rep(NA,(N-1)))
	Pvec <- c(Pea,rep(NA,(N-1)))
	
	Errors <- c()
	
	
	#Pseudocode #3: Set process covariance, P, equal to Q
	Pea <- Queue
	
	for(i in 2:length(Temp)){
		#Pseudocode #4a: Predictions
		Alpha <- Beta*Alpha + Sea*Ewe[i-1]
		AHatvec[i] <- Alpha
		Pea <- Beta*Pea*Beta + Queue
		PHatvec[i] <- Pea
		#Pseudocode #4b: Update Predictions
		Eta <- Temp[i]-Alpha
		Errors[i] <- Eta
		Eff <- Pea+Aitch
		Alpha <- Alpha + Pea/Eff*Eta
		Avec[i] <- Alpha
		Pea <- Pea - Pea*Pea/Eff
		Pvec[i] <- Pea
		}
	
	#Kalman Smoother
	Asmooth <- rep(NA,N)
	Psmooth <- rep(NA,N)
	Asmooth[N] <- Avec[N]
	Psmooth[N] <- Pvec[N]
	
	for(i in (N-1):1){
		Pstar <- Pvec[i]*Beta/PHatvec[i+1]
		Asmooth[i] <- Avec[i] + Pstar*(Asmooth[i+1]-AHatvec[i])
		}
	return(list(Asmooth, Errors))
}
