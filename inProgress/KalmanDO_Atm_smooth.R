KFsmoothDO <- function(Params, DO, Aitch, PAR, Chla, Temp, Atm=FALSE, Wind=NA, Freq=360, Zmix, Press){
	N<- length(DO)#Length of the series
	
	#!Pseudocode #1: Initial guesses for B, C, and Q t
	Beta <- 1
	LightGPPCoef <- Params[1] #PAR coeff; in Sea, for GPP
	TempRCoef <- Params[2] #log(Temp) coeff; in Sea, for R
	
  Sea <- matrix(c(LightGPPCoef,TempRCoef,1),nrow=1,ncol=3,byrow=TRUE) #Collecting coefficients to be fit into a matrix.  
    #Elements that are a part of Params will be fit with nlm()
	Queue <- Params[3]#Variance of the Process Error

	#!Pseudocode #2: Set first true value equal to first observation
	Alpha <- DO[1]#Let's give this model some starting values
	
	#!Pseudocode #3: Set process covariance, P, equal to Q
	Pea <- Queue #starting value
	
	AHatvec <- c(Alpha,rep(NA,(N-1)))
	PHatvec <- c(Pea,rep(NA,(N-1)))
	Avec <- c(Alpha,rep(NA,(N-1)))
	Pvec <- c(Pea,rep(NA,(N-1)))
	
	FluxMatrix <- matrix(nrow=length(DO), ncol=4, byrow=TRUE)
	FluxMatrix[1,] <- rep(0,4)#No metabolism estimates for first data point
	#KMat <- matrix(nrow=length(DO), ncol=4, byrow=TRUE)
#	KMat[1,] <- rep(0,4)
#	colnames(KMat) <- c("Kt", "RefT", "Alpha", "Diff")
#	
#	Errors <- c()
	
	Kt <- KO2(Temp,Freq,Wind)
    RefT <- SatdConc(Temp,Press) #I removed multiplying by .942 because it seems that this has already been done in the sonde's calculations (when the sondes are calibrated in the bubble bath, they go to about 94.2%)

	
	Ewe <- matrix(nrow=length(DO), ncol=3)
	Ewe[,1] <- PAR
	Ewe[,2] <- Temp
	#Ewe[,3] <- (Kt * (RefT-DO))/Zmix #maybe not best to use DO here
	
	Gradient <- c()
	
	#!Pseudocode #4: Starting with 2nd time step, build a Time Series of Alpha and P
	for(i in 2:length(DO)){ 
		Ewe[i,3] <- ((Kt[i-1] * (RefT[i-1]-Alpha))/Zmix[i-1]) #
#      Kt <- KO2(Temp[i],Freq,Wind[i])
#      KMat[i,1] <- Kt 
#      
#      RefT <- SatdConc(Temp[i],(0.942*760))
#      KMat[i,2] <- RefT
#      
#      Ewe <- c(PAR[i], log(Temp[i]), (Kt * (RefT-Alpha))/Zmix[i])
#      KMat[i,4] <- Ewe[3]
#      
#      KMat[i,3] <- Alpha
		Gradient[i] <- RefT[i-1] - Alpha
		
		#!Pseudocode #4a: Predictions
		#Fill in matrix to track biological and physical fluxes of DO
		FluxMatrix[i,1] <- Sea[1,1]*Ewe[i,1] #GPP
		FluxMatrix[i,2] <- Sea[1,2]*log(Ewe[i,2]) #R
		FluxMatrix[i,4] <- Sea[1,3]*Ewe[i,3] #Diffusion
		
		#OK, real predictions.  
    #Can probably just row sum FluxMatrix[i,] instead of Sea%*%Ewe.
		Alpha <- Beta*Alpha + FluxMatrix[i,1] + FluxMatrix[i,2] + FluxMatrix[i,4]
		#Alpha <- Beta*Alpha + Sea[1,1]*Ewe[1] + Sea[1,2]*Ewe[2] + Sea[1,3]*Ewe[3] #(Beta%*%Alpha) + (Sea%*%Ewe)
		AHatvec[i] <- Alpha
		Pea <- (Beta^2)*Pea*(Kt[i]^2) + Queue #(Beta%*%Pea%*%t(Beta)) + Queue
		PHatvec[i] <- Pea
		
		#!Pseudocode #4b: Update Predictions
		Eta <- DO[i] - Alpha
		#Errors[i] <- Eta
		Eff <- Pea + Aitch
		
		#Finv_Pea <- solve(Eff,Pea)
		#Finv_Eta <- solve(Eff,Eta)
		
		Alpha <- Alpha + Pea/Eff*Eta #Alpha+ (Pea%*%Finv_Eta)
		Avec[i] <- Alpha
		Pea <- Pea - Pea*Pea/Eff #Pea- (Pea%*%Finv_Pea)
		Pvec[i] <- Pea
		}
		
	#!Kalman Smoother
	Asmooth <- rep(NA,N)
	Psmooth <- rep(NA,N)
	Asmooth[N] <- Avec[N]
	Psmooth[N] <- Pvec[N]
	
	for(i in (N-1):1){
		Pstar <- Pvec[i]*Beta/PHatvec[i+1]
		Asmooth[i] <- Avec[i] + Pstar*(Asmooth[i+1]-AHatvec[i+1])
		}
	return(list(Asmooth, FluxMatrix, Gradient))	
		
}
	