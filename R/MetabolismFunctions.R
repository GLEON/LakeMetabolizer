
################################################################################
# Summary: Takes basic metabolism parameters and returns Predicted DO
#
# Authors: Gleon Student Fellows and Luke Winslow
# 
################################################################################
calcDoHat = function(iota, rho, doInit, irr, doSat, zMix, kO2){

  nObs = length(irr)
  #Set up output
  DOHat <- rep(NA,nObs)
  atmFlux <- rep(NA,nObs)  
  #Initialize DOHat
  DOHat[1] <- doInit
  
  #Calculate atmFlux and predicted DO for each time point
  #Fluxes out of lake have negative sign
  for (i in 1:(nObs-1)) {
    atmFlux[i] <- 	-kO2[i] * (DOHat[i] - doSat[i]) / zMix[i]
    DOHat[i+1] <- 	DOHat[i] + iota*irr[i] - rho + atmFlux[i]
  }
  
  return(DOHat)
}

################################################################################
# Summary: Takes metabolism parameters and DO and returns NLL model fit
#
# Authors: Gleon Student Fellows and Luke Winslow
# 
################################################################################
calcModelNLL <- function(iota, rho, doInit, irr, doSat, zMix, kO2, doObs){
  
  modeled = calcDoHat(iota, rho, doInit, irr, doSat, zMix, kO2)
  res 	= doObs - modeled
  
  nRes 	= length(res)
  SSE 	= sum(res^2)
  sigma2 	= SSE/nRes
  NLL 	= 0.5*((SSE/sigma2) + nRes*log(2*pi*sigma2))
  
  return(NLL)
}

################################################################################
# Summary:   Takes metabolism parameters and DO and returns NLL model fit
#
# Authors:   Luke Winslow
# 
# Reference: Solomon et al. 2012?
################################################################################
metab.bootstrap <- function(iota, rho, doInit, irr, doSat, zMix, kO2, doObs, timestep, n=1000, ar1.resids=FALSE){
  
  n.obs = length(doObs)
  
  doHat  = calcDoHat(iota, rho, doInit, irr, doSat, zMix, kO2)
  resids = doObs - doHat
  
  #If we are maintaining the ar1 component of the residuals, 
  # we must estimate ar1 coeff and the ar1 residual standard deviation
  if(ar1.resids){
    ar1.lm    = lm(resids[1:n.obs-1] ~ resids[2:n.obs]-1)
    ar1.coeff = ar1.lm$coefficients
    ar1.sd    = sd(ar1.lm$residuals)
  }
  
  toOptim <- function(iotaRhoDoinit){
    calcModelNLL(iotaRhoDoinit[1], iotaRhoDoinit[2], iotaRhoDoinit[3], irr, doSat, zMix, kO2, doSim)
  }
  
  #Pre-allocate the result data frame
  result <- data.frame(boot.iter = 1:n,
                       iota = rep(NA,n),
                       rho = rep(NA,n),
                       DOInit = rep(NA,n),
                       covergence = rep(NA,n),
                       nll = rep(NA,n),
                       GPP = rep(NA,n))
  
  for(i in 1:n){
    
    #Randomize the residuals using one of two methods
    if(ar1.resids){ #residual randomization keeping the ar1 data structure
      simRes = rep(NA, n.obs)
      simRes[1] = sample(resids,1)
      for(j in 2:n.obs){
        simRes[j] = ar1.coeff*simRes[j-1] + rnorm(n=1, sd=ar1.sd)
      }
       
    }else{ #Raw residual randomization
      #Randomize residuals without replacement
      simRes = sample(resids, length(resids), replace=FALSE) 
    }
    
    doSim = doHat + simRes
    
    #Run optim again with new simulated DO signal
    optimOut = optim(par = c(iota,rho,doInit), fn = toOptim)
    
    result[i,2:3] <- optimOut$par[1:2]*(1440/timestep) #iota and rho to units of mg O2 L-1 day-1 
    result[i,4] <- optimOut$par[3] #initial DO estimate 
    result[i,5] <- optimOut$convergence  #did model converge or not (0=yes, 1=no)
    result[i,6] <- optimOut$value #value of nll 
    result[i,7] <- (result[i,2]/(60*60*24))*sum(irr)*timestep*60 #GPP in units of mg O2 L-1 d-1 
    
  }
  
  return(result)
}



