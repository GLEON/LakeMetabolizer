metab.optim.ts = function(do.obs.ts, irr.ts, do.sat.ts, z.mix.ts, k.gas.ts, ...){
  
  
}

################################################################################
# Summary: Fits a simple metabolism model using Negative Log Likelihood
#
# Authors: Luke Winslow
# 
################################################################################
metab.optim = function(do.obs, irr, do.sat, z.mix, k.gas, timestep){

  n.obs = length(do.obs)
  
  if(n.obs != length(irr) || n.obs != length(do.sat) || n.obs != length(z.mix) || n.obs != length(k.gas)){
    stop('All input data to metab.optim must be the same length')
  }
  
  if(n.obs < 4){
    stop('You must supply at least 4 observations to metab.optim.')
  }
  
  
  to.optim <- function(iotaRhoDoinit){
    calc.do.nll(iotaRhoDoinit[1], iotaRhoDoinit[2], iotaRhoDoinit[3], irr, do.sat, z.mix, k.gas, do.obs)
  }
  
  rho     = range(do.obs)/n.obs
  iota    = rho
  doInit  = do.obs[1]
  
  optimOut = optim(par = c(iota,rho,doInit), fn = to.optim)
  
  iota.hat    <- optimOut$par[1]*(1440/timestep) #iota and rho to units of mg O2 L-1 day-1 
  rho.hat     <- optimOut$par[2]*(1440/timestep) #iota and rho to units of mg O2 L-1 day-1 
  do.init.hat <- optimOut$par[3] #initial DO estimate 
  convergence <- optimOut$convergence  #did model converge or not (0=yes, 1=no)
  nll.val     <- optimOut$value #value of nll 
  gpp.hat     <- (iota.hat/(60*60*24))*sum(irr*timestep*60) #GPP in units of mg O2 L-1 d-1 
  
  
  return(list(iota=iota.hat, rho=rho.hat, gpp=gpp.hat,
              mod.info=list(do.init=do.init.hat, convergence=convergence, nll=nll.val)))
  
}

################################################################################
# Summary: Takes basic metabolism parameters and returns Predicted DO
#
# Authors: Gleon Student Fellows and Luke Winslow
# 
################################################################################
calc.do.hat = function(iota, rho, doInit, irr, doSat, zMix, k.gas){

  nObs = length(irr)
  #Set up output
  DOHat <- rep(NA,nObs)
  atmFlux <- rep(NA,nObs)  
  #Initialize DOHat
  DOHat[1] <- doInit
  
  #Calculate atmFlux and predicted DO for each time point
  #Fluxes out of lake have negative sign
  for (i in 1:(nObs-1)) {
    atmFlux[i] <- 	-k.gas[i] * (DOHat[i] - doSat[i]) / zMix[i]
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
calc.do.nll <- function(iota, rho, doInit, irr, doSat, zMix, kO2, doObs){
  
  modeled = calc.do.hat(iota, rho, doInit, irr, doSat, zMix, kO2)
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
  
  doHat  = calc.do.hat(iota, rho, doInit, irr, doSat, zMix, kO2)
  resids = doObs - doHat
  
  #If we are maintaining the ar1 component of the residuals, 
  # we must estimate ar1 coeff and the ar1 residual standard deviation
  if(ar1.resids){
    ar1.lm    = lm(resids[1:n.obs-1] ~ resids[2:n.obs]-1)
    ar1.coeff = ar1.lm$coefficients
    ar1.sd    = sd(ar1.lm$residuals)
  }
  
  toOptim <- function(iotaRhoDoinit){
    calc.do.nll(iotaRhoDoinit[1], iotaRhoDoinit[2], iotaRhoDoinit[3], irr, doSat, zMix, kO2, doSim)
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
      simRes = sample(resids, length(resids), replace=FALSE) # Should replace=TRUE? -RDB
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



