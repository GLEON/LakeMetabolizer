
metab.bookkeep <- function(do.obs, do.sat, k.gas, z.mix, irr, ...){
  #do.obs     - Concentration units
  #do.sat     - concentration units
  #k.gas      - piston velocity (m/day)
  #z.mix      - depth in meters
  #date.times - in POSIXct data structure
  
  mb.args <- list(...)

  if(c("datetime", "lake.lat")%in%names(mb.args)){
    irr <- as.integer(is.day(lake.lat, datettime))
    dayI <- irr == 1L
    nightI <- irr == 0L	
  }else{
	dayI <- irr == 1L
	nightI <- irr == 0L
  }
  
  delta.do <- diff(do.obs)
  miss.delta <- sum(is.na(delta.do))
  if(miss.delta != 0){
	warning(paste(miss.delta, " missing values (", miss.delta/length(delta.do), "%) in diff(do.obs)", sep=""))
  }
  
  #gas flux out is negative
  #normalized to z.mix, del_concentration/timestep (e.g., mg/L/10min)
  gas.flux <- (do.sat - do.obs) * k.gas / z.mix
  
  #remove the component of delta.do that is due to gas flux
  delta.do.metab <- delta.do + gas.flux[1:(length(gas.flux)-1)]
  
  #normalize units to per-day
  # delta.do.meta.daily <- delta.do.metab * (60*60*24)/as.numeric(delta.times, 'secs')

  nep.day <- delta.do.metab[dayI]
  nep.night <- delta.do.metab[nightI]
  
  nobs <- length(do.obs)
  R <- mean(nep.night, na.rm=TRUE) * nobs # should be negative
  NEP <- mean(delta.do.metab, na.rm=TRUE) * nobs # can be positive or negative
  GPP <- mean(nep.day, na.rm=TRUE) * sum(dayI) - R # should be positive
  
  # metab <- matrix(c(GPP, R, NEP), nrow=1, dimnames=list(NULL, c("GPP", "R", "NEP")))
  metab <- data.frame("GPP"=GPP, "R"=R, "NEP"=NEP)
  return(metab)
  
}

