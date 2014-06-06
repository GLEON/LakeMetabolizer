
metab.bookkeep <- function(do.obs, do.sat, k.gas, z.mix, date.times, lake.lat, ...){
  #do.obs     - Concentration units
  #do.sat     - concentration units
  #k.gas      - piston velocity (m/day)
  #z.mix      - depth in meters
  #date.times - in POSIXct data structure
  
  
  delta.do <- diff(do.obs)
  mid.times <- diff(date.times)/2 + date.times[1:(length(date.times)-1)]
  delta.times <- diff(date.times)
  
  k.gas.timestep <- k.gas[1:(length(k.gas)-1)] * as.numeric(delta.times, 'secs')/(60*60*24) #convert to per-timestep
  
  dayI <- is.day(lake.lat, mid.times)
  nightI <- !dayI
  
  #gas flux out is negative
  #normalized to z.mix, del_concentration/timestep (e.g., mg/L/10min)
  gas.flux <- (do.sat - do.obs) * k.gas / z.mix
  
  #remove the component of delta.do that is due to gas flux
  delta.do.metab <- delta.do + gas.flux[1:(length(gas.flux)-1)]
  
  #normalize units to per-day
  delta.do.meta.daily <- delta.do.metab * (60*60*24)/as.numeric(delta.times, 'secs')
  
  R <- mean(delta.do.meta.daily[nightI]) #this should be negative
  NEP <- mean(delta.do.meta.daily)       #pos or negative
  GPP <- NEP - R                    #should be positive
  
  metab <- c("GPP"=GPP, "R"=R)
  return(metab)
  
}

