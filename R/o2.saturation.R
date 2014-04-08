###########################################################################
# Summary: Estimates equilibrium dissolved oxygen (DO) concentration
# 
# Input: Water temperature in deg Celsius, 
#        Barometric Pressure in millibar
#         OR
#        Altitude in meters above sea level (Default:0 m)
#        Model name to use for o2 sat calc ['weiss', 'benson']
#
# Output: Equilibrium DO concentration in mg/L in given conditions
#
# Reference: Benson & Krause (1984)
#            Weiss 1970
#            Garcia and Gordon 1992
#            USGS memo #81.11 1981
#            USGS memo #81.15 1981
#
# Initial Author: Luke Winslow
###########################################################################

o2.at.sat = function(x, ...) UseMethod("o2.at.sat")

o2.at.sat.data.frame = function(temp, baro, altitude=0, salinity=0, model='garcia'){
  if(ncol(temp) > 2){
    stop('Temp can only have two columns, "datetime" and temperature')
  }
  
  dosat = o2.at.sat(temp[,2], baro, altitude, salinity, model)
  
  return(data.frame(datetime=temp$datetime, dosat=dosat)) 
  
}


o2.at.sat.default <- function(temp, baro, altitude=0, salinity=rep(0,length(temp)), model='garcia'){
  
  
  if(!missing(baro)){#Calc using barometric pressure
    press.corr = (baro * 0.0987 - 0.0112)/100
  }else{
    press.corr = (0.0000005 * altitude^2 - 0.0118 * altitude + 99.979)/100
  }
  
  
  if(tolower(model) == 'garcia'){
  #Garcia, HE, and LI Gordon. 1992. “Oxygen Solubility in Seawater: Better Fitting Equations.” 
  #Limnology and Oceanography 37 (6). 
    
    Ts = log((298.15 - temp)/(273.15 + temp))
    
    lnC = 2.00856 + 3.22400 *Ts + 3.99063*Ts^2 + 4.80299*Ts^2 + 4.80299*Ts^3 + 9.78188e-1*Ts^4 + 
      1.71069*Ts^5 - salinity*(6.24097e-3 + 6.93498e-3*Ts + 6.90358e-3*Ts^2 + 4.29155e-3*Ts^3) - 3.1168e-7*salinity^2
    
    o2.sat = exp(lnC)  * 1.423 #convert from ml/l to mg/l
    
  }else if(tolower(model) == 'weiss'){
    tempk = temp + 273.15
    
    lnC = -173.4292 + 249.6339 * (100 / tempk) + 143.3483 *
          log(tempk / 100) - 21.8492 * (tempk / 100) + 
          salinity * (-0.033096 + 0.014259 * (tempk / 100) - 0.0017000 * (tempk / 100)^2)
                        
    o2.sat = exp(lnC) * 1.423 #convert from ml/l to mg/l
  
  }else if(tolower(model) == 'benson'){
    ## TODO: Fix this to include salinity
    if(!all(salinity==0)){
      warning('Benson model does not currently include salinity')
    }
    
    o2.sat = (-0.00006 * (temp)^3) + (0.00725 * (temp)^2) - (0.39571 * (temp)) + 14.59030
    
  }

  o2.sat = o2.sat * press.corr
  
  return(o2.sat)
}
