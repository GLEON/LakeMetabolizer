#'@title Calculates the equilibrium saturation concentration of oxygen in water 
#'  at the supplied conditions
#'@name o2.at.sat
#'@aliases o2.at.sat o2.at.sat.base
#'@usage o2.at.sat.base(temp, baro, altitude = 0, salinity = rep(0, 
#'  length(temp)), model = "garcia") o2.at.sat(ts.data, baro, altitude = 0, 
#'  salinity = 0, model = "garcia")
#'@description Used to calculate the equilibrium concentration of oxygen in 
#'  water. The equilibration concentration of oxygen in water varies with both 
#'  temperature, salinity, and the partial pressure of oxygen in contact with 
#'  the water (calculated from supplied elevation or barometric pressure).
#'@param ts.data Object of class data.frame with two named columns 
#'  \dQuote{datetime} and \dQuote{wtr} (water temp in deg C).
#'@param temp a numeric vector of water temperature in degrees Celsius.
#'@param baro barometric pressure in millibars.
#'@param altitude a numeric value indicating the elevation above mean sea level 
#'  in meters. Defaults to mean sea level. An alternative to supplying 
#'  barometric pressure.
#'@param salinity a numeric vector of salinity in PSU. Defaults to zero. Length 
#'  must be one or equal to length of temperature.
#'@param model the empirical model to be used, \code{"garcia"}, 
#'  \code{"garcia-benson"}, \code{"weiss"} and \code{"benson"} are the available
#'  options. They correspond to the references described below, where both 
#'  \code{"garcia"} and \code{"garcia-benson"} are from the first reference.
#'  
#'@return The equilibration concentration at the supplied conditions in mg/L of 
#'  oxygen.
#'@author Luke A Winslow
#'@references
#'
#'Garcia, H., and L. Gordon (1992), \emph{Oxygen solubility in seawater: Better 
#'fitting equations}, Limnol. Oceanogr., 37(6).
#'
#'Benson, B. B., & Krause, D. (1984). \emph{The concentration and isotopic 
#'fractionation of oxygen dissolved in freshwater and seawater in equilibrium 
#'with the atmosphere}. Limnology and Oceanography, 29(3), 620-632. 
#'doi:10.4319/lo.1984.29.3.0620
#'
#'Weiss, R. (1970). \emph{The solubility of nitrogen, oxygen and argon in water 
#'and seawater}. Deep Sea Research and Oceanographic Abstracts, 17(4), 721-735. 
#'doi:10.1016/0011-7471(70)90037-9
#'
#'@seealso \link{water.density}, \link{o2.at.sat.base}
#'@keywords math, methods
#'@examples
#'temp.range = 1:25
#'sal.range = 1:25
#'
#'par(mfrow=c(1,2))
#'plot(temp.range, o2.at.sat.base(temp.range), xlab='Temperature (C)', 
#'ylab='Oxygen Saturation (mg/L)')
#'plot(o2.at.sat.base(rep(20,25), salinity=sal.range), xlab='Salinity (PSU)', ylab='')
#'
#'@export
o2.at.sat <- function(ts.data, baro, altitude=0, salinity=0, model='garcia'){
	if(ncol(ts.data) > 2){
		stop('Temp can only have two columns, "datetime" and "wtr"')
	}

	dosat <- o2.at.sat.base(ts.data$wtr, baro, altitude, salinity, model)

	return(data.frame(datetime=ts.data$datetime, do.sat=dosat))  
}

#'@export
o2.at.sat.base <- function(temp, baro, altitude=0, salinity=rep(0,length(temp)), model='garcia'){
  

  # Correction for vapor pressure of water
	if(!missing(baro)){#Calc using barometric pressure
		press.corr <- (baro * 0.0987 - 0.0112)/100 # USGS memo #81.15 1981 as cited in Staehr et al. 2010
	} else {
		press.corr <- (0.0000005 * altitude^2 - 0.0118 * altitude + 99.979)/100 # USGS memo #81.15 1981 as cited in Staehr et al. 2010
	}
  
	if(tolower(model) == 'garcia'){

	  Ts <- log((298.15 - temp)/(273.15 + temp))

	  lnC <- 2.00856 + 3.22400 *Ts + 3.99063*Ts^2 + 4.80299*Ts^3 + 9.78188e-1*Ts^4 + 
	    1.71069*Ts^5 - salinity*(6.24097e-3 + 6.93498e-3*Ts + 6.90358e-3*Ts^2 + 4.29155e-3*Ts^3) - 3.1168e-7*salinity^2

	  o2.sat <- exp(lnC)  * 1.42905 #convert from ml/l to mg/l

	} else if(tolower(model) == 'garcia-benson'){
	  
	  Ts <- log((298.15 - temp)/(273.15 + temp))
	  
	  lnC <- 2.00907 + 3.22014*Ts + 4.05010*Ts^2 + 4.94457*Ts^3 + -2.56847e-1*Ts^4 + 
	    3.88767*Ts^5 - salinity*(6.24523e-3 + 7.37614e-3*Ts + 1.03410e-2*Ts^2 + 8.17083e-3*Ts^3) - 4.88682e-7*salinity^2
	  
	  o2.sat <- exp(lnC)  * 1.42905 #convert from ml/l to mg/l
	  
	} else if(tolower(model) == 'weiss'){
		tempk <- temp + 273.15

		lnC <- -173.4292 + 249.6339 * (100 / tempk) + 143.3483 *
		log(tempk / 100) - 21.8492 * (tempk / 100) + 
		salinity * (-0.033096 + 0.014259 * (tempk / 100) - 0.0017000 * (tempk / 100)^2)
          
		o2.sat <- exp(lnC) * 1.423 #convert from ml/l to mg/l

	}else if(tolower(model) == 'benson'){
		## TODO: Fix this to include salinity
		if(!all(salinity==0)){
			warning('Benson model does not currently include salinity')
		}

	o2.sat <- (-0.00006 * (temp)^3) + (0.00725 * (temp)^2) - (0.39571 * (temp)) + 14.59030

	}

	o2.sat <- o2.sat * press.corr

	return(o2.sat)
}
