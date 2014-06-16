
metab.bookkeep <- function(do.obs, do.sat, k.gas, z.mix, irr, ...){
	#do.obs     - Concentration units
	#do.sat     - concentration units
	#k.gas      - piston velocity (m/day)
	#z.mix      - depth in meters
	#datetimes - in POSIXct data structure

	nobs <- length(do.obs)  

	mb.args <- list(...)
	if("datetime"%in%names(mb.args)){ # check to see if datetime is in the ... args
		datetime <- mb.args$datetime # extract datetime
		freq <- calc.freq(datetime) # calculate sampling frequency from datetime
		if(nobs!=freq){ # nobs and freq should agree, if they don't issue a warning
			bad.date <- format.Date(datetime[1], format="%Y-%m-%d")
			warning("number of observations on ", bad.date, " (", nobs, ") ", "does not equal estimated sampling frequency", " (", freq, ")", sep="")
		}
	}else{ # if datetime is *not* in the ... args
		warning("datetime not found, inferring sampling frequency from # of observations") # issue a warning (note checks in addNAs)
		# NOTE: because of the checks in addNA's, it is unlikely a user would receive this warning via metab()
		# warning will only be seen through direct use of metab.bookkeep when datettime is not supplied
		freq <- nobs
	}


	if(all(c("datetime", "lake.lat")%in%names(mb.args))){ # check to see if datetime and lake.late are in mb.args
		irr <- as.integer(is.day(datetimes=datetime, lat=mb.args$lake.lat)) # calculate 1's and 0's for irr from datetime and lake.lat
		dayI <- irr == 1L # TRUE when the lights are on
		nightI <- irr == 0L	# TRUE when the lights are off
	}else{
		if(!all(irr==1L | irr==0L)){ # datetime/lake.lat are not found, and if all elements of irr are not integer 0/1
			stop("either supply datetime & lake.lat arguments, or supply irr as integer vector of 1's and 0's") # then that's an error
		}
		# but if datetime/lake.lat are not found in mb.args, yet irr is all integer 0/1, then ...
		dayI <- irr == 1L 
		nightI <- irr == 0L
	}

	delta.do <- diff(do.obs)
	miss.delta <- sum(is.na(delta.do)) # number of NA's
	if(miss.delta != 0){ # the number of NA's should be 0, otherwise issue a warning
		warning(paste(miss.delta, " missing values (", miss.delta/length(delta.do), "%) in diff(do.obs)", sep=""))
		# Note: it should be hard to get this warning. Should be impossible to get this warning through metab due to addNA's
	}

	#gas flux out is negative
	#normalized to z.mix, del_concentration/timestep (e.g., mg/L/10min)
	gas.flux <- (do.sat - do.obs) * (k.gas/freq) / z.mix 

	#remove the component of delta.do that is due to gas flux
	delta.do.metab <- delta.do + gas.flux[1:(length(gas.flux)-1)]

	#normalize units to per-day
	# delta.do.meta.daily <- delta.do.metab * (60*60*24)/as.numeric(delta.times, 'secs')

	nep.day <- delta.do.metab[dayI]
	nep.night <- delta.do.metab[nightI]


	R <- mean(nep.night, na.rm=TRUE) * freq # should be negative
	NEP <- mean(delta.do.metab, na.rm=TRUE) * freq # can be positive or negative
	GPP <- mean(nep.day, na.rm=TRUE) * sum(dayI) - R # should be positive

	metab <- data.frame("GPP"=GPP, "R"=R, "NEP"=NEP)
	return(metab)

}

