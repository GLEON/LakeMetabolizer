# Given these "time series data frames", calculate K and prepare a data frame for metab(). The ".ts" indicates that the first column must be in a POSIXct date format, and the second+ columns contain data. The date column should have the header "datetime".

# 1) do.ts contains observed oxygen values in the 2nd column in units of mg/L. Header name is not referenced.
# 2) wtrAll.ts contains depth profile data. The 2nd column should be at the sonde depth. All columns (after 1st) should have format of wtr_x, where x is the depth in meters, such that a 10 m thermistor would have the header wtr_10.
# 3) wnd.ts contains wind speed in meters per second (2nd column). Header name is not referenced.
# 4) irr.ts contains PAR values. Current the unit are unimportant for metabolism values (but is probably important for some of the physics...). Header name (2nd column) should be "irr".
# 5) atmPress is a scale (or numeric vector) of atmospheric pressures in units of mBar. Conversion from mmHg to mbar is "mbar = mmHg/0.75006168".
# 6) wndHeight is the height of the anemometer in meters.

# NOTE: Output from this function is suitable to be directly supplied to metab()
# NOTE: This function currently uses k.cole, and does not accept an argument for other methods. This needs to be modified.

metabData <- function(do.ts, wtrAll.ts, wnd.ts, irr.ts, atmPress=716/0.75006168, wndHeight){

	Mode <- function(x){
			ux <- unique(x)
			ux[which.max(tabulate(match(x, ux)))]
	}
	do.doy <- date2doy(do.ts[,1])
	Freq <- round(Mode(1/diff(do.doy)))
	
	# do.obs, do.sat, k.gas, z.mix, irr, wtr
	
	m1.ts <- merge(wnd.ts, wtrAll.ts[,1:2], all=TRUE)
	names(m1.ts) <- c("datetime", "wnd", "wtr")
	
	wind <- wind.scale(m1.ts[,"wnd"], wndHeight) # convert wind
	Kvec <- k600.2.kGAS(k.cole(wind), m1.ts[,"wtr"], "O2")/Freq # calculate K for relevant sampling frequency
	m2.ts <- cbind(m1.ts, "k.gas"=Kvec)
	
	do.sat <- o2.at.sat(m2.ts[,"wtr"], baro=atmPress)
	m3.ts <- cbind(m2.ts, "do.sat"=do.sat)
	
	z.mix.ts <- ts.thermo.depth(wtrAll.ts)
	names(z.mix.ts) <- c("datetime", "z.mix")
	m4.ts <- merge(m3.ts, z.mix.ts, all=TRUE)
	
	names(do.ts) <- c("datetime", "do.obs")
	m5.ts <- merge(m4.ts, do.ts, all=TRUE)
	
	m6.ts <- merge(m5.ts, irr.ts, all=TRUE)

	data <- data.frame(
		"date"=m6.ts[,"datetime"],
		"year"= as.integer(format.Date(m6.ts[,"datetime"], "%Y")),
		"doy"=date2doy(m6.ts[,"datetime"]),
		"do.obs"=m6.ts[,"do.obs"], 
		"do.sat"=m6.ts[,"do.sat"],
		"k.gas"=m6.ts[,"k.gas"], 
		"z.mix"=m6.ts[,"z.mix"], 
		"irr"=m6.ts[,"irr"], 
		"wtr"=m6.ts[,"wtr"]
	)
	
	return(data)
}


