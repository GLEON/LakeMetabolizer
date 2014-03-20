
metabData <- function(do.ts, wtrAll.ts, wnd.ts, irr.ts, atmPress=716/0.75006168, wndHeight){
	require("rLakeAnalyzer")

	Mode <- function(x){
			ux <- unique(x)
			ux[which.max(tabulate(match(x, ux)))]
	}
	do.doy <- date2doy(do.ts[,1])
	Freq <- round(Mode(1/diff(do.doy)))
	
	# do.obs, do.sat, k.gas, z.mix, irr, wtr
	
	m1.ts <- merge(wnd.ts, wtrAll.ts[,1:2], all=TRUE)
	names(m1.ts) <- c("datetime", "wnd", "wtr")
	
	wind <- scale.exp.wind(m1.ts[,"wnd"], wndHeight) # convert wind
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


