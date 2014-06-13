#rdb
addNAs <- function(x, ...){
	dL <- grepl("^[dD][oO][yY]$", names(x)) # matches doy, regardless of case
	yL <- grepl("^[yY]ear4?$", names(x))# matches Year, year, year4, Year4
	dateL <- grepl(".?date?.", names(x), ignore.case=TRUE) # matches anything with "date" in it, regardless of what else is or is not there
	if(any(dL)){
		names(x)[dL] <- "doy"
	}else{
		warning("No 'doy' column found")
	}
	if(any(yL)){
		names(x)[yL] <- "year"
	}else{
		warning("No 'year' column found")
	}
	if(any(dateL)){
		names(x)[dateL] <- "datetime"
	}else{
		warning("No 'date' column found")
	}
	if(!"POSIXct"%in%class(x[,"datetime"])){
		stop("date column must be POSIXct")
	}
	rdig <- 4
	Mode <- function(x){
			ux <- unique(x)
			ux[which.max(tabulate(match(x, ux)))]
	}
	ex <- round(Mode(1/diff(x[,"doy"])))
	mins <- 1/ex*24*60
	is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){abs(x - round(x)) < tol}
	if(!is.wholenumber(mins)){warning("Time between samples not whole number")}
	x1 <- byeShort(X=x, Expected=ex, ToCount="doy", TruncToCount=TRUE, ...)
	if(nrow(x1)==0){
		return(x1)
	}
	Range <- range(x1[,"datetime"]) #intentionally not truncating (having already used byeShort, I don't have to worry about starting and stopping at a specific time each day)
	# Ideal <- data.frame("RoundDoY"=round(seq(Range[1], Range[2], by=(1/ex)),rdig))
	Ideal <- data.frame("datetime"=seq(Range[1], Range[2], by=paste(mins, "mins")))
	# x1[,"RoundDoY"] <- round(x1[,"doy"], rdig)
	
	print(paste("NA's added to fill in time series:",dim(Ideal)[1]-dim(x1)[1], sep=" "))
	flush.console()
	x2 <- merge(x1, Ideal, all.x=TRUE, all.y=TRUE)
	if(any(yL)){x2[,"year"] <- approx(x2[,"datetime"], x2[,"year"], xout=Ideal[,1])$y}
	if(any(dL)){x2[,"doy"] <- approx(x2[,"datetime"], x2[,"doy"], xout=Ideal[,1])$y}
	
	return(x2)
}
