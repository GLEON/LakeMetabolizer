#rdb
addNAs <- function(x){
	rdig <- 4
	Mode <- function(x){
			ux <- unique(x)
			ux[which.max(tabulate(match(x, ux)))]
	}
	ex <- round(Mode(1/diff(x[,"doy"])))
	x1 <- byeShort(X=x, Expected=ex, ToCount="doy", TruncToCount=TRUE)
	if(nrow(x1)==0){
		return(x1)
	}
	Range <- range(x1[,"doy"]) #intentionally not truncating (having already used byeShort, I don't have to worry about starting and stopping at a specific time each day)
	Ideal <- data.frame("RoundDoY"=round(seq(Range[1], Range[2], by=(1/ex)),rdig))
	x1[,"RoundDoY"] <- round(x1[,"doy"], rdig)
	print(paste("NA's added to fill in time series:",dim(Ideal)[1]-dim(x1)[1], sep=" "))
	flush.console()
	x2 <- merge(x1, Ideal, all.x=TRUE, all.y=TRUE)
	x3 <- x2[,names(x2)[names(x2)!="RoundDoY"]]
	return(x3)
}
