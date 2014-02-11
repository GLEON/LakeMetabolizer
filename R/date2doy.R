date2doy <- function(x){
	stopifnot(any(grepl("^POSIX",class(wo[,1]))))
	day <- as.numeric(format(x, "%j"))
	pat <- quote("([0-9]{2}:){2}[0-9]{2}")
	midnight <- gsub(pat, "00:00:00", x)
	frac <- as.numeric(difftime(x, midnight, units="days"))
	day+frac
}
