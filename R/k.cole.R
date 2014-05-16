# ---Author: Jake Zwart, 2013-09-30 --- 
# cleaned from previous script - original author: Arianto Santoso 

k.cole <- function(x, ...) UseMethod("k.cole")


k.cole.data.frame <- function(ts.data){
	if(!has.vars(ts.data, 'wnd')){
		stop('k.cole requires a "wnd" column in the supplied data')
	}
	wnd <- get.vars(ts.data, 'wnd')
	k600 <- k.cole(wnd[,2])

	return(data.frame(datetime=ts.data$datetime, k600=k600))
}


# wnd: wind value in m/s

# OUTPUT: returns the gas exchange velocity for k600 in units of m/day

k.cole.default <- function(wnd){
  U10 <- wnd  #This function uses just the wind speed it is supplied. 
  k600 <- 2.07 + (0.215 * (U10^1.7)) # units in cm h-1
  k600 <- k600*24/100 #units in m d-1
  return(k600)
}

# -- References 
#** COLE, J.J. and CARACO, N.F. 1998. Atmospheric exchange of **
#** carbon dioxide in a low-wind oligotrophic lake measured by **
#** the addition of SF6. Limnology and Oceanography. 43(4): 647-656. 

#** WANNINKHOF, R. 1992. Relationship between gas exchange and **
#** wind speed over the ocean. J. Geophys. Res. 97: 7373-7382.
