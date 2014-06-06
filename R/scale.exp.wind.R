
scale.exp.wind <- function(ts.data, wndZ){
	wnd = get.vars(ts.data, 'wnd')
	
	if(missing(wndZ)){
		wndZ = get.offsets(wnd)
    if(is.na(wndZ)){
			stop('Unknown wind height. Must supply wndZ parameter or have offset defined in header of ts.data.')
		}
	}
	
	if(ncol(wnd) > 2 || length(wndZ) > 1){
		stop('too many columns supplied to scale.exp.wind. Please supply only one datetime and wnd columns.')
	}
	
	u10 = scale.exp.wind.base(wnd[,2], wndZ)
	
	return(data.frame(datetime=ts.data$datetime, wnd_10=u10))
}

#Reference for this
#Arya 1988 (Introduction to micrometeorology)
scale.exp.wind.base <- function(wnd, wndZ){
	U10 <- wnd * (10/wndZ)^(0.15) # I'm pretty sure that (1/7) should be 0.15.  (1/7) used to be 1.7, but I think it was incorrectly copied from the U10^1.7 of k.cole.R. I have other code that has this value as 0.15.
	return(U10)
}