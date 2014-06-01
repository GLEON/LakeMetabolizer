
scale.exp.wind <- function(ts.data, wndZ){
	wnd = get.vars(tb.data$data, 'wnd')
	
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
	U10 <- wnd * (10/wndZ)^(1/7) 
	return(U10)
}