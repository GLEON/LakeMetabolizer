
#Reference for this
#Arya 1988 (Introduction to micrometeorology)
scale.exp.wind <- function(wnd, wndHeight){
	U10 <- wnd * (10/wndHeight)^(1/7) # should this be ^(1.7)?? I've used ^0.15 ... Changing for now. -RDB
	return(U10)
}