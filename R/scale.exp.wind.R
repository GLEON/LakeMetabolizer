
#Reference for this
#Arya 1988 (Introduction to micrometeorology)
scale.exp.wind <- function(wnd, wndHeight){
	U10 <- wnd * (10/wndHeight)^1/7
	return(U10)
}