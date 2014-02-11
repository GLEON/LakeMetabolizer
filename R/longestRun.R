# RDB
# Find the longest streak of consecutive observations (non-NA, non-NaN) in a data.frame, matrix, or vector, and return the indices corresponding to the rows (or elements) of that longRun
longestRun <- function(x){
	runs <- rle(complete.cases(x))  # the length of consecutive trues / false
	rL <- runs$lengths
	rV <- runs$values
	maxR <- max(rL[rV])
	wMax <- which(rL==maxR)
	if(length(wMax)!=1){
		ties <- c("(starting at:")
		for(i in 1:length(wMax)){
			thisRun <- sum(rL[1:(wMax[i]-1)])+1
			ties <- paste(ties, " ", thisRun, sep="")
		}
		ties <- paste(ties, ")", sep="")
		warning(paste("multiple runs tied for longest", ties, " -- choosing first"))
		wMax <- wMax[1]
	}
	startMaxR <- sum(rL[1:(wMax-1)])+1
	endMaxR <- startMaxR+rL[wMax]-1
	startMaxR:endMaxR
}
