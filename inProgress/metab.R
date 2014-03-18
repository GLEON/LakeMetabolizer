# wrapper function
# supply a data set, specify a method by name, and it returns daily metabolism
# will update to accept further arguments to be passed to fitting functions

# install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", repos=NULL, type="source")
# library("LakeMetabolizer")
# source("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/inProgress/ryanData.R")
# library("LakeMetabolizer")

metab <- function(data, method){
	require("plyr")
	
	data0 <- LakeMetabolizer:::ryanData()
	data0[,"Year"] <- as.integer(format.Date(data0[,"date"], "%Y"))
	data <- addNA(data0)
	
	
}
