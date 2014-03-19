# wrapper function
# supply a data set, specify a method by name, and it returns daily metabolism
# will update to accept further arguments to be passed to fitting functions

# detach(package:LakeMetabolizer)
# install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", repos=NULL, type="source")
# library("LakeMetabolizer")
# source("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/inProgress/ryanData.R")
# library("LakeMetabolizer")

metab <- function(data, method, fitBy=c("doy","year")){
	require("plyr")
	
	data0 <- LakeMetabolizer:::ryanData()
	
	data <- LakeMetabolizer:::addNAs(data0)
	# NOTE: The above removes days that have many missing ROWS; to remove days with many NA's, do:
	# data <- addNA(data0[complete.cases(data0),])
	
	
	
	
	
	
	
	
}
