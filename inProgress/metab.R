# wrapper function
# supply a data set, specify a method by name, and it returns daily metabolism
# will update to accept further arguments to be passed to fitting functions

install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", repos=NULL, type="source")
library("LakeMetabolizer")
# source("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/inProgress/ryanData.R")
# library("LakeMetabolizer")

metab <- function(data, method){
	require("plyr")
	
	data0 <- LakeMetabolizer:::ryanData()
	
	# Removes days with many missing ROWS:
	# data <- LakeMetabolizer:::addNAs(data0)
	
	# Removes days with many NA's:
	data <- LakeMetabolizer:::addNAs(data0[complete.cases(data0),])
	
	ids <- id(list(data[,"year"],trunc(data[,"doy"]))) # ID chunks to be analyzed
	
	
	
	
	
	
	
	
	
}
