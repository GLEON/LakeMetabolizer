# ryanData <- function(max.rows=3E4, prefDay=NULL){
	origList <- ls()
	max.rows=3E4
	lake <- "troutbog"
	direc <- "/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/inst/extdata/"
	fileName <- paste(direc, paste(lake, c("wnd","par","wtr","doobs"), sep="."), sep="")
	
	# pickDay <- function(data, lake, prefDay){
	# 	myDay <- function(x,y){
	# 		x[,"DateTime"]>=prefDay & x[,"DateTime"]<(prefDay+1)
	# 	}
	# 	if(is.null(prefDay)){
	# 		stopifnot(lake%in%c("troutbog","mendota","sparkling"))
	# 		switch(lake,
	# 			troutbog=myDay(data,137),
	# 			sparkling=myDay(data,137),
	# 			mendota=myDay(data,137)
	# 		)	
	# 	}else{
	# 		myDay(data,prefDay)
	# 	}
	# }
	
	# wind data >> scale.exp.wind() >> k.cole() >> k600.2.kGAS() >>
	wo <- read.table(fileName[1], sep="\t", header=TRUE, colClasses=c("POSIXct","numeric"), nrows=max.rows)
	po <- read.table(fileName[2], sep="\t", header=TRUE, colClasses=c("POSIXct","numeric"), nrows=max.rows)
	to <- read.table(fileName[3], sep="\t", header=TRUE, nrows=max.rows)[,1:2]
	to[,1] <- as.POSIXct(to[,1])
	to[,2] <- as.numeric(to[,2])
	doo <- read.table(fileName[4], sep="\t", header=TRUE, colClasses=c("POSIXct","numeric"), nrows=max.rows)

	source("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/inProgress/longestRun.R")
	# ===========================
	# = Combine & Organize Data =
	# ===========================
	#merge
	d1 <- merge(to, doo, all=TRUE)
	d2 <- merge(wo, po, all=TRUE)
	d3 <- merge(d1,d2, all=TRUE)

	#convert to DoY format (not actually needed in this case, but conventient to have)
	d3[,1] <- LakeMetabolizer:::date2doy(d3[,1])
	# d3[,1] <- date2doy(d3[,1])

	#subset to the portion of the data set with the most consecutive observations (not necessary)
	# data0 <- d3[longestRun(d3),]
	# data0 <- d3[longestRun(d3[d3[,"DateTime"]>130&d3[,"DateTime"]<240,]),]
	# data0 <- d3[d3[,"DateTime"]>125&d3[,"DateTime"]<240,]
	# data0 <- d3[d3[,"DateTime"]>=136&d3[,"DateTime"]<137,]
	data0 <- d3[d3[,"DateTime"]>=137&d3[,"DateTime"]<138,]
	# data0 <- d3[d3[,"DateTime"]>=319&d3[,"DateTime"]<320,]
	# plot(data0[,c(1,3)], type="l")
	# plot(data0[1500:2500,c(1,3)], type="l")
	names(data0) <- c("DoY", "Temp", "DO", "Wind", "PAR") # rename columns while still data frame
	data0 <- as.matrix(data0) # convert to matrix
	row.names(data0) <- NULL # remove row names left over from d3

	data0[44, "Temp"] <- 10.635

	Freq <- median(diff(data0[,"DoY"])) # determine the sampling frequency; i have a function for mode if we are worried about it

	wind <- LakeMetabolizer:::scale.exp.wind(data0[,"Wind"], 2) # convert wind
	Kvec <- k600.2.kGAS(LakeMetabolizer:::k.cole(wind)*Freq, data0[,"Temp"], "O2") # calculate K for relevant sampling frequency


	data <- matrix(c(data0[,"DO"], LakeMetabolizer:::o2.at.sat(data0[,"Temp"], baro=716), Kvec, rep(1, dim(data0)[1]), data0[,"PAR"], data0[,"Temp"]), nrow=dim(data0)[1], dimnames=list(NULL, c("do.obs", "do.sat", "K", "Zmix", "irr", "wtr")))

	# return(data)
	rm(list=ls()[!ls()%in%c(origList,"data")])
# }

	


