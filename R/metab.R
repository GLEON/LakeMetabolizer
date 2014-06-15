
#' Calculate Metabolism
#'
#' Returns daily time series of gross primary production (GPP), respiration (R), and net ecosystem production (NEP). Depending on the method used, other information may be returned as well. Calculations are made using one of 5 statistical methods.
#'
#'
#'@param data a data.frame whose columns are "year", "doy", "datetime", "do.obs","do.sat","k.gas","z.mix", "irr", "wtr", "priors". Data columns (i.e., not year, doy, or datetime) that are not used by a particular statistical method do not need to be supplied.
#'@param method a character string specifying one of the 5 statistical methods ("bayesian", "bookkeep", "kalman", "ols", "mle")

metab <- function(data, method, wtr.name="wtr", ...){
	
	m.args <- list(...)
	
	if(wtr.name != "wtr"){
		names(data)[names(data)==wtr.name] <- "wtr"
	}
	
	# ===================
	# = Identify method =
	# ===================
	possibleMethods <- c("bayesian", "bookkeep", "kalman", "ols", "mle")
	mtd <- possibleMethods[which.min(adist(method, possibleMethods, ignore.case=TRUE))]
	if(!method%in%possibleMethods){
		warning(paste("method '",method,"' matched to '",mtd,"'. Supply perfect match to avoid warning.", sep=""))
	}
	stopifnot(length(mtd)==1)
	mtdCall <- paste("metab",mtd, sep=".")
	
	# ==============
	# = Groom data =
	# ==============
	# Removes days with many NA's:
	data1 <- addNAs(data[complete.cases(data),], percentReqd=1) # note that addNAs ALSO checks for POSIXct datetime, and adds year/doy
	data2 <- data1[complete.cases(data1),]
	
	# ==================================
	# = Prepare to apply metab to data =
	# ==================================
	ids <- id(list(data2[,"year"],trunc(data2[,"doy"]))) # ID chunks to be analyzed
	ids <- as.integer(ids - (min(ids)-1))
	nid <- length(unique(ids))
	results <- vector("list", nid)
	
	# ==================================
	# = Apply metab to subsets of data =
	# ==================================
	for(i in unique(ids)){
		
		poss.args <- c("do.obs","do.sat","k.gas","z.mix", "irr", "wtr", "datetime") # data2 columns that could correspond to arguments
		used.args <- poss.args[poss.args%in%names(data2)] # assuming arguments are used if they are in data2
		largs0 <- as.list(data2[i==ids, used.args]) # subsetting columns of data2 to used.args (just a safety check, probably not needed)
		largs <- c(largs0, m.args[!names(m.args)%in%names(largs0)]) # adding on any other arguments supplied via ...
		# note that in largs, argument supplied through data/data2/poss.args take precedent over arguments from ...
		
		# print(paste("Analyzing day #", i)); flush.console(); # Is this annoying? I'm commenting-out
		results[[i]] <- do.call(mtdCall, largs) # this is where all of the work happens
	}
	answer0 <- conquerList(results, naming=data.frame("year"=data2[!duplicated(ids),"year"], "doy"=trunc(data2[!duplicated(ids),"doy"])))
	
	
	a0.names <- names(results[[1]])
	
	# =======================================================
	# = Add non-metab list elements as attributes to output =
	# =======================================================
	# only need to add attributes if it's a list (more than 1 element, not a data frame)
	if(length(a0.names)>1 & is.list(answer0) & !is.data.frame(answer0)){
		
		names(answer0) <- a0.names
		answer <- answer0$metab
		for(i in 1:length(a0.names)){
			if(a0.names[i]=="metab"){next}
			if(a0.names[i]=="smoothDO"){ # do a little extra attribute work if smoothDO
				t.sDO <- answer0[[a0.names[i]]] # grab the element of the list that contains smoothed DO concs
				t.sDO <- t.sDO[,!names(t.sDO)%in%c("doy","year")] # remove the columns that are the year/ doy
				attr(answer, "smoothDO.vec") <- c(t(t.sDO)) # provide the smoothed DO as a long vector, instead of each row containing nobs+2 columns of smoothed DO from a given day
			}
			attr(answer, a0.names[i]) <- answer0[[a0.names[i]]] # assign non "metab" list element as attribute
		}
		
	}else{ # if the list only has one element, or if the list is really a data.frame, then answer0 is what we want
		answer <- answer0
	}
	
	return(answer)
}


