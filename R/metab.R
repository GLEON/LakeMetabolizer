
#' Calculate Metabolism
#'
#' Returns daily time series of gross primary production (GPP), respiration (R), and net ecosystem production (NEP). Depending on the method used, other information may be returned as well. Calculations are made using one of 5 statistical methods.
#'
#'
#'@param data a data.frame whose columns are "year", "doy", "datetime", "do.obs","do.sat","k.gas","z.mix", "irr", "wtr", "priors". Data columns (i.e., not year, doy, or datetime) that are not used by a particular statistical method do not need to be supplied.
#'@param method a character string specifying one of the 5 statistical methods ("bayesian", "bookkeep", "kalman", "ols", "mle")



metab <- function(data, method){
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
	# data0 <- LakeMetabolizer:::ryanData()
	
	# Removes days with many missing ROWS:
	# data <- LakeMetabolizer:::addNAs(data0)
	
	# Removes days with many NA's:
	# data1 <- LakeMetabolizer:::addNAs(data0[complete.cases(data0),], percentReqd=1)
	data1 <- addNAs(data[complete.cases(data),], percentReqd=1)
	data2 <- data1[complete.cases(data1),]
	
	# ==================================
	# = Prepare to apply metab to data =
	# ==================================
	ids <- id(list(data2[,"year"],trunc(data2[,"doy"]))) # ID chunks to be analyzed
	ids <- as.integer(ids - (min(ids)-1))
	nid <- length(unique(ids))#attributes(ids)$n
	# metabArgs <- as.list(data[1:10,c("do.obs","do.sat","k.gas","z.mix", "irr", "wtr")])
	results <- vector("list", nid)
	
	# ==================================
	# = Apply metab to subsets of data =
	# ==================================
	for(i in unique(ids)){
		
		poss.args <- c("do.obs","do.sat","k.gas","z.mix", "irr", "wtr", "priors") # column names that could correspond to arguments in a metab.xx() function
		used.args <- poss.args[poss.args%in%names(data2)] # assuming that arguments are used if they are found in column names
		# note that all metab.xx() have the ... argument, so it is OK to supply extra arguments. However, it is important to keep the order in poss.args and in the metab functions as-is
		largs <- as.list(data2[i==ids, used.args])
		
		results[[i]] <- do.call(mtdCall, largs)
	}
	answer <- conquerList(results, naming=data.frame("year"=data2[!duplicated(ids),"year"], "doy"=trunc(data2[!duplicated(ids),"doy"])))
	return(answer)
}

# ===========
# = EXAMPLE =
# ===========
# metab(data=tb.data, method="kalman")
# metab(data=tb.data, method="bayes")


