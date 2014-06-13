
#' Calculate Metabolism
#'
#' Returns daily time series of gross primary production (GPP), respiration (R), and net ecosystem production (NEP). Depending on the method used, other information may be returned as well. Calculations are made using one of 5 statistical methods.
#'
#'
#'@param data a data.frame whose columns are "year", "doy", "datetime", "do.obs","do.sat","k.gas","z.mix", "irr", "wtr", "priors". Data columns (i.e., not year, doy, or datetime) that are not used by a particular statistical method do not need to be supplied.
#'@param method a character string specifying one of the 5 statistical methods ("bayesian", "bookkeep", "kalman", "ols", "mle")



metab <- function(data, method, ...){
	
	m.args <- list(...)
	
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
	
	print(paste("data names:", paste(names(data), collapse=" "))); flush.console();
	print(paste("data1 names:", paste(names(data1), collapse=" "))); flush.console();
	print(paste("data2 names:", paste(names(data2), collapse=" "))); flush.console();
	
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
		
		poss.args <- c("do.obs","do.sat","k.gas","z.mix", "irr", "wtr", "datetime") # data2 columns that could correspond to arguments
		used.args <- poss.args[poss.args%in%names(data2)] # assuming arguments are used if they are in data2
		largs0 <- as.list(data2[i==ids, used.args]) # subsetting columns of data2 to used.args (just a safety check, probably not needed)
		print(paste("largs0 names:", paste(names(largs0), collapse=" "))); flush.console();
		largs <- c(largs0, m.args[!names(m.args)%in%names(largs0)]) # adding on any other arguments supplied via ...
		print(paste("largs names:", paste(names(largs), collapse=" "))); flush.console();
		# note that in largs, argument supplied through data/data2/poss.args take precedent over arguments from ...
		
		print(paste("Analyzing day #", i)); flush.console();
		results[[i]] <- do.call(mtdCall, largs)
	}
	answer0 <- conquerList(results, naming=data.frame("year"=data2[!duplicated(ids),"year"], "doy"=trunc(data2[!duplicated(ids),"doy"])))
	
	
	a0.names <- names(results[[1]])
	
	if(length(a0.names)>1 & is.list(answer0) & !is.data.frame(answer0)){
		
		names(answer0) <- a0.names
		answer <- answer0$metab
		for(i in 1:length(a0.names)){
			if(a0.names[i]=="metab"){next}
			if(a0.names[i]=="smoothDO"){
				attr(answer, "smoothDO.vec") <- answer0[[a0.names[i]]]
			}
			attr(answer, a0.names[i]) <- answer0[[a0.names[i]]]
		}
		
	}else{
		answer <- answer0
	}
	
	return(answer)
}

# ===========
# = EXAMPLE =
# ===========
# metab(data=tb.data, method="kalman")
# metab(data=tb.data, method="bayes")


