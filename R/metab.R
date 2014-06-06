# wrapper function
# supply a data set, specify a method by name, and it returns daily metabolism
# will update to accept further arguments to be passed to fitting functions

# install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", repos=NULL, type="source")
# library("LakeMetabolizer")
# method="kalman"
# source("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer/inProgress/ryanData.R")
# library("LakeMetabolizer")




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
	require("plyr")
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
	
	# ==================
	# = Conquer a List =
	# ==================
	conquerList <- function(x, naming=NULL){
		# If x is not a list, don't bother
		if(!is.list(x)){return(x)}
		
		s1 <- length(x)
		s2 <- length(x[[1]])
		u1 <- unlist(x, recursive=FALSE)
		stopifnot(length(u1)==s1*s2)
		
		# return value from ldply() if it will work (e.g., if each element of list x contains a row of a data frame)
		if(s2==1 & (!is.data.frame(x)|!is.list(x))){
			return(ldply(x))
		}
		
		#
		s2C <- unlist(lapply(x[[1]], class))
		cqd <- vector("list", s2)
		for(i in 1:s2){
			ti <- seq(i, s1*s2, s2)
			tl <- vector("list", s1)
			for(j in 1:s1){
				tl[[j]] <- u1[[ti[j]]]
			}
			if(is.data.frame(tl[[1]])|!is.list(tl[[1]])){
				if(!is.null(naming)){
					cqd[[i]] <- cbind(naming,ldply(tl))
				}else{
					cqd[[i]] <- ldply(tl)
				}
			}else{
				cqd[[i]] <- llply(tl)
			}
		}
		return(cqd)
	}
	# ==================================
	# = Apply metab to subsets of data =
	# ==================================
	for(i in unique(ids)){
		# do.obs, do.sat, k.gas, z.mix, irr, wtr
		
		bk.args <- do.obs, do.sat, k.gas, z.mix, date.times, lake.lat, ...
		ols.args <- do.obs, do.sat, k.gas, z.mix, irr, wtr, ...
		mle.args <- do.obs, do.sat, k.gas, z.mix, irr, wtr
		kal.args <- do.obs, do.sat, k.gas, z.mix, irr, wtr
		bayes.args <- do.obs, do.sat, k.gas, z.mix, irr, wtr
		
		
		largs <- as.list(data2[i==ids,c("do.obs","do.sat","k.gas","z.mix", "irr", "wtr")])
		
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


