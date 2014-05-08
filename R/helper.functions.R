#helper.functions 

has.vars = function(data, var.names){
  
  if(!is(data, 'data.frame')){
    stop('Data must be of class data.frame')
  }
  
  header = names(data)
  
  has.var = rep(FALSE, times=length(var.names))
  
  for(i in 1:length(var.names)){
    has.var[i] = any(grepl(paste(var.names[i], '[\\w_]?', sep=''), header, ignore.case=TRUE))
  }
  
  return(has.var)  
}


get.vars = function(data, var.names){
  
  if(!is(data, 'data.frame')){
    stop('Data must be of class data.frame')
  }
  
  header = names(data)
  
  datetimeI = grepl('datetime', header, ignore.case=TRUE)
  
  if(!any(datetimeI)){
    stop("Can't find datetime column in supplied data.")
  }
  
  varI = rep(FALSE, times=length(var.names))
  
  for(i in 1:length(var.names)){
    varI = varI | grepl(paste(var.names[i], '[\\w_]?', sep=''), header, ignore.case=TRUE)
  }
  
  if(!any(varI)){
    stop("No variable pattern matches:", paste(var.names, collapse=' '))
  }
  
  return(data[, varI | datetimeI])
}


var.indx = function(data, var.name){
  if(length(var.name) != 1){
    stop('var.indx only operates on one variable at a time')
  }
  
  header = names(data)  
  
  indx = grep(paste(var.name, '[\\w_]?', sep=''), header, ignore.case=TRUE)
  
  return(indx)
}




# =============================================
# = Function to predict dimensions of merge() =
# =============================================
# RDB 07-May-2014
# Predict merge() dimensions
# Both desired (no duplicates) and expected (how merge will behave) dimensions
# 'all' argument corresponds to the argument of the same name in merge()
pred.merge <- function(x1, x2, all=FALSE){
	common.names <- intersect(names(x1), names(x2))
	fact1 <- do.call(paste, as.list(x1[,common.names])) #'factors' from x1
	fact2 <- do.call(paste, as.list(x2[,common.names])) # factors from x2

	fInt <- intersect(fact1, fact2) # common elements of fact1 and fact2, same as desired.aF (see below)

	o1 <- table(fact1[fact1%in%fInt])
	o2 <- table(fact2[fact2%in%fInt])
	out.aF <- sum(o1*o2)

	if(all){ # if you used all=TRUE in merge()		
		just1 <- sum(fact1%in%setdiff(fact1, fact2))
		just2 <- sum(fact2%in%setdiff(fact2, fact1))
	
		out.aT <- just1 + just2 + out.aF # This is what merge will give when all=TRUE
		desired.aT <- sum(!duplicated(c(fact1, fact2))) # non-duplicated output if all=TRUE
	
		list(desired=desired.aT, merge=out.aT)
	
	}else{
		out.aF # This is what merge *will* output
		desired.oF <- length(fInt) # Length of desired output, assuming you don't want duplicated join rows
	
		list(desired=desired.oF, merge=out.aF)
	}
}
