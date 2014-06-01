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

rmv.var = function(data, var.name, ignore.missing=TRUE, ignore.offset=FALSE){
	if(ignore.offset){
		varI = var.indx(data, var.name)
	}else{
		varI = grep(var.name, names(data), ignore.case=TRUE)
	}
	
	if(length(varI) > 0){
		varI = varI * -1
		return(data[, varI])
	}else{
		if(!ignore.missing){
			stop('No variable by that name found')
		}
	}

}


var.indx = function(data, var.name){
  if(length(var.name) != 1){
    stop('var.indx only operates on one variable at a time')
  }
  
  header = names(data)  
  
  indx = grep(paste(var.name, '[\\w_]?', sep=''), header, ignore.case=TRUE)
  
  return(indx)
}



# test edit
# =============================================
# = Function to predict dimensions of merge() =
# =============================================
# RDB 07-May-2014
# Predict merge() dimensions
# Both desired (no duplicates) and expected (how merge will behave) dimensions
# 'all' argument corresponds to the argument of the same name in merge()
pred.merge <- function(x1, x2, all=FALSE){
	common.names <- intersect(names(x1), names(x2))
	
	if(length(common.names)>1){
		fact1 <- do.call(paste, as.list(x1[,common.names])) #'factors' from x1
		fact2 <- do.call(paste, as.list(x2[,common.names])) # factors from x2	
	}else{
		fact1 <- x1[,common.names]
		fact2 <- x2[,common.names]
	}


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



# =======================
# = round.time Function =
# =======================
# RDB 16May2014
# example: round.time(x, "5 minutes")
# example: round.time(x, "90 min")
# example: round.time(x, "0.5 hours")
# x is a time format – preferably POSIXct, or can be coerced w/ as.POSIXct
# if x needs to be converted to POSIX, define input.format if x currently isn't in a 'standard unambiguous format'
round.time <- function(x, units, input.format=NULL, output.format="%Y-%m-%d %H:%M:%S"){
	# x = head(t.sonde0.na2[,"date"], 20) + 120
	# units = "df.345 min"
	# Check for invalid input classes
	stopifnot(is.character(units)&is.character(output.format)&(is.null(input.format)|is.character(input.format)))
	
	# Determine time unit
	unit.choices <- c("sec", "min", "hour", "day")
	choices.or <- paste(unit.choices, collapse="|")
	unit.pattern <- paste(".*(", choices.or, ").*", sep="")
	unit <- gsub(unit.pattern, "\\1", units)
	if(is.na(unit)){stop("not a valid unit, use sec, min, hour, or day")}
	which.choice <- which(unit==unit.choices)
	
	# Determine time interval
	u.time.pattern <- "(?:[0-9]+\\.[0-9]+)|(?:[0-9]+\\.)|(?:\\.[0-9]+)|(?:[0-9]+)"
	u.time.char <- regmatches(units, regexpr(u.time.pattern, units, perl=TRUE))
	u.time <- as.numeric(u.time.char)
	u.time <- ifelse(is.na(u.time), 1, u.time)
	
	unit.cutoff <- switch(unit, sec=60, min=60, hour=24, day=1)
	
	# =========================================================================
	# = Check for invalid input (before slow [attempted] conversion to POSIX) =
	# =========================================================================
	if(sign(u.time)==-1L){
		stop("time interval must be positive")
	}
	# Deal with case where units are 1 second (or less)
	if(unit=="sec" & u.time<=1L){
		return(format.Date(x, format=output.format))
	} else
	
	# Fractional time intervals – convert to smaller unit
	if((trunc(u.time)-u.time)!=0){
		if(sign(u.time)==1L){
			while((trunc(u.time)-u.time)!=0){
				if(unit=="sec"){stop("time interval must be an integer when converted to units of seconds")}
				unit <- unit.choices[which.choice-1]
				which.choice <- which(unit==unit.choices)
				unit.cutoff <- switch(unit, sec=60, min=60, hour=24)
				u.time <- unit.cutoff*u.time
			}
		}else{
			stop("time interval must be positive")
		}
	} else 
	
	# Deal with case where units are days
	if(unit=="day"){
		if(u.time==1){
			return(format.Date(trunc.POSIXt(x + 43200, units = units), format=output.format))
		}else{
			stop("units must be <= 1 day")
		}
	} else 
	
	# Deal w/ cases where time interval is 1 unit
	if(u.time==1){
			unit <- unit.choices[which.choice-1]
			which.choice <- which(unit==unit.choices)
			unit.cutoff <- switch(unit, sec=60, min=60, hour=24)
			u.time <- unit.cutoff
	} 
	
	# Deal with cases where time interval is > 1 of a larger unit
	# Note that this follows up on case where u.time is > 1 and is not an integer
	if(u.time>unit.cutoff){
		u.time <- u.time%%unit.cutoff
		mod.mess <- paste("Rounding to units =", u.time, unit) # may or may not want to make this a warning ...
		warning(mod.mess)
	}
	
	# =============================================
	# = Convert to POSIX, or if can't, give error =
	# =============================================
	if(all(class(x)!="POSIX.ct")){
		if(is.null(input.format)){
			x <- as.POSIXct(x)
		}else{
			x <- as.POSIXct(x, format=input.format)
		}
	}
	
	# ===========================================================
	# = Matching units (e.g., min) and unit multiples (e.g., 5) =
	# ===========================================================
	
	which.choice <- which(unit==unit.choices)
	form.unit <- c("%S", "%M", "%H", "%d")[which.choice]
	mult <- as.integer(format.Date(x, format=form.unit))/u.time
	after <- round(mult, 0)*u.time
	# direction <- sign(after-before)
	
	# trunc.unit <- unit.choices[min(which.choice+1, length(unit.choices))]
	trunc.unit <- unit.choices[min(which.choice+1, length(unit.choices))]
	rounded <- trunc.POSIXt(x, trunc.unit) + switch(unit, sec = 1, min = 60, hour = 3600, day = 86400)*after
	format.Date(rounded, format=output.format)
}