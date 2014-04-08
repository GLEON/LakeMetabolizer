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