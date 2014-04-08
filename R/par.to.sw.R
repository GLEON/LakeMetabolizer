
# s3 dispatchers
par.to.sw = function(x, ...) UseMethod("par.to.sw")

sw.to.par = function(x, ...) UseMethod("sw.to.par")


#default functions
par.to.sw.default = function(par){
  
  sw = par* 2.114
}

sw.to.par.default = function(sw){
  
  par = sw * 0.473
}

#data.frame functions for timeseries handling 
sw.to.par.data.frame = function(data, sw.col='sw'){
  
  output = data
  
  indx = var.indx(data, 'sw')
  
  data[,indx] = sw.to.par(data[,indx])
  names(data)[indx] = 'par' #rename to par
  
  return(data)
}


par.to.sw.data.frame = function(data, par.col='par'){
  
  output = data
  
  indx = var.indx(data, 'par')
  
  data[,indx] = par.to.sw(data[,indx])
  names(data)[indx] = 'sw' #rename to par
  
  return(data)
}

#Papaioannou, G, N Papanikolaou, and D Retalis. 1993.
#"Relationships of Photosynthetically Active Radiation and Shortwave Irradiance." 
#Theoretical and Applied Climatology 48: 23-27. http://www.sciencedirect.com/science/article/pii/0002157176900807.