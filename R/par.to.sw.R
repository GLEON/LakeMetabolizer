


# s3 dispatchers
#par.to.sw = function(x, ...) UseMethod("par.to.sw")

#sw.to.par = function(x, ...) UseMethod("sw.to.par")


#default functions
par.to.sw.base = function(par){
  
  sw = par* 2.114
}


#'@title Convert PAR to shortwave
#'@description 
#'Returns incoming shortwave radiation by converting PAR measuremt.
#'
#'@usage
#'calc.lw.net(ts.data, lat, atm.press)
#'
#'calc.lw.net.base(dateTime, sw, Ts, lat, atm.press, airT, RH)
#'
#'@param \code{ts.data} Object of class data.frame with column name 'par'
#'
#'@return Object of class data.frame with column name 'sw' and other values from \code{ts.data}
#'
#'@keywords methods math
#'@references
#'Britton, C. M., and J. D. Dodd. \emph{Relationships of photosynthetically active radiation and shortwave irradiance.} 
#'Agricultural Meteorology 17, no. 1 (1976): 1-7.
#'@author
#'LakeMetabolizer
#'@seealso \link{sw.to.par}
#'@examples 
#'par <- 800
#'par.to.sw.base(par)
#'@export
par.to.sw = function(data, par.col='par'){
  
  output = data
  
  indx = var.indx(data, 'par')
  
  data[,indx] = par.to.sw.base(data[,indx])
  names(data)[indx] = 'sw' #rename to par
  
  return(data)
}


