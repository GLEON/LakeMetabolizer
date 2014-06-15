
#default functions
par.to.sw.base = function(par){
  
  sw = par* 2.114
}


#'@name par.to.sw
#'@aliases par.to.sw.base
#'
#'@title Convert PAR to shortwave
#'@description 
#'Returns incoming shortwave radiation by converting PAR measuremt.
#'
#'@usage
#'par.to.sw.base(par)
#'
#'par.to.sw(data, par.col='par')
#'
#'@param data Object of class data.frame with column name 'par' (units umol/m^2/sec)
#'@param par.col String of alternative name for PAR column 
#'@param par Numeric vector of PAR values
#'
#'
#'@return 
#'#For par.to.sw
#'
#'Object of class data.frame with column name 'sw' and other values from \code{ts.data}
#'
#'#For par.to.sw.base
#'
#'Numeric vector of shortwave values with units W/m^2
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
  
  indx = var.indx(data, par.col)
  
  data[,indx] = par.to.sw.base(data[,indx])
  names(data)[indx] = 'sw' #rename to par
  
  return(data)
}


