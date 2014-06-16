


sw.to.par.base = function(sw){
  
  par = sw * 0.473
}

#'@name sw.to.par
#'@aliases sw.to.par.base
#'@title Convert shortwave radiation to PAR
#'@description 
#'Returns PAR by converting incoming shortwave radiation measuremt.
#'
#'@usage
#'sw.to.par(data, sw.col='sw')
#'
#'sw.to.par.base(sw)
#'
#'@param data Object of class data.frame with column name \code{sw} (or specified alternate)
#'@param sw.col Name of column containing shortwave data (units must be W/m^2)
#'@param sw Numeric shortwave value in W/m^2
#'
#'@return 
#'#For sw.to.par
#'
#'Object of class data.frame with column name 'par' and other values from \code{ts.data}
#'
#'#for sw.to.par.base
#'
#'Numeric vector of PAR values in units umol/m^2/sec
#'
#'@keywords methods math
#'@references
#'Britton, C. M., and J. D. Dodd. \emph{Relationships of photosynthetically active radiation and shortwave irradiance.} 
#'Agricultural Meteorology 17, no. 1 (1976): 1-7.
#'@author
#'Luke Winslow and others
#'@seealso \link{par.to.sw}
#'@examples 
#'#For base function
#'sw <- 800
#'sw.to.par.base(sw)
#'
#'@export
sw.to.par = function(data, sw.col='sw'){
  
  output = data
  
  indx = var.indx(data, sw.col)
  
  data[,indx] = sw.to.par.base(data[,indx])
  names(data)[indx] = 'par' #rename to par
  
  return(data)
}

