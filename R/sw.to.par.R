


sw.to.par.base = function(sw){
  
  par = sw * 0.473
}

#'@title Convert shortwave radiation to PAR
#'@description 
#'Returns PAR by converting incoming shortwave radiation measuremt.
#'
#'@usage
#'calc.lw.net(ts.data, lat, atm.press)
#'
#'calc.lw.net.base(dateTime, sw, Ts, lat, atm.press, airT, RH)
#'
#'@param \code{ts.data} Object of class data.frame with column name 'sw'
#'
#'@return Object of class data.frame with column name 'par' and other values from \code{ts.data}
#'
#'@keywords methods math
#'@references
#'Britton, C. M., and J. D. Dodd. \emph{Relationships of photosynthetically active radiation and shortwave irradiance.} 
#'Agricultural Meteorology 17, no. 1 (1976): 1-7.
#'@author
#'LakeMetabolizer
#'@seealso \link{par.to.sw}
#'@examples 
#'sw <- 800
#'sw.to.par.base(sw)
#'@export
sw.to.par = function(data, sw.col='sw'){
  
  output = data
  
  indx = var.indx(data, 'sw')
  
  data[,indx] = sw.to.par.base(data[,indx])
  names(data)[indx] = 'par' #rename to par
  
  return(data)
}

