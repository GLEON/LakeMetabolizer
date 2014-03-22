\name{metab.bookkeep}
\alias{metab.bookkeep}
\title{Returns metabolism estimates for observed system based on Cole et al, 2010}
\description{Returns metabolism estimates for observed system based on Cole et al, 2010 
  in units of concentration/(timeStep) (i.e. 30 minute sampling interval will return 
  kO2 in units of m/(1/48) - converts to fraction of day}
\usage{
  k.read <- function(wndZ,Kd,lat,lake.area,atm.press,dateTime,wtr,airT,Uz,RH,sw,lwnet,par,lw)
}
\arguments{
  \item{dateTime}{datetime (Y-\%m-\%d \%H:\%M), POSIXct format}
}

\value{
  Numeric value of gas exchange velocity for O2 in units of m/(timeStep*min)
}
\author{
  R. Iestyn. Woolway, Hilary Dugan, Luke Winslow, Jordan S Read, GLEON fellows
}

\seealso{
  \code{\link{k.cole}}, \code{\link{k.crusius}}, \code{\link{k.macIntyre}}
}

\examples{
  library(rLakeAnalyzer)
  doobs = load.ts(system.file('extdata', 
			'Sparkling.doobs', package="LakeMetabolizer"))
  wtr = load.ts(system.file('extdata', 
			'Sparkling.wtr', package="LakeMetabolizer"))
  wnd = load.ts(system.file('extdata', 
			'Sparkling.wnd', package="LakeMetabolizer"))
  
  #Subset a day
  mod.date = as.POSIXct('2009-08-12')
  doobs = doobs[trunc(doobs$datetime, 'day') == mod.date, ]
  wtr = wtr[trunc(wtr$datetime, 'day') == mod.date, ]
  wnd = wnd[trunc(wnd$datetime, 'day') == mod.date, ]
  
  k.gas = k600.2.kGAS(k.cole(wnd[,2]), wtr[,3], 'O2')
  do.sat = o2.at.sat(wtr[,3], altitude=300)
  
  metab.bookkeep(doobs[,2], do.sat, k.gas, z.mix=1, doobs[,1], 45)
  
  
}

\references{
  Read, J.S., Hamilton, D.P., Desai, A.R., Rose, K.C., MacIntyre, S.,
  Lenters, J.D., Smyth, R.L., Hanson, P. .C., Cole, J.J., Staehr, P.A.,
  Rusak, J.A., Pierson, D.C., Brookes, J.D., Laas, A., Wu, C.H. 2012.
  \emph{Lake-size dependency of wind shear and convection as controls on gas
  exchange}. Geophysical Research Letters 39, DOI: 10.1029/2012GL051886.
}

\keyword{GLEON}
