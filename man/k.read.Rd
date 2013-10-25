\name{k.read}
\alias{k.read}
\title{Returns the gas exchange velocity for O2 based on Read et al, 2012}
\description{Returns the gas exchange velocity for O2 based on Read et al, 2012 
  in units of m/(timeStep*min) (i.e. 30 minute sampling interval will return 
  kO2 in units of m/(1/48) - converts to fraction of day}
\usage{
  k.read <- function(wndZ,Kd,lat,lake.area,atm.press,dateTime,wtr,airT,Uz,RH,sw,lwnet,par,lw)
}
\arguments{
  \item{wndZ}{Height of wind measurement, m}
  \item{Kd}{Light attenuation coefficient}
  \item{lat}{Latitude, degrees north}
  \item{lake.area}{Lake area, m2}
  \item{atm.press}{Atmospheric pressure, mb}  
  \item{dateTime}{datetime (Y-\%m-\%d \%H:\%M), POSIXct format}
  \item{wtr}{Numeric vector of water temperature, degC}
  \item{depths}{Numeric vector of depths. Must be the same length as the wtr parameter}  
  \item{airT}{Numeric value of air temperature, degC}
  \item{Uz}{Numeric value of wind speed, m/s}
  \item{RH}{Numeric value of relative humidity, \%}
  \item{sw}{Numeric value of short wave radiation, W/m2} 
  \item{lwnet}{Numeric value net long wave radiation, W/m2}  
  \item{par}{Numeric value of photosynthetically active radiation, umol m-2 s-1 (if sw not available)}  
  \item{lw}{Numeric value of incoming long wave radiation, W/m2}  
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
  wndZ <- 2
  Kd <- 2
  lat <- 54
  lake.area <- 5000 
  atm.press <- 1013
  dateTime <- as.POSIXct("2013-12-30 14:00")
  wtr <- c(22.51, 22.42, 22.4, 22.4, 22.4, 22.36, 22.3, 22.21, 22.11, 21.23, 16.42)
  depth <- c(0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8)
  airT <- 20
  Uz <- 3
  RH <- 90
  sw <- 800
  lwnet <- -55
  timeStep <- 30
  ko2_read <- k.read(wndZ,Kd,lat,lake.area,atm.press,dateTime,wtr,depth,airT,Uz,RH,sw,lwnet)
}

\references{
  Read, J.S., Hamilton, D.P., Desai, A.R., Rose, K.C., MacIntyre, S.,
  Lenters, J.D., Smyth, R.L., Hanson, P. .C., Cole, J.J., Staehr, P.A.,
  Rusak, J.A., Pierson, D.C., Brookes, J.D., Laas, A., Wu, C.H. 2012.
  \emph{Lake-size dependency of wind shear and convection as controls on gas
  exchange}. Geophysical Research Letters 39, DOI: 10.1029/2012GL051886.
}

\keyword{GLEON}
