\name{k.read.base}
\alias{k.read.base}
\alias{k.cole.base}
\alias{k.crusius.base}
\alias{k.macIntyre.base}

\title{Returns the gas exchange velocity (k600) for the specified model}

\description{Returns the gas exchange velocity based on the chosen model in units of m/day}
\usage{
  ## Method for Cole and Caraco, 1998
  k.cole.base <- function(wnd)
  
  ## Method for Crusius and Wanninkhof 2003
  k.crusius.base <- function(wnd, method='power')
  
  ## Method for Read et al 2012
  k.read.base <- function(wndZ, Kd, lat, lake.area, atm.press, dateTime, wtr, airT, Uz, RH, sw, lwnet, par, lw)
  
  ## Method for Macintyre et al 2010
  k.macIntyre.base <- function(wndZ, Kd, atm.press, dateTime, surf.temp, z.mix, airT, Uz, RH, sw, lwnet)
}
\arguments{

  \item{wnd}{Numeric value of wind speed, (Units:m/s)}
  \item{method}{Only for \code{k.crusius.base}. String of valid method . Either "linear", "bilinear", or "power"}
  
  \item{wndZ}{Height of wind measurement, (Units: m)}
  \item{Kd}{Light attenuation coefficient (Units: m^-1)}
  \item{lat}{Latitude, degrees north}
  \item{lake.area}{Lake area, m^2}
  \item{atm.press}{Atmospheric pressure, (Units: millibar)}  
  \item{dateTime}{datetime (Y-\%m-\%d \%H:\%M), (Format: \code{\link{POSIXct)}}}
  \item{wtr}{Numeric vector of water temperature, (Units(deg C)}
  \item{depths}{Numeric vector of depths. Must be the same length as the wtr parameter}  
  \item{airT}{Numeric value of air temperature, Units(deg C)}
  
  \item{RH}{Numeric value of relative humidity, \%}
  \item{sw}{Numeric value of short wave radiation, W m^-2} 
  \item{lwnet}{Numeric value net long wave radiation, W m^-2}  
  \item{par}{Numeric value of photosynthetically active radiation, umol m^-2 s^-1 (if sw not available)}  
  \item{lw}{Numeric value of incoming long wave radiation, W m^-2}  
}

\value{
  Numeric value of gas exchange velocity (k600) in units of m/day. Before use, should be converted to appropriate gas using \code{\link{k600.2.kGAS}}. 
}
\author{
  R. Iestyn. Woolway, Hilary Dugan, Luke Winslow, Jordan S Read, GLEON fellows
}

\seealso{
  \code{\link{k.cole}}, \code{\link{k.crusius}}, \code{\link{k.macIntyre}}, \code{\link{k.read}}
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
  Uz <- 6
  RH <- 90
  sw <- 800
  lwnet <- -55
  timeStep <- 30
  
  U10 <- scale.exp.wind(Uz, wndZ)
  
  k600_cole <- k.cole.base(U10)
  
  k600_crusius <- k.crusius.base(U10)
  
  k600_read <- k.read.base(wndZ, Kd, lat, lake.area, atm.press, dateTime, wtr, depth, airT, Uz, RH, sw, lwnet)
  
  k600_macInytre <- k.macIntyre.base(wndZ, Kd, atm.press, dateTime, wtr, depth, airT, Uz, RH, sw, lwnet)
  
}

\references{
  Cole, JJ and Caraco, NF 1998. \emph{Atmospheric exchange of carbon dioxide in a low-wind oligotrophic 
  lake measured by the addition of SF6.} Limnology and Oceanography. 43(4): 647-656. 

  Crusius, J, and R Wanninkhof. \emph{Gas transfer velocities measured at low wind speed over a lake.} 
  Limnology and Oceanography 48, no.3 (2003): 1010-1017.
  
  Read, J.S., Hamilton, D.P., Desai, A.R., Rose, K.C., MacIntyre, S.,
  Lenters, J.D., Smyth, R.L., Hanson, P. .C., Cole, J.J., Staehr, P.A.,
  Rusak, J.A., Pierson, D.C., Brookes, J.D., Laas, A., Wu, C.H. 2012.
  \emph{Lake-size dependency of wind shear and convection as controls on gas
  exchange}. Geophysical Research Letters 39, DOI: 10.1029/2012GL051886.
  
  MacIntyre, S, A. Jonsson, M Jansson, J Aberg, DE Turney AND SD Miller.
  2010. \emph{Buoyancy flux, turbulence, and the gas transfer coefficient in a stratified lake.}
  Geophysical Research Letters. 37: L24604. doi:10.1029/2010GL044164 

}

\keyword{methods}
\keyword{math}
