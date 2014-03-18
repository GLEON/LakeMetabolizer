\name{k.macIntyre}
\alias{k.macIntyre}
\title{k600 from MacIntyre et al. 2010}
\description{k600 from surface renewal model described in MacIntyre et al. 2010, returned in m/day}
\usage{
  k.crusius(wndZ, Kd, atm.press, dateTime, wtr, depth, airT, Uz, RH, sw, lwnet, par)
}
\arguments{
  \item{wndZ}{numeric value wind height in m}
  \item{Kd}{numeric value light attenuation coefficient m-1}
  \item{atm.press}{Atmospheric pressure in hPa or Mb. If unknown use 1013 for sea level}
  \item{dateTime}{date and time vector}
  \item{wtr}{dataframe of water temperatures}
  \item{depth}{vector of water temperature depths}
  \item{airT}{vector of air temperatures in °C}
  \item{Uz}{vector of wind speed in m/s}
  \item{RH}{vector of relative humidity in percent (0-100)}
  \item{sw}{vector of shortwave radiation in W/m2}
  \item{lwnet}{vector of longwave radiation. If missing, run calc.lw.net function first}
  \item{par}{vector of par data in umol m-2 s-1 }
}
\value{
  A vector for k600 in m/day
}
\author{
  Hilary Dugan
}

\seealso{
\code{\link{k.read}}, \code{\link{k.cole}}, \code{\link{k.crusius}}
}

\examples{
  wndZ <- 2
  Kd <- 2
  atm.press <- 1013
  dateTime <- as.POSIXct("2013-12-30 14:00")
  wtr <- c(22.51, 22.42, 22.4, 22.4, 22.4, 22.36, 22.3, 22.21, 22.11, 21.23, 16.42)
  depth <- c(0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8)
  airT <- 20
  Uz <- 6
  RH <- 90
  sw <- 800
  lwnet <- -55
  k500_macInytre <- k.macIntyre(wndZ, Kd, atm.press, dateTime, wtr, depth, airT, Uz, RH, sw, lwnet)
}

\references{
 MACINTRYE, sALLY, ANDERS JONSSON, MATS JANSSON, JAN ABERG, DAMON E. TURNEY AND SCOTT D. MILLER.
 2010. Buoyancy flux, turbulence, and the gas transfer coefficient in a stratified lake.
 Geophysical Research Letters. 37: L24604. doi:10.1029/2010GL044164 
}

\keyword{methods}
\keyword{math}
 
