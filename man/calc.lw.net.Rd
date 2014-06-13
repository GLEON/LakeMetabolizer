\name{calc.lw.net}
\alias{calc.lw.net}
\title{Estimate net long wave heat radiation}
\usage{
  calc.lw.net(ts.data, lat, atm.press)

  calc.lw.net.base(dateTime, sw, Ts, lat, atm.press, airT,
    RH)
}
\arguments{
  \item{\code{ts.data}}{Object of class data.frame}

  \item{\code{dateTime}}{vector of datetime in POSIXct
  format}

  \item{\code{sw}}{numeric value of short wave radiation,
  W/m2}

  \item{\code{Ts}}{numeric value of surface water
  temperature, degC}

  \item{\code{lat}}{latitude in degrees north}

  \item{\code{atm.press}}{atmospheric pressure in mb}

  \item{\code{airT}}{numeric value of air temperature,
  degC}

  \item{\code{RH}}{numeric value of relative humidity, \%}
}
\value{
  A numeric value of net long wave heat flux in W/m2
}
\description{
  Returns the net long wave radiation based on Crawford and
  Duchon, 1999.
}
\examples{
dateTime <- as.POSIXct("2013-12-30 23:00")
Uz <- 3
airT <- 20
RH <- 90
sw <- 800
wndZ <- 2
Kd <- 2
lat <- 54
lake.area <- 5000
atm.press <- 1013
Ts <- 22
calc.lw.net.base(dateTime,sw,Ts,lat,atm.press,airT,RH)
}
\author{
  R Iestyn Woolway Jordan S. Read Hilary Dugan
}
\references{
  Crawford, T.M., and Duchon, C.E. 1999. An improved
  parameterization for estimating effective atmospheric
  emissivity for use in calculating daytime downwelling
  longwave radiation. Journal of Applied Meteorology 38:
  474-480.
}
\seealso{
  k.read
}
\keyword{math}
\keyword{methods}

