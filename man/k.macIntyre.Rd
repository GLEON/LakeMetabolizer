\name{k.macIntyre}
\alias{k.macIntyre}
\title{k600 from MacIntyre et al. 2010}
\description{k600 from MacIntyre et al. 2010, returns in m/day}
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

\value{
  A number for k600 in m/day
}
\author{
  Hilary Dugan
}

\seealso{
\code{\link{k.read}}, \code{\link{k.cole}}, \code{\link{k.macIntyre}}
}

\examples{
	# this is incomplete...
  k.macintyre(wndZ=2, Kd=0.5, atm.press=1013, dateTime, wtr, depth, airT, Uz, RH, sw, lwnet, par)
}

\references{
  MacIntyre, Sally, Werner Eugster, and George W. Kling. The critical importance of buoyancy flux for gas flux across the air-water interface. Geophysical monograph series 127 (2002): 135-139.
}

\keyword{GLEON} 
 
