\name{o2.at.sat}
\alias{o2.at.sat}

\title{
Calculates the equilibrium saturation concentration of oxygen in water at the supplied conditions
}
\description{
Used to calculate the equilibrium concentration of oxygen in water. The equilibration concentration of oxygen 
in water varies with both temperature, salinity, and the partial pressure of oxygen in contact with the water 
(calculated from supplied elevation or barometric pressure). 
}
\usage{
o2.at.sat(temp, baro, altitude = 0, salinity = rep(0, length(temp)), model = "weiss")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{temp}{
	a numeric vector of water temperature in degrees Celsius.
}
  \item{baro}{
	barometric pressure in millibar.
}
  \item{altitude}{
	a numeric value indicating the elevation above mean sea level in meters. Defaults to mean sea level. An alternative to supplying barometric pressure.
}
  \item{salinity}{
	a numeric vector of salinity in PSU. Defaults to zero. Length must be one or equal to length of temperature.
}
  \item{model}{
	the empirical model to be used, \code{"weiss"} and \code{"benson"} are the two available options. They correspond to the references described below.
}

}%%Close arguments section.

\value{
	The equilibration concentration at the supplied conditions in mg/L of oxygen.
}
\references{
	
}
\author{
	Luke Winslow
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{water.density}}
}
\examples{
	temps = 1:25
	sals = 1:25
	
	par(mfrow=c(1,2))
	plot(temps, o2.at.sat(temps), xlab='Temperature (C)', ylab='Oxygen Saturation (mg/L)')
	plot(o2.at.sat(rep(20,25), salinity=sals), xlab='Salinity (PSU)', ylab='')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{arith}
