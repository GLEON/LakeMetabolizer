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
Benson, B. B., & Krause, D. (1984). \emph{The concentration and isotopic fractionation of oxygen dissolved in freshwater and seawater in equilibrium with the atmosphere.} Limnology and Oceanography, 29(3), 620–632. doi:10.4319/lo.1984.29.3.0620

Weiss, R. (1970). \emph{The solubility of nitrogen, oxygen and argon in water and seawater.} Deep Sea Research and Oceanographic Abstracts, 17(4), 721–735. doi:10.1016/0011-7471(70)90037-9

}
\author{
	Luke Winslow
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{water.density}}
}
\examples{
	temp.range = 1:25
	sal.range = 1:25
	
	par(mfrow=c(1,2))
	plot(temp.range, o2.at.sat(temps), xlab='Temperature (C)', ylab='Oxygen Saturation (mg/L)')
	plot(o2.at.sat(rep(20,25), salinity=sal.range), xlab='Salinity (PSU)', ylab='')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{arith}
