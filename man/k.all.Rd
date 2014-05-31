\name{k.cole}
\alias{k.cole}
\alias{k.read}
\alias{k.crusius}
\alias{k.macIntyre}

\title{Returns a timeseries of gas exchange velocity for k600 based on a specific model 
}
\description{
  Returns the gas exchange velocity for k600 based on a specific model
  in units of m/day)
}
\usage{
  ## Method for Cole and Caraco, 1998
  k.cole(ts.data)
  
  ## Method for Crusius and Wanninkhof 2003
  k.crusius(ts.data, method='power')
  
  ## Method for Read et al 2012
  k.read(ts.data, wndZ, Kd, atm.press, lat, lake.area)
  
  ## Method for from MacIntyre et al. 2010
  k.macIntyre(ts.data, wndZ, Kd, atm.press)
  
}
\arguments{
  \item{ts.data}{Object of class \code{data.frame} }
  \item{method}{Only for \code{k.crusius.base}. String of valid method . Either "linear", "bilinear", or "power"}
  \item{wndZ}{Height of wind measurement, (Units: m)}
  \item{Kd}{Light attenuation coefficient (Units: m^-1)}
  \item{atm.press}{Atmospheric pressure, (Units: millibar)}  
  \item{lat}{Latitude, degrees north}
  \item{lake.area}{Lake area, m^2}
}

\value{
## When called on a numeric vector
Numeric vector of gas exchange velocities for O2 in units of m/(timeStep*min)

## When called on a \code{data.frame}
data.frame with 
}

\author{
Hilary Dugan, Jake Zwart, Luke Winslow
}

\seealso{
\code{\link{k.read}}, \code{\link{k.crusius}}, \code{\link{k.macIntyre}}
}

\examples{

}

\references{
  Cole, J.J. and Caraco, N.F. 1998. Atmospheric exchange of carbon dioxide in a low-wind oligotrophic 
  lake measured by the addition of SF6. Limnology and Oceanography. 43(4): 647-656. 
}

\keyword{methods}
\keyword{math}
