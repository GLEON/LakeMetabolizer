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
Returns a data.frame with a \code{datetime} column and a \code{k600} column. k600 is in units of meters per day (m/d).
}

\author{
Hilary Dugan, Jake Zwart, Luke Winslow, R. Iestyn. Woolway
}

\seealso{
\code{\link{k.read}}, \code{\link{k.crusius}}, \code{\link{k.macIntyre}}
}

\examples{
  data.path = system.file('extdata', '', package="LakeMetabolizer")
  
  tb.data = load.all.data('sparkling', data.path)
  
  ts.data = tb.data$data #pull out just the timeseries data
  
  #calculate U10 and add it back onto the original 
  
  u10 = scale.exp.wind(ts.data)
  ts.data = rmv.var(ts.data, 'wnd', ignore.offset=TRUE) #drop old wind speed column
  ts.data = merge(ts.data, u10)                          #merge new u10 into big dataset  
  
  
  k600_cole = k.cole(ts.data)
  
  k600_crusius = k.crusius(ts.data)
  
  kd        = tb.data$metadata$averagekd
  wndZ      = 10   #because we converted to u10
  atm.press	= 1018
  lat       = tb.data$metadata$latitude
  lake.area = tb.data$metadata$lakearea
  
  #for k.read and k.macIntyre, we need LW_net. Calculate because we don't have obs
  lwnet = calc.lw.net(ts.data, lat, atm.press)
  
  ts.data = merge(ts.data, lwnet)
  
  k600_read = k.read(ts.data, wndZ=wndZ, Kd=kd, atm.press=atm.press, lat=lat, lake.area=lake.area)
  
  k600_macIntyre = k.macIntyre(ts.data, wndZ=wndZ, Kd=kd, atm.press=atm.press)
  
}

\references{
	See \code{\link{k.cole}} for references.
}

\keyword{methods}
\keyword{math}
