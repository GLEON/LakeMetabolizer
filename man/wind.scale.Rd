\name{wind.scale.base}
\alias{wind.scale.base}
\alias{wind.scale}

\title{
	Wind Scaling U10 - exponential conversion to 10m wind speed
}
\description{
	Scale wind speed to standard U10 (10 meters) based on height of observations
}
\usage{
	# Used for timeseries data in a data.frame 
	wind.scale(ts.data, wndZ)
	
	# Used for raw numeric data
	wind.scale.base(wnd, wndZ)
}

\arguments{
	\item{ts.data}{Object of class \code{data.frame} containing a \code{wnd} column.}
  \item{wnd}{measured wind speed (Units: typically m s-1, but it is unit agnostic)}
  \item{wndZ}{height of anemometer (Units: meters)}
}
\details{
	This transforms wind speed to the standard U10, speed at 10 meters, based on the common exponential wind profile assumption. \code{wind.scale} defaults to using the supplied \code{wndZ} value. If \code{wndZ} is not supplied, it attempts to determine the anemometer height from the suffix of the header (e.g., a header of \code{wnd_3} would mean an anemometer height of 3 meters). 
}
\value{
	## \code{wind.scale} 
	
	Returns a data frame with columns \code{datetime} and \code{wnd_10} and the same number of rows as \code{ts.data}
	
	## \code{wind.scale.base} 
	
	Returns a vector with the same length as \code{wnd}
}

\references{
	Saucier, W. 2003. \emph{Principles of Meteorological Analysis}. Dover Publications. New York. p433
}
\author{
	Aline Jaimes, Luke Winslow
}

\seealso{
	Models of gas flux \code{\link{k.cole}}, \code{\link{k.crusius}}, \code{\link{k.macIntyre}}, \code{\link{k.read}}.
}
\examples{
wndSpeed <- c(5.1,6.3,6.3,5.2,7,7.2)
wndHeight <- 2

wind.scale.base(wndSpeed, wndHeight)

}

\keyword{arith}
