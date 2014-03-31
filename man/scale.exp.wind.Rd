\name{scale.exp.wind}
\alias{scale.exp.wind}
\title{
	Wind Scaling U10 - exponential conversion to 10m wind speed
}
\description{
	Scale wind speed to standard U10 (10 meters) based on height of observations
}
\usage{
	scale.exp.wind(wnd, wndHeight)
}

\arguments{
  \item{wnd}{measured wind speed (typically m s-1, but it is unit agnostic)}
  \item{wndHeight}{anemometer height in meters}
}
\details{
	This transforms wind speed to the standard U10, speed at 10 meters, based on the common exponential wind profile assumption.
}
\value{
	U10 wind speed at 10m in ms-1
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

scale.exp.wind(wndSpeed, wndHeight)

}

\keyword{arith}
