\name{k.cole}
\alias{k.cole}
\title{Returns the gas exchange velocity for k600 based on Cole and Caraco, 1998 
}
\description{
  Returns the gas exchange velocity for k600 based on Cole and Caraco, 1998 
  in units of m/day)
}
\usage{
	k.cole(wnd)
}
\arguments{
  \item{wnd}{ Numeric vector of wind speed in m/s}
}

\value{
Numeric vector of gas exchange velocities for O2 in units of m/(timeStep*min)
}

\author{
Jake Zwart
}

\seealso{
\code{\link{k.read}}, \code{\link{k.crusius}}, \code{\link{k.macIntyre}}
}

\examples{
k.cole(wnd=4.3)
}

\references{
  Cole, J.J. and Caraco, N.F. 1998. Atmospheric exchange of carbon dioxide in a low-wind oligotrophic 
  lake measured by the addition of SF6. Limnology and Oceanography. 43(4): 647-656. 
}

\keyword{methods}
\keyword{math}
