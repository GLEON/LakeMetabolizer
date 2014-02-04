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
  \item{wnd}{ Numeric value of wind speed in m/s}
}

\value{
Numeric value of gas exchange velocity for O2 in units of m/(timeStep*min)
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

\keyword{GLEON}
