\name{k.cole}
\alias{k.cole}
\title{#short description of function#
  Returns the gas exchange velocity for O2 based on Cole and Caraco, 1998 
}
\description{#more detailed description of function# 
  Returns the gas exchange velocity for O2 based on Cole and Caraco, 1998 
  in units of m/(timeStep*min) (i.e. 30 minute sampling interval will return 
  kO2 in units of m/(1/48) - converts to fraction of day)
}
\usage{
k.cole(wnd,wtr,wndHeight)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wnd}{ Numeric value of wind speed in m/s}
  \item{wtr}{ Numeric value of water temperature at dissolved oxygen sensor depth in degrees celsius}
  \item{wndHeight}{ Numeric value of wind height measurement in meters}
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
wnd <- 4.3
wtr <- 23.4
wndHeight <- 2
timeStep <- 30 #sampling interval in minutes 
k.cole(wnd,wtr,wndHeight)
}

\keyword{GLEON}
