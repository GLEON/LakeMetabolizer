\name{k.crusius}
\alias{k.crusius}
\title{k600 from Crusius and Wanninkhof 2003}
\description{k600 from Crusius and Wanninkhof 2003, returns in m/day}
\usage{
  k.crusius(wnd,method)
}
\arguments{
  \item{wnd}{numeric value of wind speed, m/s}
  \item{method}{string value of valid method name (default is "power")}}

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

  k.crusius(wnd=2.3,method="bilinear")
}

\references{
  Crusius, John, and Rik Wanninkhof. "Gas transfer velocities measured at low wind speed over a lake." Limnology and Oceanography 48, no. 3 (2003): 1010-1017.
}

\keyword{GLEON} 
