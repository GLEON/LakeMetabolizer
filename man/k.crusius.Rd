\name{k.crusius}
\alias{k.crusius}
\title{Returns k600 from Crusius and Wanninkhof 2003}
\description{k600 from Crusius and Wanninkhof 2003, returns in m/day}
\usage{
  k.crusius(wnd,method)
}
\arguments{
  \item{wnd}{numeric vector of wind speeds, m/s}
  \item{method}{string of valid method name. Either "linear", 
  "bilinear", or "power" (default is "power")}}

\value{
  k600 vector in m/day
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
  Crusius, John, and Rik Wanninkhof. Gas transfer velocities measured at low wind speed over a lake. 
  Limnology and Oceanography 48, no.3 (2003): 1010-1017.
}

\keyword{methods}
\keyword{math}
