\name{scale.exp.wind}
\alias{scale.exp.wind}
\title{Wind Scaling U10 - true 10m wind speed  m s-1
}
\description{Limits of wind speed at 10m in m s-1
}
\usage{
scale.exp.wind(wnd, wndHeight)
}

\arguments{
  
  \item{wnd}{measured wind speed m s-1}
  \item{wndHeight}{anemometer height in m}
}
\details{
}
\value
\item{U10 wind speed at 10m in ms-1}
\references{
Saucier, W. 2003.Principles of Meteorological Analysis. Dover Publications. New York. p433
}
\author{Aline Jaimes
}
\note{
}

\seealso{

}
\examples{
wndSpeed <- c(5.1,6.3,6.3,5.2,7,7.2)
wndHeight <- 2

scale.exp.wind(wndSpeed, wndHeight)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{arith}
