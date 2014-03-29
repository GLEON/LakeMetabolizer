\name{scale.exp.wind}
\alias{scale.exp.wind}
\title{Wind Scaling U10 - true 10m wind speed  ms-1
}
\description{Limits of wind speed at 10m in ms-1
}
\usage{
scale.exp.wind(wnd, wndHeight)
}

\arguments{
  
  \item{wnd= measured wind speed ms-1}
  \item{wndHeight=anemometer height in m
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
wndHeight <-2

    U10 <- wnd * (10/wndHeight)^1/7
    return(U10)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{WindScale}
\keyword{GLEON}% __ONLY ONE__ keyword per line
