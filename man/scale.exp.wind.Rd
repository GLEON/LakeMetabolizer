\name{scale.exp.wind}
\alias{scale.exp.wind}
\title{Wind Scaling U10 - true 10m wind speed  [m/s]
}
\description{Limits of wind speed at 10m in m/s
}
\usage{
scale.exp.wind(wnd, wndHeight)
}

\arguments{
  
  \item{wnd= measured wind speed [m/s]}
  \item{wndHeight=anemometer height in [m]
}
\details{
}
\value
\item{U10 wind speed in 10 m [m/s]}
\references{
Saucier, W. 2003.Principles of Meteorological Analysis.Dover Publications. New York. p433
}
\author{Aline Jaimes
}
\note{
}

\seealso{
http://www.atmos.washington.edu/academics/classes/2013Q2/547/lect5.pdf
}
\examples{
% Ex. 1
test = read.table("scale.exp.wind_Example.csv")  % read file 
U10 <- wnd=(1:10) * (10/wndHeight=3)^1/7
% %%%%%
% Ex. 2
wnd->1:10
wndHeight->3

% Run Function
{
    U10 <- wnd * (10/wndHeight)^1/7
    return(U10)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{WindScale}
\keyword{GLEON}% __ONLY ONE__ keyword per line
