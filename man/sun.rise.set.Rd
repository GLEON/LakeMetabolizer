\name{sun.rise.set}
\alias{sun.rise.set}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sun.rise.set(lat, datetimes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lat}{
%%     ~~Describe \code{lat} here~~
}
  \item{datetimes}{
%%     ~~Describe \code{datetimes} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (lat, datetimes) 
{
    doy = as.POSIXlt(datetimes)$yday
    dayAngle = 2 * pi * (doy - 1)/365
    degToRad = 2 * pi/360
    radToDeg = 180/pi
    dec = 0.006918 - 0.399912 * cos(dayAngle) + 0.070257 * sin(dayAngle) - 
        0.006758 * cos(2 * dayAngle) + 0.000907 * sin(2 * dayAngle) - 
        0.002697 * cos(3 * dayAngle) + 0.00148 * sin(3 * dayAngle)
    latRad = lat * degToRad
    sunriseHourAngle = acos(-tan(latRad) * tan(dec)) * radToDeg
    sunriseHourAngle[is.complex(sunriseHourAngle)] = NA
    sr = 12 - sunriseHourAngle/15
    ss = 12 + sunriseHourAngle/15
    rise = trunc(datetimes, "day") + sr * 60 * 60
    set = trunc(datetimes, "day") + ss * 60 * 60
    return(as.POSIXct(matrix(c(rise, set), ncol = 2), origin = "1970-01-01"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
