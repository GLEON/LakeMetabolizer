\name{sun.rise.set}
\alias{sun.rise.set}
\title{
Calculates the time of sunrise and sunset
}

\description{
Calculates the time of sunrise and sunset based on latitude and date.
}
\usage{
sun.rise.set(lat, datetimes)
}

\arguments{
	\item{lat}{Single latitude value of site. South should be negative, north positive}
	\item{datetimes}{Vector of dates as POSIXct or POSIXlt (see \code{\link{DateTimeClasses}}) format}
}

\value{
	A 2-column matrix, first column sunrise, second column sunset, as \code{\link{POSIXct}} format.
}
\references{
	Iqbal, Muhammad. 1983. \emph{An Introduction to Solar Radiation}. Elsevier.
}
\author{
	Aline Jaimes, Luke Winslow
}

\examples{

sun.rise.set(lat=40.75,datetimes=as.POSIXlt('2013-03-31'))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
\keyword{arith}
