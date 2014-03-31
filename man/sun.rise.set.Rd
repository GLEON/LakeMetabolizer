\name{sun.rise.set}
\alias{sun.rise.set}
\title{sun.rise.set:Calculates the time of sunrise and sunset}
\description{Takes lattitude (deg south should be negative, north positive) and a
matlab datenum (days since year 0) and returns the time of sunrise as a
day fraction (e.g., 0.5 would be noon, 0 would be midnight)
\usage{
sun.rise.set(lat, datetimes)
}

\arguments{
  \item{lat}{The latitude for the desired location on earth}
}
  \item{datetimes}{Either Matlab datenum or day-of-year (Jan 1st is 1)}

}
}

\value{
\item{is.daytime}{time of sunrise and sunset}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Aline Jaimes
}
\note{

}
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
\keyword{arith}
