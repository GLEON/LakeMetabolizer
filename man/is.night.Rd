\name{is.night}
\alias{is.night}
\title{determines if measurement was taken during the night}
\usage{
  is.night(datetimes, lat)
}
\arguments{
  \item{datetimes}{Vector of dates as \code{POSIXct} or
  \code{POSIXlt} (see \code{\link{DateTimeClasses}})
  format}

  \item{lat}{Single latitude value of site. South should be
  negative, north positive}
}
\value{
  a boolean vector of same length as \code{datetimes}
}
\description{
  determines if measurement was taken during the nighttime
}
\author{
  Luke A. Winslow
}
\seealso{
  \link{is.day} \link{sun.rise.set}
}
\keyword{methods}

