\name{sw.to.par}
\alias{sw.to.par}
\title{Convert shortwave radiation to PAR}
\usage{
  calc.lw.net(ts.data, lat, atm.press)

  calc.lw.net.base(dateTime, sw, Ts, lat, atm.press, airT,
    RH)
}
\arguments{
  \item{\code{ts.data}}{Object of class data.frame with
  column name 'sw'}
}
\value{
  Object of class data.frame with column name 'par' and
  other values from \code{ts.data}
}
\description{
  Returns PAR by converting incoming shortwave radiation
  measuremt.
}
\examples{
sw <- 800
sw.to.par.base(sw)
}
\author{
  LakeMetabolizer
}
\references{
  Britton, C. M., and J. D. Dodd. \emph{Relationships of
  photosynthetically active radiation and shortwave
  irradiance.} Agricultural Meteorology 17, no. 1 (1976):
  1-7.
}
\seealso{
  \link{par.to.sw}
}
\keyword{math}
\keyword{methods}

