\name{temp.kalman}
\alias{temp.kalman}
\title{Smooth temperature time series using a Kalman filter/ smoother}
\usage{
  temp.kalman(wtr, watts, ampH=1, ...)
}
\arguments{
  \item{wtr}{Vector (regular time series) of water
  temperature in degrees C}

  \item{watts}{estimate of watts entering the layer at each
  time step, from \link{watts.in}}

  \item{ampH}{some parameter that}

  \item{...}{parameters to be passed to \link{optim}}
}
\value{
  a smoothed temperature time series
}
\description{
  Smoothes a temperature time series uses a Kalman filter/
  smoother.
}
\details{
  basic model process is \code{x[t] = beta*x[t-1] +
  c1*watts[t-1]}
}
\author{
  Ryan Batt
}
\references{
  Batt, Ryan D. and Stephen R. Carpenter. 2012.
  \emph{Free-water lake metabolism: addressing noisy time
  series with a Kalman filter}. Limnology and Oceanography:
  Methods 10: 20-30. doi: 10.4319/lom.2012.10.20
}
\seealso{
  \link{watts.in} \link{metab.kalman}
}

