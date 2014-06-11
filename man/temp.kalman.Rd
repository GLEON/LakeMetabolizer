\name{temp.kalman}
\alias{temp.kalman}

\title{
	Smooth temperature time series using a Kalman filter/ smoother
}
\description{
	Smoothes a temperature time series uses a Kalman filter/ smoother.
}
\usage{
temp.kalman(wtr, watts, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wtr}{
	Vector (regular time series) of water temperature in degrees C
}
  \item{watts}{
	estimate of watts entering the layer at each time step, from \code{\link{watts.in}}
}
  \item{\dots}{
	parameters to be passed to \code{\link{optim}}
}
}
\details{
	basic model process is x[t] = beta*x[t-1] + c1*watts[t-1]
}
\value{
	 a smoothed temperature time series
}
\references{
	Batt, Ryan D. and Stephen R. Carpenter. 2012. \emph{Free-water lake metabolism: 
	addressing noisy time series with a Kalman filter}. Limnology and Oceanography: Methods
	10: 20-30. doi: 10.4319/lom.2012.10.20
}
\author{
Ryan Batt
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{watts.in}}
\code{\link{metab.kalman}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (wtr, watts, ...) 
{
    guesses <- c(0.9, 1e-04, log(5), log(5))
    fit <- optim(guesses, fn = KFnllTemp, wtr = wtr, watts = watts, 
        ...)
    pars0 <- fit$par
    pars <- c(ar1 = pars0[1], c1 = pars0[2], Q = exp(pars0[3]), 
        H = exp(pars0[4]))
    KFsmoothTemp(pars, wtr = wtr, watts = watts)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
