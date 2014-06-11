\name{watts.in}
\alias{watts.in}
\title{
	Simple estimate of energy gained by a layer of water
}
\description{
	Estimate the amount of energy gained by a layer of water as the difference between 
	energy entering from the top of the layer and energy leaving at the bottom. 
	Energy gained/ lost is calculated from photosynthetically active radiation 
	(PAR, which is then converted to watts) and an estimate of kd (light attenuation coefficient) 
	which is derived from the depth of 1 percent surface light.
}
\usage{
watts.in(top, bot, irr, z1perc)
}

\arguments{
  \item{top}{
	Depth of the top of the layer, in meters
}
  \item{bot}{
	Depth of the bottom of the layer, in meters
}
  \item{irr}{
	PAR in uE/s (umol / m^2 / s)
}
  \item{z1perc}{
	Depth of 1 percent of surface light, in meters
}
}
\details{
	This rough estimate is used in the Kalman filter/ smoother for water temperature.
	It does not account for a variety of potentially important factors, and is made specifically for use
	with temp.kalman(), which uses maximum likelihood to fit a linear coefficient that converts this heat gain estimate
	into temperature change.
}
\value{
	numeric vector of estimates of energy gain 
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
	\code{\link{temp.smooth}}
	\code{\link{metab.kalman}}
}
\examples{
watts.in(3.2, 4, 1200, 4.5)

## The function is currently defined as
function (top, bot, irr, z1perc) 
{
    watts <- 0.2174 * irr
    kd <- log(0.01)/-z1perc
    watts * exp(-kd * top) - watts * exp(-kd * bot)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ temperature }
\keyword{ watts }% __ONLY ONE__ keyword per line
