\name{metab.kalman}
\alias{metab.kalman}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Estimate metabolism using a Kalman filter
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
metab.kalman(do.obs, do.sat, k.gas, z.mix, irr, wtr, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{do.obs}{
	numeric vector of dissolved oxygen concentration, mg/ L
}
  \item{do.sat}{
	numeric vector of the saturated concentration of dissolved oxygen given current temperature and pressure (\code{\link{o2.at.sat}})
}
  \item{k.gas}{
	Vector of kGAS values calculated from any of the gas flux models (e.g., \code{\link{k.cole}}) and converted to kGAS using \code{\link{k600.2.kGAS}}
}
  \item{z.mix}{
	Vector of mixed-layer depths in meters. To calculate, see \code{\link{rLakeAnalyzer}}
}
  \item{irr}{
	Vector of photosynthetically active radiation in umoles/m2/s
}
  \item{wtr}{
	Vector of water temperatures in deg C. Used in scaling respiration with temperature
}
  \item{\dots}{
	arguments to be passed to \code{\link{optim}}
}
}
\details{
	Uses a Kalman filter to fit parameters relating irr to GPP, log(wtr) to R, process error, and observation error. Also provides a smoothed time series of oxygen concentration.
}
\value{
\item{smoothDO }{smoothed time series of oxygen concentration, from Kalman smoother}
\item{params }{parameters estimated by the Kalman filter}
\item{metab }{daily metabolism estimates in mg O2 / L / day}

}
\references{
	Batt, Ryan D. and Stephen R. Carpenter. 2012. \emph{Free-water lake metabolism: 
	addressing noisy time series with a Kalman filter}. Limnology and Oceanography: Methods
	10: 20-30. doi: 10.4319/lom.2012.10.20
}
\author{
	Ryan Batt, Luke Winslow
}
\note{

\seealso{
\code{\link{temp.smooth}}
\code{\link{watts.in}}
\code{\link{metab}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  \dontrun{
  doobs <- load.ts(system.file('extdata', 
			'Sparkling.doobs', package="LakeMetabolizer"))
  wtr <- load.ts(system.file('extdata', 
			'Sparkling.wtr', package="LakeMetabolizer"))
  wnd <- load.ts(system.file('extdata', 
			'Sparkling.wnd', package="LakeMetabolizer"))
  irr <- load.ts(system.file('extdata', 
			'Sparkling.par', package="LakeMetabolizer"))

  #Subset a day
  Sys.setenv(TZ='GMT')
  mod.date <- as.POSIXct('2009-08-12')
  doobs <- doobs[trunc(doobs$datetime, 'day') == mod.date, ]
  wtr <- wtr[trunc(wtr$datetime, 'day') == mod.date, ]
  wnd <- wnd[trunc(wnd$datetime, 'day') == mod.date, ]
  irr <- irr[trunc(irr$datetime, 'day') == mod.date, ]

  k600 <- k.cole.base(wnd[,2])
  k.gas <- k600.2.kGAS.base(k600, wtr[,3], 'O2')
  do.sat <- o2.at.sat.base(wtr[,3], altitude=300)

  metab.kalman(irr=irr[,2], z.mix=rep(1, length(k.gas)), 
            do.sat=do.sat, wtr=wtr[,2],
            k.gas=k.gas, do.obs=doobs[,2])
  }#dontrun

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
