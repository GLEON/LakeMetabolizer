\name{metab.bookkeep}
\alias{metab.bookkeep}
\title{Metabolism model based on simple day/night summation NEP-interpreted changes in DO.}
\description{
	This model is a simple model based on the assumption that movements in DO during the day 
	are due to NEP and gas exchange. Respiration is estimated from night-time decreases. 
	GPP is calculated from the algebraic manipulation of NEP and R. Based on Cole et al 2000.
}
\usage{
	metab.bookkeep(do.obs, do.sat, k.gas, z.mix, irr, ...)
}
\arguments{
  \item{do.obs}{Vector of dissovled oxygen concentration observations, mg L^-1}
  \item{do.sat}{Vector of dissolved oxygen saturation values based on water temperature. Calculate using \code{\link{o2.at.sat}}}
  \item{k.gas}{Vector of kGAS values calculated from any of the gas flux models (e.g., \code{\link{k.cole}}) and converted to kGAS using \code{\link{k600.2.kGAS}}}
  \item{z.mix}{Vector of mixed-layer depths in meters. To calculate, see \code{\link{ts.meta.depths}}}
  \item{irr}{Vector of integers: 1 for day, 0 for night}
  \item{...}{Not used, extra input is ignored.}
}

\value{
	A named list of parameter estimates.
	
	\item{GPP}{Estimated Gross Primary Productivity}
	\item{R}{Estimated ecosystem respiration}
}
\author{
  R. Iestyn. Woolway, Hilary Dugan, Luke Winslow, Jordan S Read, GLEON fellows
}

\seealso{
  \code{\link{metab.bayesian}}, \code{\link{metab.mle}}, \code{\link{metab.bookkeep}}
}

\examples{
  library(rLakeAnalyzer)
  Sys.setenv(TZ='GMT')
  
  doobs = load.ts(system.file('extdata', 
			'Sparkling.doobs', package="LakeMetabolizer"))
  wtr = load.ts(system.file('extdata', 
			'Sparkling.wtr', package="LakeMetabolizer"))
  wnd = load.ts(system.file('extdata', 
			'Sparkling.wnd', package="LakeMetabolizer"))
  
  #Subset a day
  mod.date = as.POSIXct('2009-07-12', 'GMT')
  doobs = doobs[trunc(doobs$datetime, 'day') == mod.date, ]
  wtr = wtr[trunc(wtr$datetime, 'day') == mod.date, ]
  wnd = wnd[trunc(wnd$datetime, 'day') == mod.date, ]
  
  k.gas = k600.2.kGAS.base(k.cole.base(wnd[,2]), wtr[,3], 'O2')
  do.sat = o2.at.sat.base(wtr[,3], altitude=300)
  
  # Must supply 1 for daytime timesteps and 0 for nighttime timesteps
  irr = as.integer(is.day(doobs[,1], 45))
  
  metab.bookkeep(doobs[,2], do.sat, k.gas, z.mix=1, irr)
  
}

\references{
	Cole, Jonathan J., Michael L. Pace, Stephen R. Carpenter, and James F. Kitchell. 2000. 
	\emph{Persistence of Net Heterotrophy in Lakes during Nutrient Addition and Food Web Manipulations.} 
	Limnology and Oceanography 45 (8): 1718?1730. doi:10.4319/lo.2000.45.8.1718.
}

\keyword{methods}
\keyword{math}
