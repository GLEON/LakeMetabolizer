\name{metab.bayesian}
\alias{metab.bayesian}
\title{
	Metabolism model based on a bayesian parameter estimation framework.
}
\description{
	This function runs the bayesian metabolism model on the supplied 
	gas concentration and other supporting data. This allows for both 
	estimates of metabolism along with uncertainty around the parameters.
  }
\usage{
	metab.bayesian(do.obs, do.sat, k.gas, z.mix, date.times, irr, wtr, priors)
}
\arguments{
  \item{do.obs}{Vector of dissovled oxygen concentration observations, mg L^-1}
  \item{do.sat}{Vector of dissolved oxygen saturation values based on water temperature. Calculate using \code{\link{o2.at.sat}}}
  \item{k.gas}{Vector of kGAS values calculated from any of the gas flux models (e.g., \code{\link{k.cole}}) and converted to kGAS using \code{\link{k600.2.kGAS}}}
  \item{z.mix}{Vector of mixed-layer depths in meters. To calculate, see \code{\link{rLakeAnalyzer}}}
  \item{irr}{Vector of photosynthetically active radiation in umoles/m2/s}
  \item{wtr}{Vector of water temperatures in deg C. Used in scaling respiration with temperature}
  \item{priors}{Parameter priors supplied as a named list (example: \code{c("gppMu"=0, "gppSig2"=1E5, "rMu"=0, "rSig2"=1E5, "kSig2"=NA))}}
}

\value{
	A named list of parameter estimates.
	
	\item{GPP}{Estimated Gross Primary Productivity}
	\item{R}{Estimated ecosystem respiration}
}
\author{
	Luke Winslow, Ryan Batt
}

\seealso{
	\code{\link{metab.bootstrap}}, \code{\link{metab.mle}}, \code{\link{meta.bookkeep}}
}

\examples{
  library(rLakeAnalyzer)
  
  \dontrun{
  doobs = load.ts(system.file('extdata', 
			'Sparkling.doobs', package="LakeMetabolizer"))
  wtr = load.ts(system.file('extdata', 
			'Sparkling.wtr', package="LakeMetabolizer"))
  wnd = load.ts(system.file('extdata', 
			'Sparkling.wnd', package="LakeMetabolizer"))
  irr = load.ts(system.file('extdata', 
			'Sparkling.par', package="LakeMetabolizer"))
  
  #Subset a day
  mod.date = as.POSIXct('2009-08-12')
  doobs = doobs[trunc(doobs$datetime, 'day') == mod.date, ]
  wtr = wtr[trunc(wtr$datetime, 'day') == mod.date, ]
  wnd = wnd[trunc(wnd$datetime, 'day') == mod.date, ]
  irr = irr[trunc(irr$datetime, 'day') == mod.date, ]
  
  k600 = k.cole.base(wnd[,2])
  k.gas = k600.2.kGAS.base(k600, wtr[,3], 'O2')
  do.sat = o2.at.sat(wtr[,3], altitude=300)
  
  metab.bayesian(irr=irr[,2], z.mix=rep(1, length(k.gas)), 
            do.sat=do.sat, wtr=wtr[,2],
            k.gas=k.gas, do.obs=doobs[,2])
  }#dontrun
}

\references{
	Holtgrieve, Gordon W., Daniel E. Schindler, Trevor a. Branch, and Z. Teresa A’mar. 2010. 
	\emph{Simultaneous Quantification of Aquatic Ecosystem Metabolism and Reaeration Using a Bayesian 
	Statistical Model of Oxygen Dynamics.} Limnology and Oceanography 55 (3): 1047–1062. 
	doi:10.4319/lo.2010.55.3.1047. http://www.aslo.org/lo/toc/vol_55/issue_3/1047.html.
}

\keyword{methods}
\keyword{math}
