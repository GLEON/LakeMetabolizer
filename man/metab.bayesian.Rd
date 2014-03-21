\name{metab.bayes}
\alias{metab.bayes}
\title{
Returns metabolism estimates based on a bayesian model based on XXXX
}
\description{
  This function runs the bayesian metabolism model on the supplied 
  gas concentration and other supporting data. This allows for both 
  estimates of metabolism along with uncertainty around the parameters.
  }
\usage{
  metab.bayes(
}
\arguments{
  \item{dateTime}{datetime (Y-\%m-\%d \%H:\%M), POSIXct format}
}

\value{
  Matrix of four values. 
}
\author{
  Luke Winslow, Ryan Batt
}

\seealso{
  \code{\link{metab.bootstrap}}, \code{\link{metab.optim}}, \code{\link{meta.bookkeeping}}
}

\examples{
  library(rLakeAnalyzer)
  doobs = load.ts('inst//extdata/Sparkling.doobs')
  wtr = load.ts('inst//extdata/Sparkling.wtr')
  wnd = load.ts('inst//extdata/Sparkling.wnd')
  irr = load.ts('inst/extdata/Sparkling.par')
  
  #Subset a day
  mod.date = as.POSIXct('2009-08-12')
  doobs = doobs[trunc(doobs$datetime, 'day') == mod.date, ]
  wtr = wtr[trunc(wtr$datetime, 'day') == mod.date, ]
  wnd = wnd[trunc(wnd$datetime, 'day') == mod.date, ]
  irr = irr[trunc(irr$datetime, 'day') == mod.date, ]
  
  k.gas = k600.2.kGAS(k.cole(wnd[,2]), wtr[,3], 'O2')
  do.sat = o2.at.sat(wtr[,3], altitude=300)
  
  metab.bayes(irr=irr[,2], z.mix=rep(1, length(k.gas)), 
            do.sat=do.sat, wtr=wtr[,2],
            k.gas=k.gas, do.obs=doobs[,2])
  
  
}

\references{
  Read, J.S., Hamilton, D.P., Desai, A.R., Rose, K.C., MacIntyre, S.,
  Lenters, J.D., Smyth, R.L., Hanson, P. .C., Cole, J.J., Staehr, P.A.,
  Rusak, J.A., Pierson, D.C., Brookes, J.D., Laas, A., Wu, C.H. 2012.
  \emph{Lake-size dependency of wind shear and convection as controls on gas
  exchange}. Geophysical Research Letters 39, DOI: 10.1029/2012GL051886.
}

\keyword{GLEON}
