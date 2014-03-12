\name{calc.zeng}
\alias{calc.zeng}
\title{Estimate sensible and latent heat fluxes}
\description{Returns the sensible and latent heat fluxed based on Zeng et al, 1998}
\usage{
  calc.zeng(dateTime,Ts,airT,Uz,rh,atm.press,wndZ,airTZ,RHZ)
}
\arguments{
  \item{dateTime}{vector of datetime in POSIXct format}
  \item{Ts}{numeric value of surface water temperature, degC}    
  \item{airT}{numeric value of air temperature, degC}  
  \item{Uz}{numeric value of wind speed, m/s}
  \item{RH}{numeric value of relative humidity, \%}
  \item{atm.press}{atmospheric pressure in mb}  
  \item{wndZ}{height of wind measurement, m}
  \item{airTZ}{height of air temperature measurement, m (optional)}
  \item{RHZ}{height of relative humidity measurement, m (optional)}
}

\value{
  A data.frame including sensible and latent heat flux estimates, and 
  other variables used in calculating these fluxes. 
}
\author{
  R. Iestyn. Woolway
}

\seealso{
  \code{\link{k.read}}
  }

\examples{
  dateTime <- as.POSIXct("2013-12-30 23:00")
  Ts <- 22.51
  airT <- 20
  Uz <- 3  
  RH <- 90
  atm.press <- 1013
  wndZ <- 2
  calc.zeng(dateTime,Ts,airT,Uz,RH,atm.press,wndZ)
}

\references{
  Zeng, X., M. Zhao., and Dickinson, R.E. 1998. \emph{Intercomparison of bulk aerodynamic algorithms for the computation of sea surface fluxes using TOGA COARE and TAO data}. Journal of Climate 11: 2628-2644.
}

\keyword{iteration}
\keyword{methods}
\keyword{math}
