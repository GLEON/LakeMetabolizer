\name{k600.2.kGAS}
\alias{k600.2.kGAS}
\title{Returns the gas exchange velocity for gas of interest w/ no unit conversions}
\description{Returns the gas exchange velocity for gas of interest w/ no unit conversions}
\usage{
  k600.2.kGAS <- function(k600,temperature,gas)
}
\arguments{
  \item{k600}{k600 as vector array of numbers or single number}
  \item{temperature}{Water temperature as vector array of numbers or single number}
  \item{gas}{gas for conversion, as string (e.g., 'CO2' or 'O2')}}

\value{
  Numeric value of gas exchange velocity for gas
}
\author{
  Jordan S Read}

\seealso{
  \code{\link{k.cole}}, \code{\link{k.crusius}}, \code{\link{k.read}}, \code{\link{k.macIntyre}}
}

\examples{
	## single example
  kO2 <- k600.2.kGAS.base(k600=2.4,temperature=20.4,gas='O2')
  
  ## Timeseries example
	#load data
	data.path = system.file('extdata', '', package="LakeMetabolizer")
  sp.data = load.all.data('sparkling', data.path)
  ts.data = sp.data$data #pull out just the timeseries data
  
  
  #calculate U10 and add it back onto the original 
  u10 = scale.exp.wind(ts.data)
  ts.data = rmv.var(ts.data, 'wnd', ignore.offset=TRUE) #drop old wind speed column
  ts.data = merge(ts.data, u10)                          #merge new u10 into big dataset  
	
	k600 = k.cole(ts.data)
	ts.data = merge(k600, ts.data)
	
  k.gas = k600.2.kGAS(ts.data, 'O2')
}

\references{

}

\keyword{GLEON}
