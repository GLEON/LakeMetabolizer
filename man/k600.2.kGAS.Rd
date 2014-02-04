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
  kO2 <- k600.2.kGAS(k600=2.4,temperature=20.4,gas='O2')
}

\references{

}

\keyword{GLEON}
