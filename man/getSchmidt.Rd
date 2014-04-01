\name{getSchmidt}
\alias{getSchmidt}
\title{
Returns Schmidt number for a specific gas at a given temperature
}
\description{
Schmidt number is temperature dependant, and is the ratio of the 
kinematic viscosity of water to a diffusion coefﬁcient. Coefficients
are included for He, O2, CO2, CH4, SF6, N2O, Ar, and N2.
}

\usage{
getSchmidt(temperature, gas)
}

\arguments{
  \item{temperature}{Numeric vector of water temperatures in deg. Celsius}
  \item{gas}{String for gas code. Valid inputs include: 'He', 'O2', 'CO2', 
    'CH4', 'SF6', 'N2O', 'Ar', and 'N2'}
}

\value{Schmidt number (unitless)}

\references{
  Raymond, Peter A., Christopher J. Zappa, David Butman, Thomas L. Bott, 
  Jody Potter, Patrick Mulholland, Andrew E. Laursen, William H. McDowell,
	and Denis Newbold. \emph{Scaling the gas transfer velocity and hydraulic 
	geometry in streams and small rivers.} Limnology & Oceanography: 
	Fluids & Environments 2 (2012): 41-53.
}
\author{
Jordan Read
}
\note{
Temperature range is only valid from 4-35 deg Celsius
}


\examples{
getSchmidt(temperature=12, gas="O2")
}

\keyword{methods}
\keyword{math}
