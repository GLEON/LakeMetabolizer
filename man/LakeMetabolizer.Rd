\name{LakeMetabolizer-package}
\alias{LakeMetabolizer-package}
\alias{LakeMetabolizer}
\docType{package}
\title{
Tools for the analysis of ecosystem metabolism
}
\description{
A collection of tools for the calculation of freewater metabolism. 

}
\details{
\tabular{ll}{
Package: \tab LakeMetabolizer\cr
Version: \tab 1.0.1\cr
License: \tab file LICENSE\cr
Imports: \tab 
rLakeAnalyzer, rjags, R2jags, R2WinBUGS\cr
Depends: \tab 
R (>= 2.15.0)\cr
Packaged: \tab 2013-08-01 09:09:24 UTC; Luke\cr
Repository: \tab http://github.com/GLEON/LakeMetabolizer\cr
Date/Publication: \tab 2014-01-01 09:09:24\cr
Built: \tab R 3.0.2; ; 2014-03-18 21:17:12 UTC; windows\cr
}

Index:
\preformatted{
calc.lw.net             Estimate net long wave heat radiation
calc.zeng               Estimate sensible and latent heat fluxes
getSchmidt              Delete this text.
k.cole                  Returns the gas exchange velocity for k600
                        based on Cole and Caraco, 1998
k.crusius               Returns k600 from Crusius and Wanninkhof 2003
k.macIntyre             k600 from MacIntyre et al. 2010
k.read                  Returns the gas exchange velocity for O2 based
                        on Read et al, 2012
k600.2.kGAS             Returns the gas exchange velocity for gas of
                        interest w/ no unit conversions
manuscript.code         replace this text.
metab.bayes             Returns metabolism estimates based on a
                        bayesian model based on XXXX
metab.bookkeep          Returns metabolism estimates for observed
                        system based on Cole et al, 2010
metab.optim             replace this text.
o2.at.sat               Calculates the equilibrium saturation
                        concentration of oxygen in water at the
                        supplied conditions
scale.exp.wind          replace this text.
sun.rise.set            replace this text.
}

}
\author{
Luke Winslow, Jake Zwart, Ryan Batt, Jessica Corman, Hilary Dugan, Paul Hanson, Aline Jaimes, Jordan Read, Richard Woolway

Maintainer: Luke Winslow <lawinslow@gmail.com>

}
\references{

}

\keyword{ package }
\seealso{
 \code{\link[rLakeAnalyzer:rLakeAnalyzer-package]{rLakeAnalyzer}} 
}
\examples{

}
