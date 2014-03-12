\name{metab.optim}
\alias{metab.optim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
metab.optim(do.obs, irr, do.sat, z.mix, k.gas, timestep)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{do.obs}{
%%     ~~Describe \code{do.obs} here~~
}
  \item{irr}{
%%     ~~Describe \code{irr} here~~
}
  \item{do.sat}{
%%     ~~Describe \code{do.sat} here~~
}
  \item{z.mix}{
%%     ~~Describe \code{z.mix} here~~
}
  \item{k.gas}{
%%     ~~Describe \code{k.gas} here~~
}
  \item{timestep}{
%%     ~~Describe \code{timestep} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (do.obs, irr, do.sat, z.mix, k.gas, timestep) 
{
    n.obs = length(do.obs)
    if (n.obs != length(irr) || n.obs != length(do.sat) || n.obs != 
        length(z.mix) || n.obs != length(k.gas)) {
        stop("All input data to metab.optim must be the same length")
    }
    if (n.obs < 4) {
        stop("You must supply at least 4 observations to metab.optim.")
    }
    to.optim <- function(iotaRhoDoinit) {
        calc.do.nll(iotaRhoDoinit[1], iotaRhoDoinit[2], iotaRhoDoinit[3], 
            irr, do.sat, z.mix, k.gas, do.obs)
    }
    rho = range(do.obs)/n.obs
    iota = rho
    doInit = do.obs[1]
    optimOut = optim(par = c(iota, rho, doInit), fn = to.optim)
    iota.hat <- optimOut$par[1] * (1440/timestep)
    rho.hat <- optimOut$par[2] * (1440/timestep)
    do.init.hat <- optimOut$par[3]
    convergence <- optimOut$convergence
    nll.val <- optimOut$value
    gpp.hat <- (iota.hat/(60 * 60 * 24)) * sum(irr * timestep * 
        60)
    return(list(iota = iota.hat, rho = rho.hat, gpp = gpp.hat, 
        mod.info = list(do.init = do.init.hat, convergence = convergence, 
            nll = nll.val)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
