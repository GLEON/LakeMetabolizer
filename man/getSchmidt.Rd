\name{getSchmidt}
\alias{getSchmidt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Delete this text.%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSchmidt(temperature, gas)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{temperature}{
%%     ~~Describe \code{temperature} here~~
}
  \item{gas}{
%%     ~~Describe \code{gas} here~~
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
function (temperature, gas) 
{
    Schmidt <- data.frame(He = c(377.09, 19.154, 0.50137, 0.005669), 
        O2 = c(1800.6, 120.1, 3.7818, 0.047608), CO2 = c(1911.1, 
            118.11, 3.4527, 0.04132), CH4 = c(1897.8, 114.28, 
            3.2902, 0.039061), SF6 = c(3255.3, 217.13, 6.837, 
            0.08607), N2O = c(2055.6, 137.11, 4.3173, 0.05435), 
        Ar = c(1759.7, 117.37, 3.6959, 0.046527))
    if (!is.character(gas)) {
        stop(paste("gas must be a character. was given as", gas))
    }
    if (length(gas) > 1) {
        stop("only one gas can be specified for this version")
    }
    if (!any(names(Schmidt) == gas)) {
        stop(paste(gas, "not found in list of coded gasses"))
    }
    A <- unlist(Schmidt[gas])[1]
    B <- unlist(Schmidt[gas])[2]
    C <- unlist(Schmidt[gas])[3]
    D <- unlist(Schmidt[gas])[4]
    Sc = as.numeric(A - B * temperature + C * temperature^2 - 
        D * temperature^3)
    return(Sc)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
