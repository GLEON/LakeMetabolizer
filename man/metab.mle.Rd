\name{metab.mle}
\alias{metab.mle}

\title{
	Metabolism model based on a maximum likelihood parameter estimation framework.
}
\description{
	This function runs the maximum likelihood metabolism model on the supplied 
	gas concentration and other supporting data. This is a common approach
	that allows for the concurrent estimation of metabolism paramters from a timeseries.
}
\usage{
	metab.mle(do.obs, do.sat, k.gas, z.mix, irr, wtr, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{do.obs}{Vector of dissolved oxygen concentration observations, mg L^-1}
  \item{do.sat}{Vector of dissolved oxygen saturation values based on water temperature. Calculate using \code{\link{o2.at.sat}}}
  \item{k.gas}{Vector of kGAS values calculated from any of the gas flux models (e.g., \code{\link{k.cole}}) and converted to kGAS using \code{\link{k600.2.kGAS}}}
  \item{z.mix}{Vector of mixed-layer depths in meters. To calculate, see \code{\link{rLakeAnalyzer}}}
  \item{irr}{Vector of photosynthetically active radiation in umoles/m2/s}
  \item{wtr}{Vector of water temperatures in deg C. Used in scaling respiration with temperature}
  \item{...}{Additional parameters are passed on to the \code{\link{optim}} function}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A named list of parameter estimates.
	
	\item{GPP}{Estimated Gross Primary Productivity}
	\item{R}{Estimated ecosystem respiration}
}
\references{
Solomon CT, DA Bruesewitz, DC Richardson, KC Rose, MC Van de Bogert, PC Hanson, TK Kratz, B Larget, R Adrian, B Leroux Babin,
CY Chiu, DP Hamilton, EE Gaiser, S Hendricks, V Istvanovics, A Laas, DM O’Donnell, ML Pace, E Ryder, PA Staehr, T Torgersen, 
MJ Vanni, KC Weathers, G Zhuw 2013. 
\emph{Ecosystem Respiration: Drivers of Daily Variability and Background Respiration in Lakes around the Globe.}
Limnology and Oceanograph 58 (3): 849–866. doi:10.4319/lo.2013.58.3.0849.
}
\author{
Luke Winslow, Ryan Batt, GLEON Fellows
}

\seealso{
	\code{\link{meta.bayesian}}, \code{\link{metab.bookkeep}}, \code{\link{metab.ols}}.
}
\examples{

}


\keyword{methods}
\keyword{math}
