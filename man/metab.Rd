\name{metab}
\alias{metab}
\title{Calculate Metabolism}
\usage{
  metab(data, method)
}
\arguments{
  \item{data}{a data.frame whose columns are "year", "doy",
  "datetime", "do.obs","do.sat","k.gas","z.mix", "irr",
  "wtr", "priors". Data columns (i.e., not year, doy, or
  datetime) that are not used by a particular statistical
  method do not need to be supplied.}

  \item{method}{a character string specifying one of the 5
  statistical methods ("bayesian", "bookkeep", "kalman",
  "ols", "mle")}
}
\description{
  Returns daily time series of gross primary production
  (GPP), respiration (R), and net ecosystem production
  (NEP). Depending on the method used, other information
  may be returned as well. Calculations are made using one
  of 5 statistical methods.
}

