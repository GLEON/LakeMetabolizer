\name{rmv.var}
\alias{rmv.var}
\title{subsets data.frame according to header names}
\usage{
  rmv.vars(data, var.names)
}
\arguments{
  \item{\code{data}}{Object of class data.frame}

  \item{\code{var.names}}{A character vector of names to
  remove from \code{data}}
}
\value{
  An object of class data.frame
}
\description{
  subsets \code{data} according to header names. Excludes
  all matches to \code{var.name}
}
\author{
  Luke A. Winslow
}
\seealso{
  \link{has.vars} \link{get.vars}
}
\keyword{methods}

