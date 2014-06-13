\name{rmv.vars}
\alias{rmv.vars}
\title{subsets data.frame according to header names}
\usage{
  rmv.vars(data, var.name, ignore.missing=TRUE,
    ignore.offset=FALSE)
}
\arguments{
  \item{data}{Object of class data.frame}

  \item{var.name}{A character vector of names to remove
  from \code{data}}

  \item{ignore.missing}{Boolean, should an error be thrown
  if no matching data found}

  \item{ignore.offset}{Should the numerical offset be
  ignored in the match, (e.g. all \code{wtr} columns
  removed, or \code{wtr_0} specifically)}
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

