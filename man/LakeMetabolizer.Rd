% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LakeMetabolizer.R
\docType{package}
\name{LakeMetabolizer}
\alias{LakeMetabolizer}
\alias{LakeMetabolizer-package}
\title{Tools for the analysis of ecosystem metabolism}
\description{
A collection of tools for the calculation of freewater metabolism.
}
\details{
\tabular{ll}{
Package: \tab LakeMetabolizer\cr
License: \tab GPL (>=2 )\cr
Repository: \tab http://github.com/GLEON/LakeMetabolizer\cr
}

LakeMetabolizer contains sets of functions to handle metabolism modeling
through the free water oxygen technique.

\strong{Metabolism models}

\link{metab.bayesian}\cr
\link{metab.bookkeep}\cr
\link{metab.kalman}\cr
\link{metab.mle}\cr
\link{metab.ols}\cr
\link{metab}\cr

\strong{Gas exchange models}

\link{k.cole}\cr
\link{k.crusius}\cr
\link{k.macIntyre}\cr
\link{k.read}\cr
\link{k.heiskanen}\cr
\link{k.read.soloviev}\cr
\link{k600.2.kGAS}\cr


\strong{Various Helper functions for light and gas}

\link{calc.lw.net}\cr
\link{calc.zeng}\cr
\link{getSchmidt}\cr
\link{o2.at.sat}\cr
\link{wind.scale}\cr
\link{sun.rise.set}\cr
}
\author{
Luke Winslow, Jake Zwart, Ryan Batt, Jessica Corman, Hilary Dugan, Paul Hanson, Aline Jaimes, Jordan Read, Richard Woolway

Maintainer: Luke Winslow <lawinslow@gmail.com>
}
\references{
The original concept of freewater metabolism attributed to:

Odum, Howard T. 1956. \emph{Primary Production in Flowing Waters.}
Limnology and Oceanography 1 (2): 102-117. http://www.jstor.org/stable/2833008.
}

