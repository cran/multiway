\name{fitted}
\alias{fitted.cpd}
\alias{fitted.indscal}
\alias{fitted.mcr}
\alias{fitted.parafac}
\alias{fitted.parafac2}
\alias{fitted.sca}
\alias{fitted.tucker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract Multi-Way Fitted Values
}
\description{
Calculates fitted array (or list of arrays) from a \code{multiway} object.
}
\usage{
\method{fitted}{cpd}(object, ...)
\method{fitted}{indscal}(object, ...)
\method{fitted}{mcr}(object, type = c("X", "Y"), ...)
\method{fitted}{parafac}(object, ...)
\method{fitted}{parafac2}(object, simplify = TRUE, ...)
\method{fitted}{sca}(object, ...)
\method{fitted}{tucker}(object, ...)
}
\arguments{
  \item{object}{
  	Object of class "cpd" (output from \code{\link{cpd}}), "indscal" (output from \code{\link{indscal}}), class "mcr" (output from \code{\link{mcr}}), class "parafac" (output from \code{\link{parafac}}), class "parafac2" (output from \code{\link{parafac2}}), class "sca" (output from \code{\link{sca}}), or class "tucker" (output from \code{\link{tucker}}).
}
  \item{simplify}{
  For "parafac2", setting \code{simplify = FALSE} will always return a list of fitted arrays. Default of \code{simplify = TRUE} returns a fitted array if all levels of the nesting mode have the same number of observations (and a list of fitted arrays otherwise).
}
  \item{type}{
  For "mcr", setting \code{type = "X"} returns the fitted predictor array (default), whereas setting \code{type = "Y"} retuns the fitted response array.
}
  \item{\ldots}{
  Ignored.
}
}
\value{
  "cpd" objects: N-way array.
  
  "indscal" objects: 3-way array.
  
  "mcr" objects: 3-way (X) or 2-way (Y) array.
  
  "parafac" objects: 3-way or 4-way array.
  
  "parafac2" objects: 3-way or 4-way array (if possible and simplify=TRUE); otherwise list of 2-way or 3-way arrays.
  
  "sca" objects: list of 2-way arrays.
  
  "tucker" objects: 3-way or 4-way array.
}
\details{
	See \code{\link{cpd}}, \code{\link{indscal}}, \code{\link{mcr}}, \code{\link{parafac}}, \code{\link{parafac2}}, \code{\link{sca}}, and \code{\link{tucker}} for more details.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\examples{
# See examples for... 
#   cpd (Canonical Polyadic Decomposition)
#   indscal (INividual Differences SCALing)
#   mcr (Multiway Covariates Regression)
#   parafac (Parallel Factor Analysis-1)
#   parafac2 (Parallel Factor Analysis-2)
#   sca (Simultaneous Component Analysis)
#   tucker (Tucker Factor Analysis)
}

