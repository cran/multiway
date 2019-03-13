\name{reorder}
\alias{reorder.cpd}
\alias{reorder.indscal}
\alias{reorder.mcr}
\alias{reorder.parafac}
\alias{reorder.parafac2}
\alias{reorder.sca}
\alias{reorder.tucker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reorder Multi-Way Factors
}
\description{
Reorders factors from a \code{multiway} object.
}
\usage{
\method{reorder}{cpd}(x, neworder, ...)
\method{reorder}{indscal}(x, neworder, ...)
\method{reorder}{mcr}(x, neworder, mode = "A", ...)
\method{reorder}{parafac}(x, neworder, ...)
\method{reorder}{parafac2}(x, neworder, ...)
\method{reorder}{sca}(x, neworder, ...)
\method{reorder}{tucker}(x, neworder, mode = "A", ...)
}
\arguments{
  \item{x}{
  	Object of class "cpd" (output from \code{\link{cpd}}), "indscal" (output from \code{\link{indscal}}), class "mcr" (output from \code{\link{mcr}}), class "parafac" (output from \code{\link{parafac}}), class "parafac2" (output from \code{\link{parafac2}}), class "sca" (output from \code{\link{sca}}), or class "tucker" (output from \code{\link{tucker}}).
}
  \item{neworder}{
  Vector specifying the new factor ordering. Must be a permutation of the integers 1 to \code{nfac}.
}
  \item{mode}{
  Character indicating which mode to reorder (only for \code{tucker} models). For 3-way Tucker options include "A", "B", and "C". For 4-way Tucker, options are "A", "B", "C", and "D". 
}
  \item{...}{
  Ignored.
}
}
\value{
  Same as input.
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

