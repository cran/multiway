\name{print}
\alias{print.cpd}
\alias{print.indscal}
\alias{print.mcr}
\alias{print.parafac}
\alias{print.parafac2}
\alias{print.sca}
\alias{print.tucker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print Multi-Way Model Results
}
\description{
Prints constraint, fit, and convergence details for a fit \code{multiway} model.
}
\usage{
\method{print}{cpd}(x,...)
\method{print}{indscal}(x,...)
\method{print}{mcr}(x,...)
\method{print}{parafac}(x,...)
\method{print}{parafac2}(x,...)
\method{print}{sca}(x,...)
\method{print}{tucker}(x,...)
}
\arguments{
  \item{x}{
  	Object of class "cpd" (output from \code{\link{cpd}}), class "indscal" (output from \code{\link{indscal}}), class "mcr" (output from \code{\link{mcr}}), class "parafac" (output from \code{\link{parafac}}), class "parafac2" (output from \code{\link{parafac2}}), class "sca" (output from \code{\link{sca}}), or class "tucker" (output from \code{\link{tucker}}).
}
  \item{\ldots}{
  Ignored.
}
}
\details{
	See \code{\link{cpd}}, \code{\link{indscal}}, \code{\link{mcr}}, \code{\link{parafac}}, \code{\link{parafac2}}, \code{\link{sca}}, and \code{\link{tucker}} for examples.
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

