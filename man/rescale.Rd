\name{rescale}
\alias{rescale}
\alias{rescale.cpd}
\alias{rescale.indscal}
\alias{rescale.mcr}
\alias{rescale.parafac}
\alias{rescale.parafac2}
\alias{rescale.sca}
\alias{rescale.tucker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rescales Multi-Way Factors
}
\description{
Rescales factors from a \code{multiway} object.
}
\usage{
\method{rescale}{cpd}(x, mode = 1, newscale = 1, absorb = 3, ...)
\method{rescale}{indscal}(x, mode = "B", newscale = 1, ...)
\method{rescale}{mcr}(x, mode = "A", newscale = 1, absorb = "C", ...)
\method{rescale}{parafac}(x, mode = "A", newscale = 1, absorb = "C", ...)
\method{rescale}{parafac2}(x, mode = "A", newscale = 1, absorb = "C", ...)
\method{rescale}{sca}(x, mode = "B", newscale = 1, ...)
\method{rescale}{tucker}(x, mode = "A", newscale = 1, ...)
}
\arguments{
  \item{x}{
  	Object of class "indscal" (output from \code{\link{indscal}}), class "mcr" (output from \code{\link{mcr}}), class "parafac" (output from \code{\link{parafac}}), class "parafac2" (output from \code{\link{parafac2}}), class "sca" (output from \code{\link{sca}}), or class "tucker" (output from \code{\link{tucker}}).
}
  \item{mode}{
  Character indicating which mode to rescale. For "cpd" objects, should be an integer between 1 and N.
}
  \item{newscale}{
  Desired root mean-square for each column of rescaled mode. Can input a scalar or a vector with length equal to the number of factors for the given mode.
}
  \item{absorb}{
  Character indicating which mode should absorb the inverse of the rescalings applied to \code{mode} (cannot be equal to \code{mode}). For "cpd" objects, should be an integer between 1 and N.
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
\references{
Helwig, N. E. (2013). The special sign indeterminacy of the direct-fitting Parafac2 model: Some implications, cautions, and recommendations, for Simultaneous Component Analysis. \emph{Psychometrika, 78}, 725-739. \doi{10.1007/S11336-013-9331-7}
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

