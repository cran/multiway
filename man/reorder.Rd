\name{reorder}
\alias{reorder.indscal}
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
\method{reorder}{indscal}(x, neworder, ...)
\method{reorder}{parafac}(x, neworder, ...)
\method{reorder}{parafac2}(x, neworder, ...)
\method{reorder}{sca}(x, neworder, ...)
\method{reorder}{tucker}(x, neworder, mode="A", ...)
}
\arguments{
  \item{x}{
  	Object of class "indscal" (output from \code{\link{indscal}}), class "parafac" (output from \code{\link{parafac}}), class "parafac2" (output from \code{\link{parafac2}}), class "sca" (output from \code{\link{sca}}), or class "tucker" (output from \code{\link{tucker}}).
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
	See \code{\link{indscal}}, \code{\link{parafac}}, \code{\link{parafac2}}, \code{\link{sca}}, and \code{\link{tucker}} for more details.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\examples{

### see examples for indscal, parafac, parafac2, sca, and tucker

}

