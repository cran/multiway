\name{fitted}
\alias{fitted.indscal}
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
\method{fitted}{indscal}(object,...)
\method{fitted}{parafac}(object,...)
\method{fitted}{parafac2}(object,simplify=TRUE,...)
\method{fitted}{sca}(object,...)
\method{fitted}{tucker}(object,...)
}
\arguments{
  \item{object}{
  	Object of class "indscal" (output from \code{\link{indscal}}), class "parafac" (output from \code{\link{parafac}}), class "parafac2" (output from \code{\link{parafac2}}), class "sca" (output from \code{\link{sca}}), or class "tucker" (output from \code{\link{tucker}}).
}
  \item{simplify}{
  For "parafac2", setting \code{simplify=FALSE} will always return a list of fitted arrays. Default of \code{simplify=TRUE} returns a fitted array if all levels of the nesting mode have the same number of observations (and a list of fitted arrays otherwise).
}
  \item{\ldots}{
  Ignored.
}
}
\value{
  "indscal" objects: 3-way array.
  
  "parafac" objects: 3-way or 4-way array.
  
  "parafac2" objects: 3-way or 4-way array (if possible and simplify=TRUE); otherwise list of 2-way or 3-way arrays.
  
  "sca" objects: list of 2-way arrays.
  
  "tucker" objects: 3-way or 4-way array.
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

