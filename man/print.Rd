\name{print}
\alias{print.indscal}
\alias{print.parafac}
\alias{print.parafac2}
\alias{print.sca}
\alias{print.tucker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print Multi-Way Model
}
\description{
Prints constraint, fit, and convergence details for a fit \code{multiway} model.
}
\usage{
\method{print}{indscal}(x,...)
\method{print}{parafac}(x,...)
\method{print}{parafac2}(x,...)
\method{print}{sca}(x,...)
\method{print}{tucker}(x,...)
}
\arguments{
  \item{x}{
  	Object of class "indscal" (output from \code{\link{indscal}}), class "parafac" (output from \code{\link{parafac}}), class "parafac2" (output from \code{\link{parafac2}}), class "sca" (output from \code{\link{sca}}), or class "tucker" (output from \code{\link{tucker}}).
}
  \item{\ldots}{
  Ignored.
}
}
\details{
	See \code{\link{indscal}}, \code{\link{parafac}}, \code{\link{parafac2}}, \code{\link{sca}}, and \code{\link{tucker}} for examples.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\examples{

### see examples for indscal, parafac, parafac2, sca, and tucker

}

