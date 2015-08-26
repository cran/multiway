\name{resign}
\alias{resign}
\alias{resign.indscal}
\alias{resign.parafac}
\alias{resign.parafac2}
\alias{resign.sca}
\alias{resign.tucker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Resigns Multi-Way Factors
}
\description{
Resigns factors from a \code{multiway} object.
}
\usage{
\method{resign}{indscal}(x, mode="B", newsign=1, ...)
\method{resign}{parafac}(x, mode="A", newsign=1, absorb="C", ...)
\method{resign}{parafac2}(x, mode="A", newsign=1, absorb="C", method="pearson", ...)
\method{resign}{sca}(x, mode="B", newsign=1, ...)
\method{resign}{tucker}(x, mode="A",newsign=1, ...)
}
\arguments{
  \item{x}{
  	Object of class "indscal" (output from \code{\link{indscal}}), class "parafac" (output from \code{\link{parafac}}), class "parafac2" (output from \code{\link{parafac2}}), class "sca" (output from \code{\link{sca}}), or class "tucker" (output from \code{\link{tucker}}).
}
  \item{mode}{
  Character indicating which mode to resign.
}
  \item{newsign}{
  Desired sign for each column of resigned mode. Can input a scalar or a vector with length equal to the number of factors for the given mode. If \code{x} is of class "parafac2" and \code{mode="A"} you can input a list of covariates (see Details).
}
  \item{absorb}{
  Character indicating which mode should absorb the inverse of the rescalings applied to \code{mode} (cannot be equal to \code{mode}).
}
  \item{method}{
  Correlation method to use if \code{newsign} is a list input (see Details).
}
  \item{...}{
  Ignored.
}
}
\value{
  Same as input.
}
\details{
	If \code{x} is of class "parafac2" and \code{mode="A"}, the input \code{newsign} can be a list where each element contains a covariate vector for resigning Mode A. You need \code{length(newsign[[k]]) = nrow(x$A$H[[k]])} for all \code{k} when \code{newsign} is a list. In this case, the resigning is implemented according to the sign of \code{cor(newsign[[k]], x$A$H[[k]]\%*\%x$A$G[,1], method)}. See Helwig (2013) for details.
  
  See \code{\link{indscal}}, \code{\link{parafac}}, \code{\link{parafac2}}, \code{\link{sca}}, and \code{\link{tucker}} for more details.
}
\references{
Helwig, N. E. (2013). The special sign indeterminacy of the direct-fitting Parafac2 model: Some implications, cautions, and recommendations, for Simultaneous Component Analysis. \emph{Psychometrika, 78}, 725-739.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\examples{

### see examples for indscal, parafac, parafac2, sca, and tucker

}

