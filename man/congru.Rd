\name{congru}
\alias{congru}
\title{
Tucker's Congruence Coefficient
}
\description{
Calculates Tucker's congruence coefficient (uncentered correlation) between \code{x} and \code{y} if these are vectors. If \code{x} and \code{y} are matrices then the congruence between the columns of \code{x} and \code{y} are computed.
}
\usage{
congru(x, y = NULL)
}
\arguments{
  \item{x}{
Numeric vector, matrix or data frame.
}
  \item{y}{
NULL (default) or a vector, matrix or data frame with compatible dimensions to \code{x}. The default is equivalent to \code{y = x} (but more efficient).
}
}
\value{
Returns a scalar or matrix with congruence coefficient(s).
}
\references{
Tucker, L.R. (1951). \emph{A method for synthesis of factor analysis studies} (Personnel Research Section Report No. 984). Washington, DC: Department of the Army.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\details{
Tucker's congruence coefficient is defined as \deqn{r = \frac{\sum_{i=1}^{n}x_{i}y_{i}}{ \sqrt{\sum_{i=1}^{n}x_{i}^{2}\sum_{i=1}^{n}y_{i}^{2}} } }
where \eqn{x_{i}} and \eqn{y_{i}} denote the \eqn{i}-th elements of \code{x} and \code{y}.
}
\note{
If \code{x} is a vector, you must also enter \code{y}.
}
\examples{
##########   EXAMPLE 1   ##########

set.seed(1)
A <- rnorm(100)
B <- rnorm(100)
C <- A*5
D <- A*(-0.5)
congru(A,B)
congru(A,C)
congru(A,D)



##########   EXAMPLE 2   ##########

set.seed(1)
A <- cbind(rnorm(20),rnorm(20))
B <- cbind(A[,1]*-0.5,rnorm(20))
congru(A)
congru(A,B)

}

\keyword{ array }% use one of  RShowDoc("KEYWORDS")
\keyword{ algebra }% use one of  RShowDoc("KEYWORDS")
