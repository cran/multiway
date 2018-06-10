\name{mpinv}
\alias{mpinv}
\title{
Moore-Penrose Pseudoinverse
}
\description{
Calculates the Moore-Penrose pseudoinverse of the input matrix using a truncated singular value decomposition.
}
\usage{
mpinv(X, tol = NULL)
}
\arguments{
  \item{X}{
	Real-valued matrix.
}
  \item{tol}{
  Stability tolerance for singular values.
}

}
\note{
	Default tolerance is \code{tol = max(dim(X)) * .Machine$double.eps}.
}
\details{
	Basically returns \code{Y$v \%*\% diag(1/Y$d) \%*\% t(Y$u)} where \code{Y = svd(X)}.
}
\value{
	Returns pseudoinverse of \code{X}.
}
\references{
Moore, E. H. (1920). On the reciprocal of the general algebraic matrix. Bulletin of the American Mathematical Society 26, 394-395.

Penrose, R. (1950). A generalized inverse for matrices. Mathematical Proceedings of the Cambridge Philosophical Society 51, 406-413.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}


\examples{
##########   EXAMPLE   ##########

set.seed(1)
X <- matrix(rnorm(2000),100,20)
Xi <- mpinv(X)
sum( ( X - X \%*\% Xi \%*\% X )^2 )
sum( ( Xi - Xi \%*\% X \%*\% Xi )^2 )
isSymmetric(X \%*\% Xi)
isSymmetric(Xi \%*\% X)

}

\keyword{ array }% use one of  RShowDoc("KEYWORDS")
\keyword{ algebra }% use one of  RShowDoc("KEYWORDS")
