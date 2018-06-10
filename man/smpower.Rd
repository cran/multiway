\name{smpower}
\alias{smpower}
\title{
Symmetric Matrix Power
}
\description{
Raise symmetric matrix to specified power. Default calculates symmetric square root.
}
\usage{
smpower(X, power = 0.5, tol = NULL)
}
\arguments{
  \item{X}{
	Symmetric real-valued matrix.
}
  \item{power}{
	Power to apply to eigenvalues of \code{X}.
}
  \item{tol}{
  Stability tolerance for eigenvalues.
}

}
\note{
	Default tolerance is \code{tol = max(dim(X)) * .Machine$double.eps}.
}
\details{
	Basically returns \code{tcrossprod(Y$vec\%*\%diag(Y$val^power),Y$vec)} where \code{Y = eigen(X,symmetric=TRUE)}.
}
\value{
	Returns \code{X} raised to specified power.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

\examples{
##########   EXAMPLE   ##########

X <- crossprod(matrix(rnorm(2000),100,20))
Xsqrt <- smpower(X)         # square root
Xinv <- smpower(X,-1)       # inverse
Xisqrt <- smpower(X,-0.5)   # inverse square root

}

\keyword{ array }% use one of  RShowDoc("KEYWORDS")
\keyword{ algebra }% use one of  RShowDoc("KEYWORDS")
