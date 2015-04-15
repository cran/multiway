\name{smpower}
\alias{smpower}
\title{
Symmetric Matrix Power
}
\description{
Raise symmetric matrix to specified power. Default calculates symmetric square-root.
}
\usage{
smpower(X,power=0.5,tol=.Machine$double.eps)
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
\value{
	Returns \code{X} raised to specified power.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
	Basically returns \code{tcrossprod(Y$vec\%*\%diag(Y$val^power),Y$vec)} where \code{Y = eigen(X,symmetric=TRUE)}.
}

\examples{
##########   EXAMPLE   ##########

X <- crossprod(matrix(rnorm(2000),100,20))
Xsqrt <- smpower(X)         # square root
Xinv <- smpower(X,-1)       # inverse
Xisqrt <- smpower(X,-0.5)   # inverse square root

}

