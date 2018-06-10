\name{fnnls}
\alias{fnnls}
\title{
Fast Non-Negative Least Squares
}
\description{
 Finds the vector \code{b} minimizing 
 \tabular{c}{
 \code{sum( (y - X \%*\% b)^2 )}
 }
 subject to \code{b[j] >= 0} for all \code{j}.
}
\usage{
fnnls(XtX, Xty, ntol = NULL)
}
\arguments{
  \item{XtX}{
	Crossproduct matrix \code{crossprod(X)} of dimension p-by-p.
} 
  \item{Xty}{
  Crossproduct vector \code{crossprod(X,y)} of length p-by-1.
}
  \item{ntol}{
  Tolerance for non-negativity.
}
}
\value{
	The vector \code{b} such that \code{b[j] >= 0} for all \code{j}.
}
\references{
Bro, R., & De Jong, S. (1997). A fast non-negativity-constrained least squares algorithm. \emph{Journal of Chemometrics, 11}, 393-401.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
  Default non-negativity tolerance: \code{ntol=10*(.Machine$double.eps)*max(colSums(abs(XtX)))*p}.
}

\examples{
##########   EXAMPLE 1   ##########
X <- matrix(1:100,50,2)
y <- matrix(101:150,50,1)
beta <- solve(crossprod(X))\%*\%crossprod(X,y)
beta
beta <- fnnls(crossprod(X),crossprod(X,y))
beta


##########   EXAMPLE 2   ##########
X <- cbind(-(1:50),51:100)
y <- matrix(101:150,50,1)
beta <- solve(crossprod(X))\%*\%crossprod(X,y)
beta
beta <- fnnls(crossprod(X),crossprod(X,y))
beta


##########   EXAMPLE 3   ##########
X <- matrix(rnorm(400),100,4)
btrue <- c(1,2,0,7)
y <- X\%*\%btrue + rnorm(100)
fnnls(crossprod(X),crossprod(X,y))


##########   EXAMPLE 4   ##########
X <- matrix(rnorm(2000),100,20)
btrue <- runif(20)
y <- X\%*\%btrue + rnorm(100)
beta <- fnnls(crossprod(X),crossprod(X,y))
crossprod(btrue-beta)/20

}

\keyword{ optimize }% __ONLY ONE__ keyword per line
