\name{krprod}
\alias{krprod}
\title{
Khatri-Rao Product
}
\description{
 Given \code{X} (n-by-p) and \code{Y} (m-by-p), the Khatri-Rao product \code{Z = krprod(X,Y)} is defined as 
 \tabular{c}{
 \code{ Z[,j] = kronecker(X[,j],Y[,j]) }
 }
 which is the mn-by-p matrix containing Kronecker products of corresponding columns of \code{X} and \code{Y}.
}
\usage{
krprod(X,Y)
}
\arguments{
  \item{X}{
	Matrix of order n-by-p.
}
  \item{Y}{
	Matrix of order m-by-p.
}
}
\value{
	The mn-by-p matrix of columnwise Kronecker products.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
	\code{X} and \code{Y} must have the same number of columns.
}

\examples{
##########   EXAMPLE 1   ##########
X <- matrix(1,4,2)
Y <- matrix(1:4,2,2)
krprod(X,Y)


##########   EXAMPLE 2   ##########
X <- matrix(1:2,4,2)
Y <- matrix(1:4,2,2)
krprod(X,Y)


##########   EXAMPLE 3   ##########
X <- matrix(1:2,4,2,byrow=TRUE)
Y <- matrix(1:4,2,2)
krprod(X,Y)

}

