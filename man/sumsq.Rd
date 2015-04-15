\name{sumsq}
\alias{sumsq}
\title{
Sum-of-Squares of Given Object
}
\description{
 Calculates the sum-of-squares of \code{X}.
}
\usage{
sumsq(X)
}
\arguments{
  \item{X}{
	Numeric scalar, vector, list, matrix, or array.
}
}
\value{
	Sum-of-squares of \code{X}.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

\examples{
##########   EXAMPLE 1   ##########
X <- 10
sumsq(X)


##########   EXAMPLE 2   ##########
X <- 1:10
sumsq(X)


##########   EXAMPLE 3   ##########
X <- matrix(1:10,5,2)
sumsq(X)


##########   EXAMPLE 4   ##########
X <- array(matrix(1:10,5,2),dim=c(5,2,2))
sumsq(X)


##########   EXAMPLE 5   ##########
X <- vector("list",5)
for(k in 1:5){ X[[k]] <- matrix(1:10,5,2) }
sumsq(X)

}

