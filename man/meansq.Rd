\name{meansq}
\alias{meansq}
\title{
Mean Square of Given Object
}
\description{
 Calculates the mean square of \code{X}.
}
\usage{
meansq(X, na.rm = FALSE)
}
\arguments{
  \item{X}{
	Numeric scalar, vector, list, matrix, or array.
} 
  \item{na.rm}{
	logical. Should missing values (including \code{NaN}) be removed?
}
}
\value{
	Mean square of \code{X}.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}

\examples{
##########   EXAMPLE 1   ##########
X <- 10
meansq(X)


##########   EXAMPLE 2   ##########
X <- 1:10
meansq(X)


##########   EXAMPLE 3   ##########
X <- matrix(1:10,5,2)
meansq(X)


##########   EXAMPLE 4   ##########
X <- array(matrix(1:10,5,2),dim=c(5,2,2))
meansq(X)


##########   EXAMPLE 5   ##########
X <- vector("list",5)
for(k in 1:5){ X[[k]] <- matrix(1:10,5,2) }
meansq(X)

}

