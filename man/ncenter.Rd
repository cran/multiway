\name{ncenter}
\alias{ncenter}
\title{
Center n-th Dimension of Array
}
\description{
 Fiber-center across the levels of the specified mode. Can input 2-way, 3-way, and 4-way arrays, or input a list containing array elements.
 
 With \code{X} a matrix (I-by-J) there are two options:
 \tabular{rcl}{
 \code{ nmode=1: } \tab \tab \code{ x[i,j] - mean(x[,j]) } \cr
 \code{ nmode=2: } \tab \tab \code{ x[i,j] - mean(x[i,]) } \cr
 }
 
 With \code{X} a 3-way array (I-by-J-by-K) there are three options:
 \tabular{rcl}{
 \code{ nmode=1: } \tab \tab \code{ x[i,j,k] - mean(x[,j,k]) } \cr
 \code{ nmode=2: } \tab \tab \code{ x[i,j,k] - mean(x[i,,k]) } \cr
 \code{ nmode=3: } \tab \tab \code{ x[i,j,k] - mean(x[i,j,]) } \cr
 }
 
 With \code{X} a 4-way array (I-by-J-by-K-by-L) there are four options:
 \tabular{rcl}{
 \code{ nmode=1: } \tab \tab \code{ x[i,j,k,l] - mean(x[,j,k,l]) } \cr
 \code{ nmode=2: } \tab \tab \code{ x[i,j,k,l] - mean(x[i,,k,l]) } \cr
 \code{ nmode=3: } \tab \tab \code{ x[i,j,k,l] - mean(x[i,j,,l]) } \cr
 \code{ nmode=4: } \tab \tab \code{ x[i,j,k,l] - mean(x[i,j,k,]) } \cr
 }

}
\usage{
ncenter(X,nmode=1)
}
\arguments{
  \item{X}{
	Array (2-way, 3-way, or 4-way) or a list containing array elements.
}
  \item{nmode}{
	Mode to center across.
}
}
\value{
	Returns centered version of \code{X}.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
	When entering a list with array elements, each element must (a) be an array (2-way, 3-way, or 4-way), and (b) have the same dimension (i.e, array size).
}

\examples{
##########   EXAMPLE 1   ##########
X <- matrix(rnorm(2000),100,20)
Xc <- ncenter(X)          # center across rows
sum(colSums(Xc))
Xc <- ncenter(Xc,nmode=2) # recenter across columns
sum(colSums(Xc)) 
sum(rowSums(Xc))


##########   EXAMPLE 2   ##########
X <- array(rnorm(20000),dim=c(100,20,10))
Xc <- ncenter(X,nmode=2)   # center across columns
sum(rowSums(Xc))


##########   EXAMPLE 3   ##########
X <- array(rnorm(100000),dim=c(100,20,10,5))
Xc <- ncenter(X,nmode=4)   # center across 4-th mode
sum(rowSums(Xc))


##########   EXAMPLE 4   ##########
X <- replicate(5,array(rnorm(20000),dim=c(100,20,10)),simplify=FALSE)
Xc <- ncenter(X)
sum(colSums(Xc[[1]]))

}
