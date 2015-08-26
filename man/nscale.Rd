\name{nscale}
\alias{nscale}
\alias{nrescale}
\title{
Scale n-th Dimension of Array
}
\description{
 Slab-scale within each level of the specified mode. Can input 2-way, 3-way, and 4-way arrays, or input a list containing array elements.
 
 With \code{X} a matrix (I-by-J) there are two options:
 \tabular{rcl}{
 \code{ mode=1: } \tab \tab \code{ x[i,j]*sqrt(ssnew/sumsq(x[i,])) } \cr
 \code{ mode=2: } \tab \tab \code{ x[i,j]*sqrt(ssnew/sumsq(x[,j])) } \cr
 }
 
 With \code{X} a 3-way array (I-by-J-by-K) there are three options:
 \tabular{rcl}{
 \code{ mode=1: } \tab \tab \code{ x[i,j,k]*sqrt(ssnew/sumsq(x[i,,])) } \cr
 \code{ mode=2: } \tab \tab \code{ x[i,j,k]*sqrt(ssnew/sumsq(x[,j,]))) } \cr
 \code{ mode=3: } \tab \tab \code{ x[i,j,k]*sqrt(ssnew/sumsq(x[,,k])) } \cr
 }
 
 With \code{X} a 4-way array (I-by-J-by-K-by-L) there are four options:
 \tabular{rcl}{
 \code{ mode=1: } \tab \tab \code{ x[i,j,k,l]*sqrt(ssnew/sumsq(x[i,,,])) } \cr
 \code{ mode=2: } \tab \tab \code{ x[i,j,k,l]*sqrt(ssnew/sumsq(x[,j,,])) } \cr
 \code{ mode=3: } \tab \tab \code{ x[i,j,k,l]*sqrt(ssnew/sumsq(x[,,k,])) } \cr
 \code{ mode=4: } \tab \tab \code{ x[i,j,k,l]*sqrt(ssnew/sumsq(x[,,,l])) } \cr
 }
 
}
\usage{
nscale(X,mode=1,ssnew=1)
nrescale(X,mode=1,ssnew=1)
}
\arguments{
  \item{X}{
	Array (2-way, 3-way, or 4-way) or a list containing array elements.
}
  \item{mode}{
	Mode to scale within (set \code{mode=0} to scale across all modes).
}
  \item{ssnew}{
  Desired sum-of-squares for each level of scaled mode.
}
}
\value{
	Returns scaled version of \code{X}.
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
Xr <- nscale(X,mode=2)          # scale columns to ssnew=1
colSums(Xr^2)
Xr <- nscale(X,mode=2,ssnew=2)  # scale columns to ssnew=2
colSums(Xr^2)


##########   EXAMPLE 2   ##########
Xold <- X <- matrix(rnorm(400),20,20)
iter <- 0
chk <- 1
# iterative rescaling of modes 1 and 2
while(iter<500 & chk>=10^-9){
  Xr <- nscale(Xold,mode=1)
  Xr <- nscale(Xr,mode=2)
  chk <- sum((Xold-Xr)^2)
  Xold <- Xr
  iter <- iter + 1
}
iter
rowSums(Xr^2)
colSums(Xr^2)


##########   EXAMPLE 3   ##########
X <- array(rnorm(20000),dim=c(100,20,10))
Xc <- nscale(X,mode=2)   # scale within columns
rowSums(aperm(Xc,perm=c(2,1,3))^2)


##########   EXAMPLE 4   ##########
X <- array(rnorm(100000),dim=c(100,20,10,5))
Xc <- nscale(X,mode=4)   # scale across 4-th mode
rowSums(aperm(Xc,perm=c(4,1,2,3))^2)


##########   EXAMPLE 5   ##########
X <- replicate(5,array(rnorm(20000),dim=c(100,20,10)),simplify=FALSE)
Xc <- nscale(X,ssnew=200)   # mean square of 1
rowSums(Xc[[1]]^2)
rowMeans(Xc[[1]]^2)

}
