\name{nscale}
\alias{nscale}
\title{
Scale n-th Dimension of Array
}
\description{
Slab-scale within each level of the specified mode. Can input 2-way, 3-way, and 4-way arrays, or input a list containing array elements (see Note).
}
\usage{
nscale(X, mode = 1, ssnew = NULL, newscale = 1)
}
\arguments{
  \item{X}{
	Array (2-way, 3-way, or 4-way) or a list containing array elements.
}
  \item{mode}{
	Mode to scale within (set \code{mode = 0} to scale across all modes).
}
  \item{ssnew}{
  Desired sum-of-squares for each level of scaled mode.
}
  \item{newscale}{
  Desired root-mean-square for each level of scaled mode. Ignored if \code{ssnew} is supplied.
}
}
\details{

Default (as of ver 1.0-5) uses \code{newscale} argument...

With \code{X} a matrix (I-by-J) there are two options:
 \tabular{rcl}{
 \code{ mode=1: } \tab \tab \code{ x[i,j] * newscale / sqrt(meansq(x[i,])) } \cr
 \code{ mode=2: } \tab \tab \code{ x[i,j] * newscale / sqrt(meansq(x[,j])) } \cr
 }
 
 With \code{X} a 3-way array (I-by-J-by-K) there are three options:
 \tabular{rcl}{
 \code{ mode=1: } \tab \tab \code{ x[i,j,k] * newscale / sqrt(meansq(x[i,,])) } \cr
 \code{ mode=2: } \tab \tab \code{ x[i,j,k] * newscale / sqrt(meansq(x[,j,]))) } \cr
 \code{ mode=3: } \tab \tab \code{ x[i,j,k] * newscale / sqrt(meansq(x[,,k])) } \cr
 }
 
 With \code{X} a 4-way array (I-by-J-by-K-by-L) there are four options:
 \tabular{rcl}{
 \code{ mode=1: } \tab \tab \code{ x[i,j,k,l] * newscale / sqrt(meansq(x[i,,,])) } \cr
 \code{ mode=2: } \tab \tab \code{ x[i,j,k,l] * newscale / sqrt(meansq(x[,j,,])) } \cr
 \code{ mode=3: } \tab \tab \code{ x[i,j,k,l] * newscale / sqrt(meansq(x[,,k,])) } \cr
 \code{ mode=4: } \tab \tab \code{ x[i,j,k,l] * newscale / sqrt(meansq(x[,,,l])) } \cr
 }

If argument \code{ssnew} is provided...

 With \code{X} a matrix (I-by-J) there are two options:
 \tabular{rcl}{
 \code{ mode=1: } \tab \tab \code{ x[i,j] * sqrt(ssnew / sumsq(x[i,])) } \cr
 \code{ mode=2: } \tab \tab \code{ x[i,j] * sqrt(ssnew / sumsq(x[,j])) } \cr
 }
 
 With \code{X} a 3-way array (I-by-J-by-K) there are three options:
 \tabular{rcl}{
 \code{ mode=1: } \tab \tab \code{ x[i,j,k] * sqrt(ssnew / sumsq(x[i,,])) } \cr
 \code{ mode=2: } \tab \tab \code{ x[i,j,k] * sqrt(ssnew / sumsq(x[,j,]))) } \cr
 \code{ mode=3: } \tab \tab \code{ x[i,j,k] * sqrt(ssnew / sumsq(x[,,k])) } \cr
 }
 
 With \code{X} a 4-way array (I-by-J-by-K-by-L) there are four options:
 \tabular{rcl}{
 \code{ mode=1: } \tab \tab \code{ x[i,j,k,l] * sqrt(ssnew / sumsq(x[i,,,])) } \cr
 \code{ mode=2: } \tab \tab \code{ x[i,j,k,l] * sqrt(ssnew / sumsq(x[,j,,])) } \cr
 \code{ mode=3: } \tab \tab \code{ x[i,j,k,l] * sqrt(ssnew / sumsq(x[,,k,])) } \cr
 \code{ mode=4: } \tab \tab \code{ x[i,j,k,l] * sqrt(ssnew / sumsq(x[,,,l])) } \cr
 }
 
}
\value{
	Returns scaled version of \code{X}.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
	When entering a list with array elements, each element must be a 2-way or 3-way array. The list elements are treated as the 3rd mode (for list of 2-way arrays) or the 4th mode (for list of 3-way arrays) in the formulas provided in the Description.
}

\examples{
##########   EXAMPLE 1   ##########
X <- matrix(rnorm(2000), nrow = 100, ncol = 20)
Xr <- nscale(X, mode = 2)                # scale columns to newscale=1
sqrt(colMeans(Xr^2))
Xr <- nscale(X, mode = 2, newscale = 2)  # scale columns to newscale=2
sqrt(colMeans(Xr^2))


##########   EXAMPLE 2   ##########
Xold <- X <- matrix(rnorm(400), nrow = 20, ncol = 20)
iter <- 0
chk <- 1
# iterative scaling of modes 1 and 2
while(iter<500 & chk>=10^-9){
  Xr <- nscale(Xold, mode = 1)
  Xr <- nscale(Xr, mode = 2)
  chk <- sum((Xold-Xr)^2)
  Xold <- Xr
  iter <- iter + 1
}
iter
sqrt(rowMeans(Xr^2))
sqrt(colMeans(Xr^2))


##########   EXAMPLE 3   ##########
X <- array(rnorm(20000), dim = c(100,20,10))
Xc <- nscale(X, mode = 2)   # scale within columns
sqrt(rowMeans(aperm(Xc, perm = c(2,1,3))^2))


##########   EXAMPLE 4   ##########
X <- array(rnorm(100000), dim = c(100,20,10,5))
Xc <- nscale(X, mode = 4)   # scale across 4-th mode
sqrt(rowMeans(aperm(Xc, perm = c(4,1,2,3))^2))


##########   EXAMPLE 5   ##########
X <- replicate(5, array(rnorm(20000), dim = c(100,20,10)), simplify = FALSE)
# mean square of 1 (new way)
Xc <- nscale(X)
rowSums(sapply(Xc, function(x) rowSums(x^2))) / (20*10*5)
# mean square of 1 (old way)
Xc <- nscale(X, ssnew = (20*10*5))
rowSums(sapply(Xc, function(x) rowSums(x^2))) / (20*10*5)



}

\keyword{ array }% use one of  RShowDoc("KEYWORDS")
\keyword{ algebra }% use one of  RShowDoc("KEYWORDS")
