\name{const.control}
\alias{const.control}
\title{
Auxiliary for Controlling Multi-Way Constraints
}
\description{
Auxiliary function for controlling the \code{const} argument of the \code{\link{parafac}} and \code{\link{parafac2}} functions. Applicable when using constraints 3 (unimodal), 4 (monotonic), 5 (periodic), and/or 6 (smooth).
}
\usage{
const.control(const, df = NULL, degree = NULL, nonneg = NULL)
}
\arguments{
  \item{const}{
Constraints for each mode. Vector of length 3 or 4 with entries:
  0 = unconstrained (default),
  1 = orthogonal,
  2 = non-negative,
  3 = unimodal,
  4 = monotonic,
  5 = periodic,
  6 = smooth.
}
  \item{df}{
Integer vector of length 3 or 4 giving the degrees of freedom to use for the spline basis in each mode. Can also input a single number giving the common degrees of freedom to use for each mode. Defaults to 7 degrees of freedom for each applicable mode.
}
  \item{degree}{
Integer vector of length 3 or 4 giving the polynomial degree to use for the spline basis in each mode. Can also input a single number giving the common polynomial degree to use for each mode. Defaults to degree 3 (cubic) polynomials for each applicable mode.
}
  \item{nonneg}{
Logical vector of length 3 or 4 indicating whether the weights in each mode should be constrained to be non-negative. Can also input a single logical giving the common non-negativity constraints to use for each mode. Defaults to \code{FALSE} for each applicable mode.
}
}
\value{
Returns a list with elements: \code{const}, \code{df}, \code{degree}, and \code{nonneg}.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\details{
The \code{\link{parafac}} and \code{\link{parafac2}} functions pass the input \code{control} to this function to determine the fitting options when using constraints 3-6.
}
\examples{
##########   EXAMPLE   ##########

# create random data array with Parafac structure
set.seed(4)
mydim <- c(30,10,8,10)
nf <- 4
aseq <- seq(-3,3,length=mydim[1])
Amat <- cbind(dnorm(aseq), dchisq(aseq+3.1, df=3),
              dt(aseq-2, df=4), dgamma(aseq+3.1, shape=3, rate=1))
Bmat <- svd(matrix(runif(mydim[2]*nf),mydim[2],nf))$u
Cmat <- matrix(runif(mydim[3]*nf),mydim[3],nf)
Dmat <- matrix(runif(mydim[4]*nf),mydim[4],nf)
Xmat <- array(tcrossprod(Amat,krprod(Dmat,krprod(Cmat,Bmat))),dim=mydim)
Emat <- array(rnorm(prod(mydim)),dim=mydim)
Emat <- nscale(Emat,0,sumsq(Xmat))   # SNR=1
X <- Xmat + Emat

# fit Parafac model (unimodal A, orthogonal B, non-negative C, non-negative D)
pfac <- parafac(X,nfac=nf,nstart=1,const=c(3,1,2,2))
pfac

# same as before, but add some options to the unimodality contraints...
# fewer knots (df=5), quadratic splines (degree=2), and enforce non-negativity
cc <- const.control(c(3,1,2,2), df=5, degree=2, nonneg=TRUE)
pfac <- parafac(X,nfac=nf,nstart=1,const=c(3,1,2,2),control=cc)
pfac

}
