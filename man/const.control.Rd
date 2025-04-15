\name{const.control}
\alias{const.control}
\title{
Auxiliary for Controlling Multi-Way Constraints
}
\description{
Auxiliary function for controlling the \code{const} argument of the \code{\link{mcr}}, \code{\link{parafac}}, and \code{\link{parafac2}} functions. Applicable when using smoothness constraints.
}
\usage{
const.control(const, df = NULL, degree = NULL, intercept = NULL)
}
\arguments{
  \item{const}{
Character vector of length 3 or 4 giving the constraints for each mode. See \code{\link[CMLS]{const}} for the 24 available options.
}
  \item{df}{
Integer vector of length 3 or 4 giving the degrees of freedom to use for the spline basis in each mode. Can also input a single number giving the common degrees of freedom to use for each mode. Defaults to 7 degrees of freedom for each applicable mode.
}
  \item{degree}{
Integer vector of length 3 or 4 giving the polynomial degree to use for the spline basis in each mode. Can also input a single number giving the common polynomial degree to use for each mode. Defaults to degree 3 (cubic) polynomials for each applicable mode.
}
  \item{intercept}{
Logical vector of length 3 or 4 indicating whether the spline basis should contain an intercept. Can also input a single logical giving the common intercept indicator to use for each mode. Defaults to \code{TRUE} for each applicable mode.
}
}
\value{
Returns a list with elements: \code{const}, \code{df}, \code{degree}, and \code{intercept}.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\details{
The \code{\link{mcr}}, \code{\link{parafac}}, and \code{\link{parafac2}} functions pass the input \code{const} to this function to determine the fitting options when using smoothness constraints.

The \code{\link[CMLS]{const}} function (from \bold{CMLS} package) describes the available constraint options.
}
\examples{
##########   EXAMPLE   ##########

# create random data array with Parafac structure
set.seed(4)
mydim <- c(30, 10, 8, 10)
nf <- 4
aseq <- seq(-3, 3, length.out = mydim[1])
Amat <- cbind(dnorm(aseq), dchisq(aseq+3.1, df=3),
              dt(aseq-2, df=4), dgamma(aseq+3.1, shape=3, rate=1))
Bmat <- svd(matrix(runif(mydim[2]*nf), nrow = mydim[2], ncol = nf), nv = 0)$u
Cmat <- matrix(runif(mydim[3]*nf), nrow = mydim[3], ncol = nf)
Cstruc <- Cmat > 0.5
Cmat <- Cmat * Cstruc
Dmat <- matrix(runif(mydim[4]*nf), nrow = mydim[4], ncol = nf)
Xmat <- tcrossprod(Amat, krprod(Dmat, krprod(Cmat, Bmat)))
Xmat <- array(Xmat, dim = mydim)
Emat <- array(rnorm(prod(mydim)), dim = mydim)
Emat <- nscale(Emat, 0, ssnew = sumsq(Xmat))   # SNR = 1
X <- Xmat + Emat

# fit Parafac model (unimodal and smooth A, orthogonal B, 
#                    non-negative and structured C, non-negative D)
set.seed(123)
pfac <- parafac(X, nfac = nf, nstart = 1, Cstruc = Cstruc, 
                const = c("unismo", "orthog", "nonneg", "nonneg"))
pfac

# same as before, but add some options to the unimodality contraints...
# more knots (df=10), quadratic splines (degree=2), and enforce non-negativity
cvec <- c("unsmno", "orthog", "nonneg", "nonneg")
ctrl <- const.control(cvec, df = 10, degree = 2)
set.seed(123)
pfac <- parafac(X, nfac = nf, nstart = 1, Cstruc = Cstruc,
                const = cvec, control = ctrl)
pfac

}
