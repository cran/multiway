\name{corcondia}
\alias{corcondia}
\title{
Core Consistency Diagnostic
}
\description{
Calculates Bro and Kiers's core consistency diagnostic (CORCONDIA) for a fit \code{\link{parafac}}  or \code{\link{parafac2}} model. For Parafac2, the diagnostic is calculated after transforming the data by \code{pfac$A[[k]] \%*\% smpower(pfac$Phi, -0.5)}.
}
\usage{
corcondia(X, object, divisor=c("nfac","core"))
}
\arguments{
  \item{X}{
Three-way data array with \code{dim=c(I,J,K)} or four-way data array with \code{dim=c(I,J,K,L)}. Can also input a list of two-way or three-way arrays (for Parafac2).
}
  \item{object}{
Object of class "parafac" (output from \code{\link{parafac}}) or class "parafac2" (output from \code{\link{parafac2}}).
}
  \item{divisor}{
Divide by number of factors (default) or core sum of squares.
}
}
\value{
Returns CORCONDIA value.
}
\references{
Bro, R., & Kiers, H.A.L. (2003). A new efficient method for determining the number of components in PARAFAC models. \emph{Journal of Chemometrics, 17}, 274-286.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\details{
The core consistency diagnostic is defined as
\tabular{c}{
\code{ 100 * ( 1 - sum( (G-S)^2 ) / divisor ) }
}
where \code{G} is the least squares estimate of the Tucker core array, \code{S} is a super-diagonal core array, and \code{divisor} is the sum of squares of either \code{S} ("nfac") or \code{G} ("core"). A value of 100 indiciates a perfect multilinear structure, and smaller values indicate greater violations of multilinear structure.
}
\examples{
##########   EXAMPLE   ##########

# create random data array with Parafac structure
set.seed(3)
mydim <- c(50,20,5)
nf <- 2
Amat <- matrix(rnorm(mydim[1]*nf),mydim[1],nf)
Bmat <- matrix(runif(mydim[2]*nf),mydim[2],nf)
Cmat <- matrix(runif(mydim[3]*nf),mydim[3],nf)
Xmat <- array(tcrossprod(Amat,krprod(Cmat,Bmat)),dim=mydim)
Emat <- array(rnorm(prod(mydim)),dim=mydim)
Emat <- nscale(Emat,0,sumsq(Xmat))   # SNR=1
X <- Xmat + Emat

# fit Parafac model (1-4 factors)
pfac1 <- parafac(X,nfac=1,nstart=1)
pfac2 <- parafac(X,nfac=2,nstart=1)
pfac3 <- parafac(X,nfac=3,nstart=1)
pfac4 <- parafac(X,nfac=4,nstart=1)

# check corcondia
corcondia(X, pfac1)
corcondia(X, pfac2)
corcondia(X, pfac3)
corcondia(X, pfac4)

}
