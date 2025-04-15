\name{cpd}
\alias{cpd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
N-way Canonical Polyadic Decomposition
}
\description{
Fits Frank L. Hitchcock's Canonical Polyadic Decomposition (CPD) to N-way data arrays. Parameters are estimated via alternating least squares.
}
\usage{
cpd(X, nfac, nstart = 10, maxit = 500, 
    ctol = 1e-4, parallel = FALSE, cl = NULL, 
    output = "best", verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  N-way data array. Missing data are allowed (see Note).
}
\item{nfac}{
  Number of factors.
}
  \item{nstart}{
  Number of random starts.
}
  \item{maxit}{
  Maximum number of iterations.
}
  \item{ctol}{
  Convergence tolerance (R^2 change).
}
  \item{parallel}{
  Logical indicating if \code{\link{parLapply}} should be used. See Examples.
}
  \item{cl}{
  Cluster created by \code{\link{makeCluster}}. Only used when \code{parallel=TRUE}.
}
  \item{output}{
  Output the best solution (default) or output all \code{nstart} solutions.
}
  \item{verbose}{
  If \code{TRUE}, fitting progress is printed via \code{\link{txtProgressBar}}. Ignored if \code{parallel=TRUE}.
}
}
\details{
This is an N-way extension of the \code{\link{parafac}} function without constraints. The form of the CPD for 3-way and 4-way data is given in the documentation for the \code{\link{parafac}} function. For N > 4, the CPD has the form
\tabular{c}{
\code{ X[i.1, ..., i.N] = sum A.1[i.1,r] * ... * A.N[i.N,r] + E[i.1, ..., i.N] }
} 
where \code{A.n} are the n-th mode's weights for n = 1, ..., N, and \code{E} is the N-way residual array. The summation is for \code{r = seq(1,R)}.
}
\value{
\item{A}{List of length N containing the weights for each mode.}
\item{SSE}{Sum of Squared Errors.}
\item{Rsq}{R-squared value.}
\item{GCV}{Generalized Cross-Validation.}
\item{edf}{Effective degrees of freedom.}
\item{iter}{Number of iterations.}
\item{cflag}{Convergence flag. See Note.}
}
\references{
Harshman, R. A. (1970). Foundations of the PARAFAC procedure: Models and conditions for an "explanatory" multimodal factor analysis. \emph{UCLA Working Papers in Phonetics, 16}, 1-84.

Harshman, R. A., & Lundy, M. E. (1994). PARAFAC: Parallel factor analysis. \emph{Computational Statistics and Data Analysis, 18}, 39-72. \doi{10.1016/0167-9473(94)90132-5}

Hitchcock, F. L. (1927). The expression of a tensor or a polyadic as a sum of products. \emph{Journal of Mathematics and Physics, 6}, 164-189. \doi{10.1002/sapm192761164}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Missing data should be specified as \code{NA} values in the input \code{X}. The missing data are randomly initialized and then iteratively imputed as a part of the algorithm.

Output \code{cflag} gives convergence information: \code{cflag = 0} if algorithm converged normally and \code{cflag = 1} if maximum iteration limit was reached before convergence.
}

\section{Warnings }{
The algorithm can perform poorly if the number of factors \code{nfac} is set too large.
}

\seealso{
The \code{\link{parafac}} function provides a more flexible implemention for 3-way and 4-way arrays.

The \code{\link{fitted.cpd}} function creates the model-implied fitted values from a fit \code{"cpd"} object.

The \code{\link{resign.cpd}} function can be used to resign factors from a fit \code{"cpd"} object.

The \code{\link{rescale.cpd}} function can be used to rescale factors from a fit \code{"cpd"} object.

The \code{\link{reorder.cpd}} function can be used to reorder factors from a fit \code{"cpd"} object.
}
\examples{
##########   3-way example   ##########

# create random data array with CPD/Parafac structure
set.seed(3)
mydim <- c(50, 20, 5)
nf <- 3
Amat <- matrix(rnorm(mydim[1]*nf), nrow = mydim[1], ncol = nf)
Bmat <- matrix(runif(mydim[2]*nf), nrow = mydim[2], ncol = nf)
Cmat <- matrix(runif(mydim[3]*nf), nrow = mydim[3], ncol = nf)
Xmat <- tcrossprod(Amat, krprod(Cmat, Bmat))
Xmat <- array(Xmat, dim = mydim)
Emat <- array(rnorm(prod(mydim)), dim = mydim)
Emat <- nscale(Emat, 0, ssnew = sumsq(Xmat))   # SNR = 1
X <- Xmat + Emat

# fit CPD model
set.seed(0)
cano <- cpd(X, nfac = nf, nstart = 1)
cano

# fit Parafac model
set.seed(0)
pfac <- parafac(X, nfac = nf, nstart = 1)
pfac


##########   4-way example   ##########

# create random data array with CPD/Parafac structure
set.seed(4)
mydim <- c(30,10,8,10)
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

# fit CPD model
set.seed(0)
cano <- cpd(X, nfac = nf, nstart = 1)
cano

# fit Parafac model
set.seed(0)
pfac <- parafac(X, nfac = nf, nstart = 1)
pfac


##########   5-way example   ##########

# create random data array with CPD/Parafac structure
set.seed(5)
mydim <- c(5, 6, 7, 8, 9)
nmode <- length(mydim)
nf <- 3
Amat <- vector("list", nmode)
for(n in 1:nmode) {
  Amat[[n]] <- matrix(rnorm(mydim[n] * nf), mydim[n], nf)
}
Zmat <- krprod(Amat[[3]], Amat[[2]])
for(n in 4:5) Zmat <- krprod(Amat[[n]], Zmat)
Xmat <- tcrossprod(Amat[[1]], Zmat)
Xmat <- array(Xmat, dim = mydim)
Emat <- array(rnorm(prod(mydim)), dim = mydim)
Emat <- nscale(Emat, 0, ssnew = sumsq(Xmat))   # SNR = 1
X <- Xmat + Emat

# fit CPD model
set.seed(0)
cano <- cpd(X, nfac = nf, nstart = 1)
cano


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
\keyword{ multivariate }% __ONLY ONE__ keyword per line
