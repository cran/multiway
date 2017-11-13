\name{parafac}
\alias{parafac}
\title{
Parallel Factor Analysis-1
}
\description{
Given a 3-way array \code{X = array(x,dim=c(I,J,K))}, the 3-way Parafac model can be written as 
\tabular{c}{
\code{ X[i,j,k] = sum A[i,r]*B[j,r]*C[k,r] + E[i,j,k] }
} 
where \code{A = matrix(a,I,R)} are the Mode A (first mode) weights, \code{B = matrix(b,J,R)} are the Mode B (second mode) weights, \code{C = matrix(c,K,R)} are the Mode C (third mode) weights, and \code{E = array(e,dim=c(I,J,K))} is the 3-way residual array. The summation is for \code{r = seq(1,R)}.

Given a 4-way array \code{X = array(x,dim=c(I,J,K,L))}, the 4-way Parafac model can be written as 
\tabular{c}{
\code{ X[i,j,k,l] = sum A[i,r]*B[j,r]*C[k,r]*D[l,r] + E[i,j,k,l] }
} 
where \code{D = matrix(d,L,R)} are the Mode D (fourth mode) weights, \code{E = array(e,dim=c(I,J,K,L))} is the 4-way residual array, and the other terms can be interprered as previously described.

Weight matrices are estimated using an alternating least squares algorithm with optional constraints.
}
\usage{
parafac(X, nfac, nstart = 10, const = NULL, control = NULL,
        Bfixed = NULL, Cfixed = NULL, Dfixed = NULL,
        Bstart = NULL, Cstart = NULL, Dstart = NULL,
        Bstruc = NULL, Cstruc = NULL, Dstruc = NULL,
        maxit = 500, ctol = 1e-4, parallel = FALSE, cl = NULL,
        output = c("best", "all"), verbose = FALSE)
}
\arguments{
  \item{X}{
  Three-way data array with \code{dim=c(I,J,K)} or four-way data array with \code{dim=c(I,J,K,L)}. Missing data are allowed (see Note).
}
  \item{nfac}{
  Number of factors.
}
  \item{nstart}{
  Number of random starts.
}
  \item{const}{
  Constraints for each mode. Vector of length 3 or 4 with entries:
  0 = unconstrained (default),
  1 = orthogonal,
  2 = non-negative,
  3 = unimodal,
  4 = monotonic,
  5 = periodic,
  6 = smooth.
  Use \code{control} argument to adjust options for constraints 3-6.
}
\item{control}{
  List of parameters controlling options for constraints 3-6. This is passed to \code{\link{const.control}}, which describes the available options.
}
  \item{Bfixed}{
  Fixed Mode B weights. Only used to fit model with fixed Mode B weights.
}
  \item{Cfixed}{
  Fixed Mode C weights. Only used to fit model with fixed Mode C weights.
}
  \item{Dfixed}{
  Fixed Mode D weights. Only used to fit model with fixed Mode D weights. 
}
  \item{Bstart}{
  Starting Mode B weights for ALS algorithm. Default uses random weights.
}
  \item{Cstart}{
  Starting Mode C weights for ALS algorithm. Default uses random weights.
}
  \item{Dstart}{
  Starting Mode D weights for ALS algorithm. Default uses random weights.
}
  \item{Bstruc}{
  Structure constraints for Mode B weights. Default uses unstructured weights.
}
  \item{Cstruc}{
  Structure constraints for Mode C weights. Default uses unstructured weights.
}
  \item{Dstruc}{
  Structure constraints for Mode D weights. Default uses unstructured weights.
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
  Logical indicating if extra information on progress should be reported. Ignored if \code{parallel=TRUE}.
}
}
\value{
If \code{output="best"}, returns an object of class \code{"parafac"} with the following elements:
\item{A}{Mode A weight matrix.}
\item{B}{Mode B weight matrix.}
\item{C}{Mode C weight matrix.}
\item{D}{Mode D weight matrix.}
\item{SSE}{Sum of Squared Errors.}
\item{Rsq}{R-squared value.}
\item{GCV}{Generalized Cross-Validation.}
\item{edf}{Effective degrees of freedom.}
\item{iter}{Number of iterations.}
\item{cflag}{Convergence flag.}
\item{const}{See argument \code{const}.}
\item{control}{See argument \code{control}.}
\item{fixed}{Logical vector indicating whether 'fixed' weights were used for each mode.}
\item{struc}{Logical vector indicating whether 'struc' constraints were used for each mode.}

Otherwise returns a list of length \code{nstart} where each element is an object of class \code{"parafac"}.
}
\references{
Bro, R., & De Jong, S. (1997). A fast non-negativity-constrained least squares algorithm. \emph{Journal of Chemometrics, 11}, 393-401.

Harshman, R. A. (1970). Foundations of the PARAFAC procedure: Models and conditions for an "explanatory" multimodal factor analysis. \emph{UCLA Working Papers in Phonetics, 16}, 1-84.

Harshman, R. A., & Lundy, M. E. (1994). PARAFAC: Parallel factor analysis. \emph{Computational Statistics and Data Analysis, 18}, 39-72.

Helwig, N. E. (2017). Estimating latent trends in multivariate longitudinal data via Parafac2 with functional and structural constraints. \emph{Biometrical Journal, 59}(4), 783-803.

Ramsay, J. O. (1988). Monotone regression splines in action. \emph{Statistical Science, 3}, 425-441.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Default use is 10 random strarts (\code{nstart=10}) with 500 maximum iterations of the ALS algorithm for each start (\code{maxit=500}) using a convergence tolerance of 1e-4 (\code{ctol=1e-4}). The algorithm is determined to have converged once the change in R^2 is less than or equal to \code{ctol}.

Output \code{cflag} gives convergence information: \code{cflag=0} if ALS algorithm converged normally, \code{cflag=1} if maximum iteration limit was reached before convergence, and \code{cflag=2} if ALS algorithm terminated abnormally due to a problem with the constraints.

Constraints 3 (unimodality) and 4 (monotonicity) are implemented using I-splines, and constraints 5 (periodicity) and 6 (smoothness) are implemented using M-splines (see Ramsay, 1988).

Missing data should be specified as \code{NA} values in the input \code{X}. The missing data are randomly initialized and then iteratively imputed as a part of the ALS algorithm.
}
\section{Warnings }{
The ALS algorithm can perform poorly if the number of factors \code{nfac} is set too large.

Non-negativity constraints can be sensitive to local optima. 

Non-negativity constraints can result in slower performance.

Structure constraints for override constraints in \code{const} input.
}

\examples{
##########   3-way example   ##########

# create random data array with Parafac structure
set.seed(3)
mydim <- c(50,20,5)
nf <- 3
Amat <- matrix(rnorm(mydim[1]*nf),mydim[1],nf)
Bmat <- matrix(runif(mydim[2]*nf),mydim[2],nf)
Cmat <- matrix(runif(mydim[3]*nf),mydim[3],nf)
Xmat <- array(tcrossprod(Amat,krprod(Cmat,Bmat)),dim=mydim)
Emat <- array(rnorm(prod(mydim)),dim=mydim)
Emat <- nscale(Emat,0,sumsq(Xmat))   # SNR=1
X <- Xmat + Emat

# fit Parafac model (unconstrained)
pfac <- parafac(X,nfac=nf,nstart=1)
pfac

# fit Parafac model (non-negativity on Modes B and C)
pfacNN <- parafac(X,nfac=nf,nstart=1,const=c(0,2,2))
pfacNN

# check solution
Xhat <- fitted(pfac)
sum((Xmat-Xhat)^2)/prod(mydim)

# reorder and resign factors
pfac$B[1:4,]
pfac <- reorder(pfac, c(3,1,2))
pfac$B[1:4,]
pfac <- resign(pfac, mode="B")
pfac$B[1:4,]
Xhat <- fitted(pfac)
sum((Xmat-Xhat)^2)/prod(mydim)

# rescale factors
colSums(pfac$B^2)
colSums(pfac$C^2)
pfac <- rescale(pfac, mode="C", absorb="B")
colSums(pfac$B^2)
colSums(pfac$C^2)
Xhat <- fitted(pfac)
sum((Xmat-Xhat)^2)/prod(mydim)


##########   4-way example   ##########

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

# check solution
Xhat <- fitted(pfac)
sum((Xmat-Xhat)^2)/prod(mydim)
congru(Amat, pfac$A)
crossprod(pfac$B)
pfac$C


\dontrun{

##########   parallel computation   ##########

# create random data array with Parafac structure
set.seed(3)
mydim <- c(50,20,5)
nf <- 3
Amat <- matrix(rnorm(mydim[1]*nf),mydim[1],nf)
Bmat <- matrix(runif(mydim[2]*nf),mydim[2],nf)
Cmat <- matrix(runif(mydim[3]*nf),mydim[3],nf)
Xmat <- array(tcrossprod(Amat,krprod(Cmat,Bmat)),dim=mydim)
Emat <- array(rnorm(prod(mydim)),dim=mydim)
Emat <- nscale(Emat,0,sumsq(Xmat))   # SNR=1
X <- Xmat + Emat

# fit Parafac model (10 random starts -- sequential computation)
set.seed(1)
system.time({pfac <- parafac(X,nfac=nf)})
pfac

# fit Parafac model (10 random starts -- parallel computation)
set.seed(1)
cl <- makeCluster(detectCores())
ce <- clusterEvalQ(cl,library(multiway))
system.time({pfac <- parafac(X,nfac=nf,parallel=TRUE,cl=cl)})
pfac
stopCluster(cl)
}

}
