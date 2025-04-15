\name{parafac}
\alias{parafac}
\title{
Parallel Factor Analysis-1
}
\description{
Fits Richard A. Harshman's Parallel Factors (Parafac) model to 3-way or 4-way data arrays. Parameters are estimated via alternating least squares with optional constraints.
}
\usage{
parafac(X, nfac, nstart = 10, const = NULL, control = NULL,
        Afixed = NULL, Bfixed = NULL, Cfixed = NULL, Dfixed = NULL,
        Astart = NULL, Bstart = NULL, Cstart = NULL, Dstart = NULL,
        Astruc = NULL, Bstruc = NULL, Cstruc = NULL, Dstruc = NULL,
        Amodes = NULL, Bmodes = NULL, Cmodes = NULL, Dmodes = NULL,
        maxit = 500, ctol = 1e-4, parallel = FALSE, cl = NULL,
        output = c("best", "all"), verbose = TRUE, backfit = FALSE)
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
  Character vector of length 3 or 4 giving the constraints for each mode (defaults to unconstrained). See \code{\link[CMLS]{const}} for the 24 available options.
}
\item{control}{
  List of parameters controlling options for smoothness constraints. This is passed to \code{\link{const.control}}, which describes the available options.
}
  \item{Afixed}{
  Used to fit model with fixed Mode A weights.
}
  \item{Bfixed}{
  Used to fit model with fixed Mode B weights.
}
  \item{Cfixed}{
  Used to fit model with fixed Mode C weights.
}
  \item{Dfixed}{
  Used to fit model with fixed Mode D weights. 
}
  \item{Astart}{
  Starting Mode A weights. Default uses random weights.
}
  \item{Bstart}{
  Starting Mode B weights. Default uses random weights.
}
  \item{Cstart}{
  Starting Mode C weights. Default uses random weights.
}
  \item{Dstart}{
  Starting Mode D weights. Default uses random weights.
}
  \item{Astruc}{
  Structure constraints for Mode A weights. See Note.
}
  \item{Bstruc}{
  Structure constraints for Mode B weights. See Note.
}
  \item{Cstruc}{
  Structure constraints for Mode C weights. See Note.
}
  \item{Dstruc}{
  Structure constraints for Mode D weights. See Note.
}
  \item{Amodes}{
  Mode ranges for Mode A weights (for unimodality constraints). See Note.
}
  \item{Bmodes}{
  Mode ranges for Mode B weights (for unimodality constraints). See Note.
}
  \item{Cmodes}{
  Mode ranges for Mode C weights (for unimodality constraints). See Note.
}
  \item{Dmodes}{
  Mode ranges for Mode D weights (for unimodality constraints). See Note.
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
  \item{backfit}{
  Should backfitting algorithm be used for \code{\link[CMLS]{cmls}}?
}
}

\details{
Given a 3-way array \code{X = array(x, dim = c(I,J,K))}, the 3-way Parafac model can be written as 
\tabular{c}{
\code{ X[i,j,k] = sum A[i,r]*B[j,r]*C[k,r] + E[i,j,k] }
} 
where \code{A = matrix(a,I,R)} are the Mode A (first mode) weights, \code{B = matrix(b,J,R)} are the Mode B (second mode) weights, \code{C = matrix(c,K,R)} are the Mode C (third mode) weights, and \code{E = array(e,dim=c(I,J,K))} is the 3-way residual array. The summation is for \code{r = seq(1,R)}.

Given a 4-way array \code{X = array(x, dim = c(I,J,K,L))}, the 4-way Parafac model can be written as 
\tabular{c}{
\code{ X[i,j,k,l] = sum A[i,r]*B[j,r]*C[k,r]*D[l,r] + E[i,j,k,l] }
} 
where \code{D = matrix(d,L,R)} are the Mode D (fourth mode) weights, \code{E = array(e,dim=c(I,J,K,L))} is the 4-way residual array, and the other terms can be interprered as previously described.

Weight matrices are estimated using an alternating least squares algorithm with optional constraints.
}

\value{
If \code{output = "best"}, returns an object of class \code{"parafac"} with the following elements:
\item{A}{Mode A weight matrix.}
\item{B}{Mode B weight matrix.}
\item{C}{Mode C weight matrix.}
\item{D}{Mode D weight matrix.}
\item{SSE}{Sum of Squared Errors.}
\item{Rsq}{R-squared value.}
\item{GCV}{Generalized Cross-Validation.}
\item{edf}{Effective degrees of freedom.}
\item{iter}{Number of iterations.}
\item{cflag}{Convergence flag. See Note.}
\item{const}{See argument \code{const}.}
\item{control}{See argument \code{control}.}
\item{fixed}{Logical vector indicating whether 'fixed' weights were used for each mode.}
\item{struc}{Logical vector indicating whether 'struc' constraints were used for each mode.}

Otherwise returns a list of length \code{nstart} where each element is an object of class \code{"parafac"}.
}
\references{
Harshman, R. A. (1970). Foundations of the PARAFAC procedure: Models and conditions for an "explanatory" multimodal factor analysis. \emph{UCLA Working Papers in Phonetics, 16}, 1-84.

Harshman, R. A., & Lundy, M. E. (1994). PARAFAC: Parallel factor analysis. \emph{Computational Statistics and Data Analysis, 18}, 39-72. \doi{10.1016/0167-9473(94)90132-5}

Helwig, N. E. (2017). Estimating latent trends in multivariate longitudinal data via Parafac2 with functional and structural constraints. \emph{Biometrical Journal, 59}(4), 783-803. \doi{10.1002/bimj.201600045}

Hitchcock, F. L. (1927). The expression of a tensor or a polyadic as a sum of products. \emph{Journal of Mathematics and Physics, 6}, 164-189. \doi{10.1002/sapm192761164}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Missing data should be specified as \code{NA} values in the input \code{X}. The missing data are randomly initialized and then iteratively imputed as a part of the algorithm.

Structure constraints should be specified with a matrix of logicals (TRUE/FALSE), such that FALSE elements indicate a weight should be constrained to be zero. Default uses unstructured weights, i.e., a matrix of all TRUE values. 

When using unimodal constraints, the \code{*modes} inputs can be used to specify the mode search range for each factor. These inputs should be matrices with dimension \code{c(2,nfac)} where the first row gives the minimum mode value and the second row gives the maximum mode value (with respect to the indicies of the corresponding weight matrix).

Output \code{cflag} gives convergence information: \code{cflag = 0} if algorithm converged normally, \code{cflag = 1} if maximum iteration limit was reached before convergence, and \code{cflag = 2} if algorithm terminated abnormally due to a problem with the constraints.
}
\section{Warnings }{
The algorithm can perform poorly if the number of factors \code{nfac} is set too large.
}

\seealso{
The \code{\link{cpd}} function implements an N-way extension without constraints.

The \code{\link{fitted.parafac}} function creates the model-implied fitted values from a fit \code{"parafac"} object.

The \code{\link{resign.parafac}} function can be used to resign factors from a fit \code{"parafac"} object.

The \code{\link{rescale.parafac}} function can be used to rescale factors from a fit \code{"parafac"} object.

The \code{\link{reorder.parafac}} function can be used to reorder factors from a fit \code{"parafac"} object.

The \code{\link[CMLS]{cmls}} function (from \bold{CMLS} package) is called as a part of the alternating least squares algorithm.
}

\examples{
##########   3-way example   ##########

# create random data array with Parafac structure
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

# fit Parafac model (unconstrained)
pfac <- parafac(X, nfac = nf, nstart = 1)
pfac

# fit Parafac model (non-negativity on Modes B and C)
pfacNN <- parafac(X, nfac = nf, nstart = 1, 
                  const = c("uncons", "nonneg", "nonneg"))
pfacNN

# check solution
Xhat <- fitted(pfac)
sum((Xmat - Xhat)^2) / prod(mydim)

# reorder and resign factors
pfac$B[1:4,]
pfac <- reorder(pfac, c(3,1,2))
pfac$B[1:4,]
pfac <- resign(pfac, mode="B")
pfac$B[1:4,]
Xhat <- fitted(pfac)
sum((Xmat - Xhat)^2) / prod(mydim)

# rescale factors
colSums(pfac$B^2)
colSums(pfac$C^2)
pfac <- rescale(pfac, mode = "C", absorb = "B")
colSums(pfac$B^2)
colSums(pfac$C^2)
Xhat <- fitted(pfac)
sum((Xmat - Xhat)^2) / prod(mydim)


##########   4-way example   ##########

# create random data array with Parafac structure
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

# fit Parafac model (unimodal and smooth A, orthogonal B, 
#                    non-negative and structured C, non-negative D)
pfac <- parafac(X, nfac = nf, nstart = 1, Cstruc = Cstruc, 
                const = c("unismo", "orthog", "nonneg", "nonneg"))
pfac

# check solution
Xhat <- fitted(pfac)
sum((Xmat - Xhat)^2) / prod(mydim)
congru(Amat, pfac$A)
crossprod(pfac$B)
pfac$C
Cstruc

\dontrun{

##########   parallel computation   ##########

# create random data array with Parafac structure
set.seed(3)
mydim <- c(50,20,5)
nf <- 3
Amat <- matrix(rnorm(mydim[1]*nf), nrow = mydim[1], ncol = nf)
Bmat <- matrix(runif(mydim[2]*nf), nrow = mydim[2], ncol = nf)
Cmat <- matrix(runif(mydim[3]*nf), nrow = mydim[3], ncol = nf)
Xmat <- tcrossprod(Amat, krprod(Cmat, Bmat))
Xmat <- array(Xmat, dim = mydim)
Emat <- array(rnorm(prod(mydim)), dim = mydim)
Emat <- nscale(Emat, 0, ssnew = sumsq(Xmat))   # SNR = 1
X <- Xmat + Emat

# fit Parafac model (10 random starts -- sequential computation)
set.seed(1)
system.time({pfac <- parafac(X, nfac = nf)})
pfac

# fit Parafac model (10 random starts -- parallel computation)
cl <- makeCluster(detectCores())
ce <- clusterEvalQ(cl, library(multiway))
clusterSetRNGStream(cl, 1)
system.time({pfac <- parafac(X, nfac = nf, parallel = TRUE, cl = cl)})
pfac
stopCluster(cl)
}

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")
