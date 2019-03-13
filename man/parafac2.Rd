\name{parafac2}
\alias{parafac2}
\title{
Parallel Factor Analysis-2
}
\description{
Fits Richard A. Harshman's Parallel Factors-2 (Parafac2) model to 3-way or 4-way ragged data arrays. Parameters are estimated via alternating least squares with optional constraints.
}
\usage{
parafac2(X, nfac, nstart = 10, const = NULL, control = NULL,
         Gfixed = NULL, Bfixed = NULL, Cfixed = NULL, Dfixed = NULL,
         Gstart = NULL, Bstart = NULL, Cstart = NULL, Dstart = NULL,
         Gstruc = NULL, Bstruc = NULL, Cstruc = NULL, Dstruc = NULL,
         Gmodes = NULL, Bmodes = NULL, Cmodes = NULL, Dmodes = NULL,
         maxit = 500, ctol = 1e-4, parallel = FALSE, cl = NULL,
         output = c("best", "all"), verbose = TRUE, backfit = FALSE)
}
\arguments{
  \item{X}{
  For 3-way Parafac2: list of length \code{K} where \code{k}-th element is \code{I[k]}-by-\code{J} matrix or three-way data array with \code{dim=c(I,J,K)}. For 4-way Parafac2: list of length \code{L} where \code{l}-th element is \code{I[l]}-by-\code{J}-by-\code{K} array or four-way data array with \code{dim=c(I,J,K,L)}. Missing data are allowed (see Note).
}
  \item{nfac}{
  Number of factors.
}
  \item{nstart}{
  Number of random starts.
}
  \item{const}{
  Character vector of length 3 or 4 giving the constraints for each mode (defaults to unconstrained). See \code{\link{const}} for the 24 available options. 
}
\item{control}{
  List of parameters controlling options for smoothness constraints. This is passed to \code{\link{const.control}}, which describes the available options.
}
  \item{Gfixed}{
  Used to fit model with fixed Phi matrix: \code{crossprod(Gfixed) = Phi}.
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
  \item{Gstart}{
  Starting Mode A crossproduct matrix: \code{crossprod(Gstart) = Phi}. Default uses random weights.
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
  \item{Gstruc}{
  Structure constraints for Mode A crossproduct matrix: \code{crossprod(Gstruc) = Phistruc}. See Note.
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
  \item{Gmodes}{
  Mode ranges for Mode A weights (for unimodality constraints). Ignored.
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
  Should backfitting algorithm be used for \code{\link{cmls}}?
}
}

\details{
Given a list of matrices \code{X[[k]] = matrix(xk,I[k],J)} for \code{k = seq(1,K)}, the 3-way Parafac2 model (with Mode A nested in Mode C) can be written as
\tabular{r}{
\code{X[[k]] = tcrossprod(A[[k]] \%*\% diag(C[k,]), B) + E[[k]]} \cr 
subject to \code{crossprod(A[[k]]) = Phi} \cr 
}
where \code{A[[k]] = matrix(ak,I[k],R)} are the Mode A (first mode) weights for the \code{k}-th level of Mode C (third mode), \code{Phi} is the common crossproduct matrix shared by all \code{K} levels of Mode C, \code{B = matrix(b,J,R)} are the Mode B (second mode) weights, \code{C = matrix(c,K,R)} are the Mode C (third mode) weights, and \code{E[[k]] = matrix(ek,I[k],J)} is the residual matrix corresponding to \code{k}-th level of Mode C.

Given a list of arrays \code{X[[l]] = array(xl, dim = c(I[l],J,K))} for \code{l = seq(1,L)}, the 4-way Parafac2 model (with Mode A nested in Mode D) can be written as
\tabular{r}{
\code{X[[l]][,,k] = tcrossprod(A[[l]] \%*\% diag(D[l,]*C[k,]), B) + E[[l]][,,k]} \cr 
subject to \code{crossprod(A[[l]]) = Phi} \cr 
}
where \code{A[[l]] = matrix(al,I[l],R)} are the Mode A (first mode) weights for the \code{l}-th level of Mode D (fourth mode), \code{Phi} is the common crossproduct matrix shared by all \code{L} levels of Mode D, \code{D = matrix(d,L,R)} are the Mode D (fourth mode) weights, and \code{E[[l]] = array(el, dim = c(I[l],J,K))} is the residual array corresponding to \code{l}-th level of Mode D.

Weight matrices are estimated using an alternating least squares algorithm with optional constraints.
}

\value{
If \code{output = "best"}, returns an object of class \code{"parafac2"} with the following elements:
\item{A}{List of Mode A weight matrices.}
\item{B}{Mode B weight matrix.}
\item{C}{Mode C weight matrix.}
\item{D}{Mode D weight matrix.}
\item{Phi}{Mode A crossproduct matrix.}
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

Otherwise returns a list of length \code{nstart} where each element is an object of class \code{"parafac2"}.
}
\references{
Harshman, R. A. (1972). PARAFAC2: Mathematical and technical notes. \emph{UCLA Working Papers in Phonetics, 22}, 30-44.

Helwig, N. E. (2013). The special sign indeterminacy of the direct-fitting Parafac2 model: Some implications, cautions, and recommendations, for Simultaneous Component Analysis. \emph{Psychometrika, 78}, 725-739.

Helwig, N. E. (2017). Estimating latent trends in multivariate longitudinal data via Parafac2 with functional and structural constraints. \emph{Biometrical Journal, 59}(4), 783-803.

Helwig, N. E. (in prep). Constrained parallel factor analysis via the R package multiway.

Kiers, H. A. L., ten Berge, J. M. F., & Bro, R. (1999). PARAFAC2-part I: A direct-fitting algorithm for the PARAFAC2 model. \emph{Journal of Chemometrics, 13}, 275-294.
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
The \code{\link{fitted.parafac2}} function creates the model-implied fitted values from a fit \code{"parafac2"} object.

The \code{\link{resign.parafac2}} function can be used to resign factors from a fit \code{"parafac2"} object.

The \code{\link{rescale.parafac2}} function can be used to rescale factors from a fit \code{"parafac2"} object.

The \code{\link{reorder.parafac2}} function can be used to reorder factors from a fit \code{"parafac2"} object.

The \code{\link{cmls}} function (from \bold{CMLS} package) is called as a part of the alternating least squares algorithm.
}

\examples{
##########   3-way example   ##########

# create random data list with Parafac2 structure
set.seed(3)
mydim <- c(NA, 10, 20)
nf <- 2
nk <- rep(c(50, 100, 200), length.out = mydim[3])
Gmat <- matrix(rnorm(nf^2), nrow = nf, ncol = nf)
Bmat <- matrix(runif(mydim[2]*nf), nrow = mydim[2], ncol = nf)
Cmat <- matrix(runif(mydim[3]*nf), nrow = mydim[3], ncol = nf)
Xmat <- Emat <- Amat <- vector("list", mydim[3])
for(k in 1:mydim[3]){
  Amat[[k]] <- matrix(rnorm(nk[k]*nf), nrow = nk[k], ncol = nf)
  Amat[[k]] <- svd(Amat[[k]], nv = 0)$u \%*\% Gmat
  Xmat[[k]] <- tcrossprod(Amat[[k]] \%*\% diag(Cmat[k,]), Bmat)
  Emat[[k]] <- matrix(rnorm(nk[k]*mydim[2]), nrow = nk[k], ncol = mydim[2])
}
Emat <- nscale(Emat, 0, ssnew = sumsq(Xmat))   # SNR = 1
X <- mapply("+", Xmat, Emat)

# fit Parafac2 model (unconstrained)
pfac <- parafac2(X, nfac = nf, nstart = 1)
pfac

# check solution
Xhat <- fitted(pfac)
sse <- sumsq(mapply("-", Xmat, Xhat))
sse / (sum(nk) * mydim[2])
crossprod(pfac$A[[1]])
crossprod(pfac$A[[2]])
pfac$Phi

# reorder and resign factors
pfac$B[1:4,]
pfac <- reorder(pfac, 2:1)
pfac$B[1:4,]
pfac <- resign(pfac, mode="B")
pfac$B[1:4,]
Xhat <- fitted(pfac)
sse <- sumsq(mapply("-", Xmat, Xhat))
sse / (sum(nk) * mydim[2])

# rescale factors
colSums(pfac$B^2)
colSums(pfac$C^2)
pfac <- rescale(pfac, mode = "C", absorb = "B")
colSums(pfac$B^2)
colSums(pfac$C^2)
Xhat <- fitted(pfac)
sse <- sumsq(mapply("-", Xmat, Xhat))
sse / (sum(nk) * mydim[2])


##########   4-way example   ##########

# create random data list with Parafac2 structure
set.seed(4)
mydim <- c(NA, 10, 20, 5)
nf <- 3
nk <- rep(c(50,100,200), length.out = mydim[4])
Gmat <- matrix(rnorm(nf^2), nrow = nf, ncol = nf)
Bmat <- scale(matrix(rnorm(mydim[2]*nf), nrow = mydim[2], ncol = nf), center = FALSE)
cseq <- seq(-3, 3, length=mydim[3])
Cmat <- cbind(pnorm(cseq), pgamma(cseq+3.1, shape=1, rate=1)*(3/4), pt(cseq-2, df=4)*2)
Dmat <- scale(matrix(runif(mydim[4]*nf)*2, nrow = mydim[4], ncol = nf), center = FALSE)
Xmat <- Emat <- Amat <- vector("list",mydim[4])
for(k in 1:mydim[4]){
  aseq <- seq(-3, 3, length.out = nk[k])
  Amat[[k]] <- cbind(sin(aseq), sin(abs(aseq)), exp(-aseq^2))
  Amat[[k]] <- svd(Amat[[k]], nv = 0)$u \%*\% Gmat
  Xmat[[k]] <- array(tcrossprod(Amat[[k]] \%*\% diag(Dmat[k,]),
                                krprod(Cmat, Bmat)), dim = c(nk[k], mydim[2], mydim[3]))
  Emat[[k]] <- array(rnorm(nk[k] * mydim[2] * mydim[3]), dim = c(nk[k], mydim[2], mydim[3]))
}
Emat <- nscale(Emat, 0, ssnew = sumsq(Xmat))   # SNR = 1
X <- mapply("+", Xmat, Emat)

# fit Parafac model (smooth A, unconstrained B, monotonic C, non-negative D)
pfac <- parafac2(X, nfac = nf, nstart = 1, 
                 const = c("smooth", "uncons", "moninc", "nonneg"))
pfac

# check solution
Xhat <- fitted(pfac)
sse <- sumsq(mapply("-", Xmat, Xhat))
sse / (sum(nk) * mydim[2] * mydim[3])
crossprod(pfac$A[[1]])
crossprod(pfac$A[[2]])
pfac$Phi


\dontrun{

##########   parallel computation   ##########

# create random data list with Parafac2 structure
set.seed(3)
mydim <- c(NA, 10, 20)
nf <- 2
nk <- rep(c(50, 100, 200), length.out = mydim[3])
Gmat <- matrix(rnorm(nf^2), nrow = nf, ncol = nf)
Bmat <- matrix(runif(mydim[2]*nf), nrow = mydim[2], ncol = nf)
Cmat <- matrix(runif(mydim[3]*nf), nrow = mydim[3], ncol = nf)
Xmat <- Emat <- Hmat <- vector("list", mydim[3])
for(k in 1:mydim[3]){
  Hmat[[k]] <- svd(matrix(rnorm(nk[k] * nf), nrow = nk[k], ncol = nf), nv = 0)$u
  Xmat[[k]] <- tcrossprod(Hmat[[k]] \%*\% Gmat \%*\% diag(Cmat[k,]), Bmat)
  Emat[[k]] <- matrix(rnorm(nk[k] * mydim[2]), nrow = nk[k], mydim[2])
}
Emat <- nscale(Emat, 0, ssnew = sumsq(Xmat))   # SNR = 1
X <- mapply("+", Xmat, Emat)

# fit Parafac2 model (10 random starts -- sequential computation)
set.seed(1)
system.time({pfac <- parafac2(X, nfac = nf)})
pfac

# fit Parafac2 model (10 random starts -- parallel computation)
cl <- makeCluster(detectCores())
ce <- clusterEvalQ(cl, library(multiway))
clusterSetRNGStream(cl, 1)
system.time({pfac <- parafac2(X, nfac = nf, parallel = TRUE, cl = cl)})
pfac
stopCluster(cl)
}

}

\keyword{ models }
\keyword{ multivariate }
