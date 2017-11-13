\name{parafac2}
\alias{parafac2}
\title{
Parallel Factor Analysis-2
}
\description{
Given a list of matrices \code{X[[k]] = matrix(xk,I[k],J)} for \code{k = seq(1,K)}, the 3-way Parafac2 model (with Mode A nested in Mode C) can be written as
\tabular{r}{
\code{X[[k]] = tcrossprod(A[[k]]\%*\%diag(C[k,]),B) + E[[k]]} \cr 
subject to \code{crossprod(A[[k]]) = Phi} \cr 
}
where \code{A[[k]] = matrix(ak,I[k],R)} are the Mode A (first mode) weights for the \code{k}-th level of Mode C (third mode), \code{Phi} is the common crossproduct matrix shared by all \code{K} levels of Mode C, \code{B = matrix(b,J,R)} are the Mode B (second mode) weights, \code{C = matrix(c,K,R)} are the Mode C (third mode) weights, and \code{E[[k]] = matrix(ek,I[k],J)} is the residual matrix corresponding to \code{k}-th level of Mode C.

Given a list of arrays \code{X[[l]] = array(xl,dim=c(I[l],J,K))} for \code{l = seq(1,L)}, the 4-way Parafac2 model (with Mode A nested in Mode D) can be written as
\tabular{r}{
\code{X[[l]][,,k] = tcrossprod(A[[l]]\%*\%diag(D[l,]*C[k,]),B) + E[[k]]} \cr 
subject to \code{crossprod(A[[l]]) = Phi} \cr 
}
\code{A[[l]] = matrix(al,I[l],R)} are the Mode A (first mode) weights for the \code{l}-th level of Mode D (fourth mode), \code{Phi} is the common crossproduct matrix shared by all \code{L} levels of Mode D, \code{D = matrix(d,L,R)} are the Mode D (fourth mode) weights, and \code{E[[l]] = matrix(el,I[l],J,K)} is the residual array corresponding to \code{l}-th level of Mode D.

Weight matrices are estimated using an alternating least squares algorithm with optional constraints.
}
\usage{
parafac2(X, nfac, nstart = 10, const = NULL, control = NULL,
         Gfixed = NULL, Bfixed = NULL, Cfixed = NULL, Dfixed = NULL,
         Gstart = NULL, Bstart = NULL, Cstart = NULL, Dstart = NULL,
         Gstruc = NULL, Bstruc = NULL, Cstruc = NULL, Dstruc = NULL,
         maxit = 500, ctol = 1e-4, parallel = FALSE, cl = NULL,
         output = c("best", "all"), verbose = FALSE)
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
  Constraints for each mode. Vector of length 3 or 4 with entries:
  0 = unconstrained (default),
  1 = orthogonal,
  2 = non-negative,
  3 = unimodal,
  4 = monotonic,
  5 = periodic,
  6 = smooth.
  Use \code{control} argument to adjust options for constraints 3-6.
  Note: constraints 2-4 cannot be applied on Mode A.
}
\item{control}{
  List of parameters controlling options for constraints 3-6. This is passed to \code{\link{const.control}}, which describes the available options.
}
  \item{Gfixed}{
  Fixed Mode A crossproducts (\code{crossprod(Gfixed)=Phi}). Only used to fit model with fixed Phi matrix.
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
  \item{Gstart}{
  Starting Mode A crossproduct matrix for ALS algorithm (\code{crossprod(Gstart)=Phi}). Default uses random weights.
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
  \item{Gstruc}{
  Structure constraints for Mode A crossproduct matrix (\code{crossprod(Gstruc) = Phi structure}). Default uses unstructured crossproducts.
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
  Convergence tolerance.
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
If \code{output="best"}, returns an object of class \code{"parafac2"} with the following elements:
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
\item{cflag}{Convergence flag.}
\item{const}{See argument \code{const}.}
\item{control}{See argument \code{control}.}
\item{fixed}{Logical vector indicating whether 'fixed' weights were used for each mode.}
\item{struc}{Logical vector indicating whether 'struc' constraints were used for each mode.}

Otherwise returns a list of length \code{nstart} where each element is an object of class \code{"parafac2"}.
}
\references{
Bro, R., & De Jong, S. (1997). A fast non-negativity-constrained least squares algorithm. \emph{Journal of Chemometrics, 11}, 393-401.

Harshman, R. A. (1972). PARAFAC2: Mathematical and technical notes. \emph{UCLA Working Papers in Phonetics, 22}, 30-44.

Helwig, N. E. (2013). The special sign indeterminacy of the direct-fitting Parafac2 model: Some implications, cautions, and recommendations, for Simultaneous Component Analysis. \emph{Psychometrika, 78}, 725-739.

Helwig, N. E. (2017). Estimating latent trends in multivariate longitudinal data via Parafac2 with functional and structural constraints. \emph{Biometrical Journal, 59}(4), 783-803.

Kiers, H. A. L., ten Berge, J. M. F., & Bro, R. (1999). PARAFAC2-part I: A direct-fitting algorithm for the PARAFAC2 model. \emph{Journal of Chemometrics, 13}, 275-294.

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

# create random data list with Parafac2 structure
set.seed(3)
mydim <- c(NA,10,20)
nf <- 2
nk <- sample(c(50,100,200),mydim[3],replace=TRUE)
Gmat <- matrix(rnorm(nf^2),nf,nf)
Bmat <- matrix(runif(mydim[2]*nf),mydim[2],nf)
Cmat <- matrix(runif(mydim[3]*nf),mydim[3],nf)
Xmat <- Emat <- Hmat <- vector("list",mydim[3])
for(k in 1:mydim[3]){
  Hmat[[k]] <- svd(matrix(rnorm(nk[k]*nf),nk[k],nf),nv=0)$u
  Xmat[[k]] <- tcrossprod(Hmat[[k]]\%*\%Gmat\%*\%diag(Cmat[k,]),Bmat)
  Emat[[k]] <- matrix(rnorm(nk[k]*mydim[2]),nk[k],mydim[2])
}
Emat <- nscale(Emat,0,sumsq(Xmat))   # SNR=1
X <- mapply("+",Xmat,Emat)

# fit Parafac2 model (unconstrained)
pfac <- parafac2(X,nfac=nf,nstart=1)
pfac

# check solution
Xhat <- fitted(pfac)
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2])
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
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2])

# rescale factors
colSums(pfac$B^2)
colSums(pfac$C^2)
pfac <- rescale(pfac, mode="C", absorb="B")
colSums(pfac$B^2)
colSums(pfac$C^2)
Xhat <- fitted(pfac)
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2])


##########   4-way example   ##########

# create random data list with Parafac2 structure
set.seed(4)
mydim <- c(NA,10,20,5)
nf <- 3
nk <- sample(c(50,100,200),mydim[4],replace=TRUE)
Gmat <- matrix(rnorm(nf^2),nf,nf)
Bmat <- scale(matrix(rnorm(mydim[2]*nf),mydim[2],nf), center=FALSE)
cseq <- seq(-3, 3, length=mydim[3])
Cmat <- cbind(pnorm(cseq), pgamma(cseq+3.1, shape=1, rate=1)*(3/4), pt(cseq-2, df=4)*2)
Dmat <- scale(matrix(runif(mydim[4]*nf)*2,mydim[4],nf), center=FALSE)
Xmat <- Emat <- Hmat <- vector("list",mydim[4])
for(k in 1:mydim[4]){
  aseq <- seq(-3,3,length=nk[k])
  Hmat[[k]] <- svd(cbind(sin(aseq), sin(abs(aseq)), exp(-aseq^2)),nv=0)$u
  Xmat[[k]] <- array(tcrossprod(Hmat[[k]]\%*\%Gmat\%*\%diag(Dmat[k,]),
                             krprod(Cmat,Bmat)),dim=c(nk[k],mydim[2],mydim[3]))
  Emat[[k]] <- array(rnorm(nk[k]*mydim[2]*mydim[3]),dim=c(nk[k],mydim[2],mydim[3]))
}
Emat <- nscale(Emat,0,sumsq(Xmat))   # SNR=1
X <- mapply("+",Xmat,Emat)

# fit Parafac model (smooth A, unconstrained B, monotonic C, non-negative D)
pfac <- parafac2(X,nfac=nf,nstart=1,const=c(6,0,4,2))
pfac

# check solution
Xhat <- fitted(pfac)
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2]*mydim[3])
crossprod(pfac$A[[1]])
crossprod(pfac$A[[2]])
pfac$Phi


\dontrun{

##########   parallel computation   ##########

# create random data list with Parafac2 structure
set.seed(3)
mydim <- c(NA,10,20)
nf <- 2
nk <- sample(c(50,100,200),mydim[3],replace=TRUE)
Gmat <- matrix(rnorm(nf^2),nf,nf)
Bmat <- matrix(runif(mydim[2]*nf),mydim[2],nf)
Cmat <- matrix(runif(mydim[3]*nf),mydim[3],nf)
Xmat <- Emat <- Hmat <- vector("list",mydim[3])
for(k in 1:mydim[3]){
  Hmat[[k]] <- svd(matrix(rnorm(nk[k]*nf),nk[k],nf),nv=0)$u
  Xmat[[k]] <- tcrossprod(Hmat[[k]]\%*\%Gmat\%*\%diag(Cmat[k,]),Bmat)
  Emat[[k]] <- matrix(rnorm(nk[k]*mydim[2]),nk[k],mydim[2])
}
Emat <- nscale(Emat,0,sumsq(Xmat))   # SNR=1
X <- mapply("+",Xmat,Emat)

# fit Parafac2 model (10 random starts -- sequential computation)
set.seed(1)
system.time({pfac <- parafac2(X,nfac=nf)})
pfac

# fit Parafac2 model (10 random starts -- parallel computation)
set.seed(1)
cl <- makeCluster(detectCores())
ce <- clusterEvalQ(cl,library(multiway))
system.time({pfac <- parafac2(X,nfac=nf,parallel=TRUE,cl=cl)})
pfac
stopCluster(cl)
}

}
