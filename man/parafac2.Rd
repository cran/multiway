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
parafac2(X,nfac,nstart=10,const=NULL,maxit=500,
         Bfixed=NULL,Cfixed=NULL,Dfixed=NULL,
         Bstart=NULL,Cstart=NULL,Dstart=NULL,
         ctol=10^-7,parallel=FALSE,cl=NULL)
}
\arguments{
  \item{X}{
  For 3-way Parafac2: list of length \code{K} where \code{k}-th element is \code{I[k]}-by-\code{J} matrix or three-way data array with \code{dim=c(I,J,K)}. For 4-way Parafac2: list of length \code{L} where \code{l}-th element is \code{I[l]}-by-\code{J}-by-\code{K} array or four-way data array with \code{dim=c(I,J,K,L)}.
}
  \item{nfac}{
  Number of factors.
}
  \item{nstart}{
  Number of random starts.
}
  \item{const}{
  Constraints for each mode. See Examples.
} 
  \item{maxit}{
  Maximum number of iterations.
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
  \item{ctol}{
  Convergence tolerance.
}
  \item{parallel}{
  Logical indicating if \code{\link{parLapply}} should be used. See Examples.
}
  \item{cl}{
  Cluster created by \code{\link{makeCluster}}. Only used when \code{parallel=TRUE}.
}
}
\value{
\item{A}{List with 2 elements (see Note).}
\item{B}{Mode B weight matrix.}
\item{C}{Mode C weight matrix.}
\item{D}{Mode D weight matrix.}
\item{Rsq}{R-squared value.}
\item{iter}{Number of iterations.}
\item{cflag}{Convergence flag.}
\item{const}{Same as input \code{const}.}
}
\references{
Bro, R., & De Jong, S. (1997). A fast non-negativity-constrained least squares algorithm. \emph{Journal of Chemometrics, 11}, 393-401.

Harshman, R. A. (1972). PARAFAC2: Mathematical and technical notes. \emph{UCLA Working Papers in Phonetics, 22}, 30-44.

Helwig, N. E. (2013). The special sign indeterminacy of the direct-fitting Parafac2 model: Some implications, cautions, and recommendations, for Simultaneous Component Analysis. \emph{Psychometrika, 78}, 725-739.

Kiers, H. A. L., ten Berge, J. M. F., & Bro, R. (1999). PARAFAC2-part I: A direct-fitting algorithm for the PARAFAC2 model. \emph{Journal of Chemometrics, 13}, 275-294.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Output: \code{A$H} is list of orthogonal matrices, and \code{A$G} is \code{R}-by-\code{R} matrix such that \code{crossprod(A$G) = Phi}. For 3-way case \code{A$H[[k]]\%*\%A$G} gives Mode A weights for \code{k}-th level of Mode C. For 4-way case \code{A$H[[l]]\%*\%A$G} gives Mode A weights for \code{l}-th level of Mode D. 

Default use is 10 random strarts (\code{nstart=10}) with 500 maximum iterations of the ALS algorithm for each start (\code{maxit=500}) using a relative convergence tolerace of 10^-7 (\code{ctol=10^-7}).

Input \code{const} should be a three or four element integer vector. Set \code{const[j]=0} for unconstrained update in j-th mode weight matrix, \code{const[j]=1} for orthogonal update in j-th mode weight matrix, or \code{const[j]=2} for non-negative update in j-th mode. Default is unconstrained update for all modes. Note: non-negativity in Mode A is not allowed.

Output \code{cflag} gives convergence information: \code{cflag=0} if ALS algorithm converged normally, \code{cflag=1} if maximum iteration limit was reached before convergence, and \code{cflag=2} if ALS algorithm terminated abnormally due to problem with non-negativity constraints.
}
\section{Warnings }{
The ALS algorithm can perform poorly if the number of factors \code{nfac} is set too large.

Non-negativity constraints can be sensitive to local optima. 

Non-negativity constraints can result in slower performance.
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
Emat <- nrescale(Emat,0,sumsq(Xmat))   # SNR=1
X <- mapply("+",Xmat,Emat)

# fit Parafac2 model (unconstrained)
pfac <- parafac2(X,nfac=nf,nstart=1)
pfac$Rsq

# check solution
Xhat <- fitted(pfac)
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2])
crossprod(pfac$A$H[[1]])
crossprod(pfac$A$G)


##########   4-way example   ##########

# create random data list with Parafac2 structure
set.seed(4)
mydim <- c(NA,10,20,5)
nf <- 3
nk <- sample(c(50,100,200),mydim[4],replace=TRUE)
Gmat <- matrix(rnorm(nf^2),nf,nf)
Bmat <- matrix(runif(mydim[2]*nf),mydim[2],nf)
Cmat <- matrix(runif(mydim[3]*nf),mydim[3],nf)
Dmat <- matrix(runif(mydim[4]*nf),mydim[4],nf)
Xmat <- Emat <- Hmat <- vector("list",mydim[4])
for(k in 1:mydim[4]){
  Hmat[[k]] <- svd(matrix(rnorm(nk[k]*nf),nk[k],nf),nv=0)$u
  Xmat[[k]] <- array(tcrossprod(Hmat[[k]]\%*\%Gmat\%*\%diag(Dmat[k,]),
                             krprod(Cmat,Bmat)),dim=c(nk[k],mydim[2],mydim[3]))
  Emat[[k]] <- array(rnorm(nk[k]*mydim[2]*mydim[3]),dim=c(nk[k],mydim[2],mydim[3]))
}
Emat <- nrescale(Emat,0,sumsq(Xmat))   # SNR=1
X <- mapply("+",Xmat,Emat)

# fit Parafac2 model (unconstrained)
pfac <- parafac2(X,nfac=nf,nstart=1)
pfac$Rsq

# check solution
Xhat <- fitted(pfac)
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2]*mydim[3])
crossprod(pfac$A$H[[1]])
crossprod(pfac$A$G)


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
Emat <- nrescale(Emat,0,sumsq(Xmat))   # SNR=1
X <- mapply("+",Xmat,Emat)

# fit Parafac2 model (10 random starts -- sequential computation)
set.seed(1)
system.time({pfac <- parafac2(X,nfac=nf)})
pfac$Rsq

# fit Parafac2 model (10 random starts -- parallel computation)
set.seed(1)
cl <- makeCluster(detectCores())
ce <- clusterEvalQ(cl,library(multiway))
system.time({pfac <- parafac2(X,nfac=nf,parallel=TRUE,cl=cl)})
pfac$Rsq
stopCluster(cl)
}

}
