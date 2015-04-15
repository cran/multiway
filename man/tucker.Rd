\name{tucker}
\alias{tucker}
\title{
Tucker Factor Analysis
}
\description{
Given a 3-way array \code{X = array(x,dim=c(I,J,K))}, the 3-way Tucker model can be written as 
\tabular{c}{
\code{ X[i,j,k] = sum sum sum A[i,p]*B[j,q]*C[k,r]*G[p,q,r] + E[i,j,k] }
} 
where \code{A = matrix(a,I,P)} are the Mode A (first mode) weights, \code{B = matrix(b,J,Q)} are the Mode B (second mode) weights, \code{C = matrix(c,K,R)} are the Mode C (third mode) weights, \code{G = array(g,dim=c(P,Q,R))} is the 3-way core array, and \code{E = array(e,dim=c(I,J,K))} is the 3-way residual array. The summations are for \code{p = seq(1,P)}, \code{q = seq(1,Q)}, and \code{r = seq(1,R)}.

Given a 4-way array \code{X = array(x,dim=c(I,J,K,L))}, the 4-way Tucker model can be written as 
\tabular{c}{
\code{ X[i,j,k,l] = sum sum sum sum A[i,p]*B[j,q]*C[k,r]*D[l,s]*G[p,q,r,s] + E[i,j,k,l] }
} 
where \code{D = matrix(d,L,S)} are the Mode D (fourth mode) weights, \code{G = array(g,dim=c(P,Q,R,S))} is the 4-way residual array, \code{E = array(e,dim=c(I,J,K,L))} is the 4-way residual array, and the other terms can be interprered as previously described.

Weight matrices are estimated using an alternating least squares algorithm.

}
\usage{
tucker(X,nfac,nstart=10,maxit=500,Afixed=NULL,
       Bfixed=NULL,Cfixed=NULL,Dfixed=NULL,
       Bstart=NULL,Cstart=NULL,Dstart=NULL,
       ctol=10^-7,parallel=FALSE,cl=NULL)
}
\arguments{
  \item{X}{
  Three-way data array with \code{dim=c(I,J,K)} or four-way data array with \code{dim=c(I,J,K,L)}.
}
  \item{nfac}{
  Number of factors in each mode.
}
  \item{nstart}{
  Number of random starts.
}
  \item{maxit}{
  Maximum number of iterations.
}
  \item{Afixed}{
  Fixed Mode A weights. Only used to fit model with fixed weights in Mode A.
}
  \item{Bfixed}{
  Fixed Mode B weights. Only used to fit model with fixed weights in Mode B.
}
  \item{Cfixed}{
  Fixed Mode C weights. Only used to fit model with fixed weights in Mode C.
}
  \item{Dfixed}{
  Fixed Mode D weights. Only used to fit model with fixed weights in Mode D.
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
\item{A}{Mode A weight matrix.}
\item{B}{Mode B weight matrix.}
\item{C}{Mode C weight matrix.}
\item{D}{Mode D weight matrix.}
\item{G}{Core array.}
\item{Rsq}{R-squared value.}
\item{iter}{Number of iterations.}
\item{cflag}{Convergence flag.}
}
\references{
Kroonenberg, P. M., & de Leeuw, J. (1980). Principal component analysis of three-mode data by means of alternating least squares algorithms. \emph{Psychometrika, 45}, 69-97.
 
Tucker, L. R. (1966). Some mathematical notes on three-mode factor analysis. \emph{Psychometrika, 31}, 279-311.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Default use is 10 random strarts (\code{nstart=10}) with 500 maximum iterations of the ALS algorithm for each start (\code{maxit=500}) using a relative convergence tolerace of 10^-7 (\code{ctol=10^-7}).

Output \code{cflag} gives convergence information: \code{cflag=0} if ALS algorithm converged normally, and \code{cflag=1} if maximum iteration limit was reached before convergence.
}
\section{Warnings }{
The ALS algorithm can perform poorly if the number of factors \code{nfac} is set too large.

Input matrices in \code{Afixed}, \code{Bfixed}, \code{Cfixed}, \code{Dfixed}, \code{Bstart}, \code{Cstart},  and \code{Dstart} must be columnwise orthonormal.
}

\examples{
##########   3-way example   ##########

# create random data array with Tucker structure
set.seed(3)
mydim <- c(50,20,5)
nf <- c(3,2,3)
Amat <- svd(matrix(rnorm(mydim[1]*nf[1]),mydim[1],nf[1]),nu=nf[1])$u
Bmat <- svd(matrix(rnorm(mydim[2]*nf[2]),mydim[2],nf[2]),nu=nf[2])$u
Cmat <- svd(matrix(rnorm(mydim[3]*nf[3]),mydim[3],nf[3]),nu=nf[3])$u
Gmat <- array(rnorm(prod(nf)),dim=nf)
Xmat <- array(tcrossprod(Amat\%*\%matrix(Gmat,nf[1],nf[2]*nf[3]),kronecker(Cmat,Bmat)),dim=mydim)
Emat <- array(rnorm(prod(mydim)),dim=mydim)
Emat <- nrescale(Emat,0,sumsq(Xmat))   # SNR=1
X <- Xmat + Emat

# fit Tucker model
tuck <- tucker(X,nfac=nf,nstart=1)
tuck$Rsq

# check solution
Xhat <- fitted(tuck)
sum((Xmat-Xhat)^2)/prod(mydim)


##########   4-way example   ##########

# create random data array with Tucker structure
set.seed(4)
mydim <- c(30,10,8,10)
nf <- c(2,3,4,3)
Amat <- svd(matrix(rnorm(mydim[1]*nf[1]),mydim[1],nf[1]),nu=nf[1])$u
Bmat <- svd(matrix(rnorm(mydim[2]*nf[2]),mydim[2],nf[2]),nu=nf[2])$u
Cmat <- svd(matrix(rnorm(mydim[3]*nf[3]),mydim[3],nf[3]),nu=nf[3])$u
Dmat <- svd(matrix(rnorm(mydim[4]*nf[4]),mydim[4],nf[4]),nu=nf[4])$u
Gmat <- array(rnorm(prod(nf)),dim=nf)
Xmat <- array(tcrossprod(Amat\%*\%matrix(Gmat,nf[1],prod(nf[2:4])),
                      kronecker(Dmat,kronecker(Cmat,Bmat))),dim=mydim)
Emat <- array(rnorm(prod(mydim)),dim=mydim)
Emat <- nrescale(Emat,0,sumsq(Xmat))   # SNR=1
X <- Xmat + Emat

# fit Tucker model
tuck <- tucker(X,nfac=nf,nstart=1)
tuck$Rsq

# check solution
Xhat <- fitted(tuck)
sum((Xmat-Xhat)^2)/prod(mydim)


\dontrun{

##########   parallel computation   ##########

# create random data array with Tucker structure
set.seed(3)
mydim <- c(50,20,5)
nf <- c(3,2,3)
Amat <- svd(matrix(rnorm(mydim[1]*nf[1]),mydim[1],nf[1]),nu=nf[1])$u
Bmat <- svd(matrix(rnorm(mydim[2]*nf[2]),mydim[2],nf[2]),nu=nf[2])$u
Cmat <- svd(matrix(rnorm(mydim[3]*nf[3]),mydim[3],nf[3]),nu=nf[3])$u
Gmat <- array(rnorm(prod(nf)),dim=nf)
Xmat <- array(tcrossprod(Amat\%*\%matrix(Gmat,nf[1],nf[2]*nf[3]),kronecker(Cmat,Bmat)),dim=mydim)
Emat <- array(rnorm(prod(mydim)),dim=mydim)
Emat <- nrescale(Emat,0,sumsq(Xmat))   # SNR=1
X <- Xmat + Emat

# fit Tucker model (10 random starts -- sequential computation)
set.seed(1)
system.time({tuck <- tucker(X,nfac=nf)})
tuck$Rsq

# fit Tucker model (10 random starts -- parallel computation)
set.seed(1)
cl <- makeCluster(detectCores())
ce <- clusterEvalQ(cl,library(multiway))
system.time({tuck <- tucker(X,nfac=nf,parallel=TRUE,cl=cl)})
tuck$Rsq
stopCluster(cl)
}

}
