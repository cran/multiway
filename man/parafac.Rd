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
parafac(X,nfac,nstart=10,const=NULL,maxit=500,
        Bfixed=NULL,Cfixed=NULL,Dfixed=NULL,
        Bstart=NULL,Cstart=NULL,Dstart=NULL,
        ctol=10^-4,parallel=FALSE,cl=NULL,
        output=c("best","all"))
}
\arguments{
  \item{X}{
  Three-way data array with \code{dim=c(I,J,K)} or four-way data array with \code{dim=c(I,J,K,L)}.
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
}
\value{
If \code{output="best"}, returns an object of class \code{"parafac"} with the following elements:
\item{A}{Mode A weight matrix.}
\item{B}{Mode B weight matrix.}
\item{C}{Mode C weight matrix.}
\item{D}{Mode D weight matrix.}
\item{Rsq}{R-squared value.}
\item{iter}{Number of iterations.}
\item{cflag}{Convergence flag.}
\item{const}{Same as input \code{const}.}

Otherwise returns a list of length \code{nstart} where each element is an object of class \code{"parafac"}.
}
\references{
Bro, R., & De Jong, S. (1997). A fast non-negativity-constrained least squares algorithm. \emph{Journal of Chemometrics, 11}, 393-401.

Harshman, R. A. (1970). Foundations of the PARAFAC procedure: Models and conditions for an "explanatory" multimodal factor analysis. \emph{UCLA Working Papers in Phonetics, 16}, 1-84.

Harshman, R. A., & Lundy, M. E. (1994). PARAFAC: Parallel factor analysis. \emph{Computational Statistics and Data Analysis, 18}, 39-72.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Default use is 10 random strarts (\code{nstart=10}) with 500 maximum iterations of the ALS algorithm for each start (\code{maxit=500}) using a convergence tolerance of 10^-4 (\code{ctol=10^-4}). The algorithm is determined to have converged once the change in R^2 is less than or equal to \code{ctol}.

Input \code{const} should be a three or four element integer vector. Set \code{const[j]=0} for unconstrained update in j-th mode weight matrix, \code{const[j]=1} for orthogonal update in j-th mode weight matrix, or \code{const[j]=2} for non-negative update in j-th mode. Default is unconstrained update for all modes.

Output \code{cflag} gives convergence information: \code{cflag=0} if ALS algorithm converged normally, \code{cflag=1} if maximum iteration limit was reached before convergence, and \code{cflag=2} if ALS algorithm terminated abnormally due to problem with non-negativity constraints.
}
\section{Warnings }{
The ALS algorithm can perform poorly if the number of factors \code{nfac} is set too large.

Non-negativity constraints can be sensitive to local optima. 

Non-negativity constraints can result in slower performance.
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
pfac$Rsq

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
Amat <- matrix(rnorm(mydim[1]*nf),mydim[1],nf)
Bmat <- svd(matrix(runif(mydim[2]*nf),mydim[2],nf))$u
Cmat <- matrix(runif(mydim[3]*nf),mydim[3],nf)
Dmat <- matrix(runif(mydim[4]*nf),mydim[4],nf)
Xmat <- array(tcrossprod(Amat,krprod(Dmat,krprod(Cmat,Bmat))),dim=mydim)
Emat <- array(rnorm(prod(mydim)),dim=mydim)
Emat <- nscale(Emat,0,sumsq(Xmat))   # SNR=1
X <- Xmat + Emat

# fit Parafac model (unconstrained A, orthogonal B, non-negative C, unconstrained D)
pfac <- parafac(X,nfac=nf,nstart=1,const=c(0,1,2,0))
pfac$Rsq

# check solution
Xhat <- fitted(pfac)
sum((Xmat-Xhat)^2)/prod(mydim)
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
pfac$Rsq

# fit Parafac model (10 random starts -- parallel computation)
set.seed(1)
cl <- makeCluster(detectCores())
ce <- clusterEvalQ(cl,library(multiway))
system.time({pfac <- parafac(X,nfac=nf,parallel=TRUE,cl=cl)})
pfac$Rsq
stopCluster(cl)
}

}
