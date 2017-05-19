\name{sca}
\alias{sca}
\title{
Simultaneous Component Analysis
}
\description{
Given a list of matrices \code{X[[k]] = matrix(xk,I[k],J)} for \code{k = seq(1,K)}, the SCA model is 
\tabular{c}{
X[[k]] = tcrossprod(D[[k]],B) + E[[k]]
}
where \code{D[[k]] = matrix(dk,I[k],R)} are the Mode A (first mode) weights for the \code{k}-th level of Mode C (third mode), \code{B = matrix(b,J,R)} are the Mode B (second mode) weights, and \code{E[[k]] = matrix(ek,I[k],J)} is the residual matrix corresponding to \code{k}-th level of Mode C. 

There are four different versions of the SCA model: SCA with invariant pattern (SCA-P), SCA with Parafac2 constraints (SCA-PF2), SCA with INDSCAL constraints (SCA-IND), and SCA with equal average crossproducts (SCA-ECP). These four models differ with respect to the assumed crossproduct structure of the \code{D[[k]]} weights:
\tabular{rcl}{
  SCA-P: \tab \tab \code{crossprod(D[[k]])/I[k] = Phi[[k]]} \cr
  SCA-PF2: \tab \tab \code{crossprod(D[[k]])/I[k] = diag(C[k,])\%*\%Phi\%*\%diag(C[k,])} \cr
  SCA-IND: \tab \tab \code{crossprod(D[[k]])/I[k] = diag(C[k,]*C[k,])} \cr
  SCA-ECP: \tab \tab \code{crossprod(D[[k]])/I[k] = Phi} \cr
  }
where \code{Phi[[k]]} is specific to the \code{k}-th level of Mode C, \code{Phi} is common to all \code{K} levels of Mode C, and \code{C = matrix(c,K,R)} are the Mode C (third mode) weights. This function estimates the weight matrices \code{D[[k]]} and \code{B} (and \code{C} if applicable) using alternating least squares. 
}
\usage{
sca(X,nfac,nstart=10,maxit=500,
    type=c("sca-p","sca-pf2","sca-ind","sca-ecp"),
    rotation=c("none","varimax","promax"),
    ctol=1e-4,parallel=FALSE,cl=NULL)
}
\arguments{
  \item{X}{
  List of length \code{K} where the \code{k}-th element contains the \code{I[k]}-by-\code{J} data matrix \code{X[[k]]}. If \code{I[k]=I[1]} for all \code{k}, can input 3-way data array with \code{dim=c(I,J,K)}.
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
  \item{type}{
  Type of SCA model to fit. 
}
  \item{rotation}{
  Rotation to use for \code{type="sca-p"} or \code{type="sca-ecp"}.
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
\item{D}{List of length \code{K} where \code{k}-th element contains \code{D[[k]]}.}
\item{B}{Mode B weight matrix.}
\item{C}{Mode C weight matrix.}
\item{Phi}{Mode A common crossproduct matrix (if \code{type!="sca-p"}).}
\item{SSE}{Sum of Squared Errors.}
\item{Rsq}{R-squared value.}
\item{GCV}{Generalized Cross-Validation.}
\item{edf}{Effective degrees of freedom.}
\item{iter}{Number of iterations.}
\item{cflag}{Convergence flag.}
\item{type}{Same as input \code{type}.}
\item{rotation}{Same as input \code{rotation}.}
}
\references{
Helwig, N. E. (2013). The special sign indeterminacy of the direct-fitting Parafac2 model: Some implications, cautions, and recommendations, for Simultaneous Component Analysis. \emph{Psychometrika, 78}, 725-739.

Timmerman, M. E., & Kiers, H. A. L. (2003). Four simultaneous component models for the analysis of multivariate time series from more than one subject to model intraindividual and interindividual differences. \emph{Psychometrika, 68}, 105-121.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Default use is 10 random strarts (\code{nstart=10}) with 500 maximum iterations of the ALS algorithm for each start (\code{maxit=500}) using a convergence tolerance of 1e-4 (\code{ctol=1e-4}). The algorithm is determined to have converged once the change in R^2 is less than or equal to \code{ctol}.

Output \code{cflag} gives convergence information: \code{cflag=0} if ALS algorithm converged normally, \code{cflag=1} if maximum iteration limit was reached before convergence, and \code{cflag=2} if ALS algorithm terminated abnormally due to problem with non-negativity constraints.
}
\section{Warnings }{
The ALS algorithm can perform poorly if the number of factors \code{nfac} is set too large.
}
\section{Computational Details }{
  The least squares SCA-P solution can be obtained from the singular value decomposition of the stacked matrix \code{rbind(X[[1]],...,X[[K]])}.
  
  The least squares SCA-PF2 solution can be obtained using the uncontrained Parafac2 ALS algorithm (see \code{\link{parafac2}}).
  
  The least squares SCA-IND solution can be obtained using the Parafac2 ALS algorithm with orthogonality constraints on Mode A.
  
  The least squares SCA-ECP solution can be obtained using the Parafac2 ALS algorithm with orthogonality constraints on Mode A and the Mode C weights fixed at \code{C[k,] = rep(I[k]^0.5,R)}.
}

\examples{
##########   sca-p   ##########

# create random data list with SCA-P structure
set.seed(3)
mydim <- c(NA,10,20)
nf <- 2
nk <- sample(c(50,100,200),mydim[3],replace=TRUE)
Dmat <- matrix(rnorm(sum(nk)*nf),sum(nk),nf)
Bmat <- matrix(runif(mydim[2]*nf),mydim[2],nf)
Dmats <- vector("list",mydim[3])
Xmat <- Emat <- vector("list",mydim[3])
dfc <- 0
for(k in 1:mydim[3]){
  dinds <- 1:nk[k] + dfc
  Dmats[[k]] <- Dmat[dinds,]
  dfc <- dfc + nk[k]
  Xmat[[k]] <- tcrossprod(Dmats[[k]],Bmat)
  Emat[[k]] <- matrix(rnorm(nk[k]*mydim[2]),nk[k],mydim[2])
}
rm(Dmat)
Emat <- nscale(Emat,0,sumsq(Xmat))   # SNR=1
X <- mapply("+",Xmat,Emat)

# fit SCA-P model (no rotation)
scamod <- sca(X,nfac=nf,nstart=1)
scamod

# check solution
crossprod(scamod$D[[1]] \%*\% diag(scamod$C[1,]^-1) ) / nk[1]
crossprod(scamod$D[[5]] \%*\% diag(scamod$C[5,]^-1) ) / nk[5]
Xhat <- fitted(scamod)
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2])

# reorder and resign factors
scamod$B[1:4,]
scamod <- reorder(scamod, 2:1)
scamod$B[1:4,]
scamod <- resign(scamod, mode="B", newsign=c(1,-1))
scamod$B[1:4,]
Xhat <- fitted(scamod)
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2])

# rescale factors
colSums(scamod$B^2)
colSums(scamod$C^2)
scamod <- rescale(scamod, mode="C")
colSums(scamod$B^2)
colSums(scamod$C^2)
Xhat <- fitted(scamod)
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2])


##########   sca-pf2   ##########

# create random data list with SCA-PF2 (Parafac2) structure
set.seed(3)
mydim <- c(NA,10,20)
nf <- 2
nk <- sample(c(50,100,200),mydim[3],replace=TRUE)
Gmat <- 10*matrix(rnorm(nf^2),nf,nf)
Bmat <- matrix(runif(mydim[2]*nf),mydim[2],nf)
Cmat <- matrix(runif(mydim[3]*nf),mydim[3],nf)
Xmat <- Emat <- Fmat <- vector("list",mydim[3])
for(k in 1:mydim[3]){
  Fmat[[k]] <- svd(matrix(rnorm(nk[k]*nf),nk[k],nf),nv=0)$u
  Xmat[[k]] <- tcrossprod(Fmat[[k]]\%*\%Gmat\%*\%diag(Cmat[k,]),Bmat)
  Emat[[k]] <- matrix(rnorm(nk[k]*mydim[2]),nk[k],mydim[2])
}
Emat <- nscale(Emat,0,sumsq(Xmat))   # SNR=1
X <- mapply("+",Xmat,Emat)

# fit SCA-PF2 model
scamod <- sca(X,nfac=nf,nstart=1,type="sca-pf2")
scamod

# check solution
scamod$Phi
crossprod(scamod$D[[1]] \%*\% diag(scamod$C[1,]^-1) ) / nk[1]
crossprod(scamod$D[[5]] \%*\% diag(scamod$C[5,]^-1) ) / nk[5]
Xhat <- fitted(scamod)
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2])

# reorder and resign factors
scamod$B[1:4,]
scamod <- reorder(scamod, 2:1)
scamod$B[1:4,]
scamod <- resign(scamod, mode="B", newsign=c(1,-1))
scamod$B[1:4,]
Xhat <- fitted(scamod)
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2])

# rescale factors
colSums(scamod$B^2)
colSums(scamod$C^2)
scamod <- rescale(scamod, mode="C")
colSums(scamod$B^2)
colSums(scamod$C^2)
Xhat <- fitted(scamod)
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2])


##########   sca-ind   ##########

# create random data list with SCA-IND structure
set.seed(3)
mydim <- c(NA,10,20)
nf <- 2
nk <- sample(c(50,100,200),mydim[3],replace=TRUE)
Gmat <- diag(nf)  # SCA-IND is Parafac2 with Gmat=identity
Bmat <- matrix(runif(mydim[2]*nf),mydim[2],nf)
Cmat <- 10*matrix(runif(mydim[3]*nf),mydim[3],nf)
Xmat <- Emat <- Fmat <- vector("list",mydim[3])
for(k in 1:mydim[3]){
  Fmat[[k]] <- svd(matrix(rnorm(nk[k]*nf),nk[k],nf),nv=0)$u
  Xmat[[k]] <- tcrossprod(Fmat[[k]]\%*\%Gmat\%*\%diag(Cmat[k,]),Bmat)
  Emat[[k]] <- matrix(rnorm(nk[k]*mydim[2]),nk[k],mydim[2])
}
Emat <- nscale(Emat,0,sumsq(Xmat))   # SNR=1
X <- mapply("+",Xmat,Emat)

# fit SCA-IND model
scamod <- sca(X,nfac=nf,nstart=1,type="sca-ind")
scamod

# check solution
scamod$Phi
crossprod(scamod$D[[1]] \%*\% diag(scamod$C[1,]^-1) ) / nk[1]
crossprod(scamod$D[[5]] \%*\% diag(scamod$C[5,]^-1) ) / nk[5]
Xhat <- fitted(scamod)
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2])

# reorder and resign factors
scamod$B[1:4,]
scamod <- reorder(scamod, 2:1)
scamod$B[1:4,]
scamod <- resign(scamod, mode="B", newsign=c(1,-1))
scamod$B[1:4,]
Xhat <- fitted(scamod)
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2])

# rescale factors
colSums(scamod$B^2)
colSums(scamod$C^2)
scamod <- rescale(scamod, mode="C")
colSums(scamod$B^2)
colSums(scamod$C^2)
Xhat <- fitted(scamod)
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2])


##########   sca-ecp   ##########

# create random data list with SCA-ECP structure
set.seed(3)
mydim <- c(NA,10,20)
nf <- 2
nk <- sample(c(50,100,200),mydim[3],replace=TRUE)
Gmat <- diag(nf)
Bmat <- matrix(runif(mydim[2]*nf),mydim[2],nf)
Cmat <- matrix(sqrt(nk),mydim[3],nf)
Xmat <- Emat <- Fmat <- vector("list",mydim[3])
for(k in 1:mydim[3]){
  Fmat[[k]] <- svd(matrix(rnorm(nk[k]*nf),nk[k],nf),nv=0)$u
  Xmat[[k]] <- tcrossprod(Fmat[[k]]\%*\%Gmat\%*\%diag(Cmat[k,]),Bmat)
  Emat[[k]] <- matrix(rnorm(nk[k]*mydim[2]),nk[k],mydim[2])
}
Emat <- nscale(Emat,0,sumsq(Xmat))   # SNR=1
X <- mapply("+",Xmat,Emat)

# fit SCA-ECP model
scamod <- sca(X,nfac=nf,nstart=1,type="sca-ecp")
scamod

# check solution
scamod$Phi
crossprod(scamod$D[[1]] \%*\% diag(scamod$C[1,]^-1) ) / nk[1]
crossprod(scamod$D[[5]] \%*\% diag(scamod$C[5,]^-1) ) / nk[5]
Xhat <- fitted(scamod)
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2])

# reorder and resign factors
scamod$B[1:4,]
scamod <- reorder(scamod, 2:1)
scamod$B[1:4,]
scamod <- resign(scamod, mode="B", newsign=c(-1,1))
scamod$B[1:4,]
Xhat <- fitted(scamod)
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2])

# rescale factors
colSums(scamod$B^2)
colSums(scamod$C^2)
scamod <- rescale(scamod, mode="B")
colSums(scamod$B^2)
colSums(scamod$C^2)
Xhat <- fitted(scamod)
sse <- sumsq(mapply("-",Xmat,Xhat))
sse/(sum(nk)*mydim[2])


\dontrun{

##########   parallel computation   ##########

# create random data list with SCA-IND structure
set.seed(3)
mydim <- c(NA,10,20)
nf <- 2
nk <- sample(c(50,100,200),mydim[3],replace=TRUE)
Gmat <- diag(nf)  # SCA-IND is Parafac2 with Gmat=identity
Bmat <- matrix(runif(mydim[2]*nf),mydim[2],nf)
Cmat <- 10*matrix(runif(mydim[3]*nf),mydim[3],nf)
Xmat <- Emat <- Fmat <- vector("list",mydim[3])
for(k in 1:mydim[3]){
  Fmat[[k]] <- svd(matrix(rnorm(nk[k]*nf),nk[k],nf),nv=0)$u
  Xmat[[k]] <- tcrossprod(Fmat[[k]]\%*\%Gmat\%*\%diag(Cmat[k,]),Bmat)
  Emat[[k]] <- matrix(rnorm(nk[k]*mydim[2]),nk[k],mydim[2])
}
Emat <- nscale(Emat,0,sumsq(Xmat))   # SNR=1
X <- mapply("+",Xmat,Emat)

# fit SCA-PF2 model (10 random starts -- sequential computation)
set.seed(1)
system.time({scamod <- sca(X,nfac=nf,type="sca-pf2")})
scamod

# fit SCA-PF2 model (10 random starts -- parallel computation)
set.seed(1)
cl <- makeCluster(detectCores())
ce <- clusterEvalQ(cl,library(multiway))
system.time({scamod <- sca(X,nfac=nf,type="sca-pf2",parallel=TRUE,cl=cl)})
scamod
stopCluster(cl)

# fit SCA-IND model (10 random starts -- sequential computation)
set.seed(1)
system.time({scamod <- sca(X,nfac=nf,type="sca-ind")})
scamod

# fit SCA-IND model (10 random starts -- parallel computation)
set.seed(1)
cl <- makeCluster(detectCores())
ce <- clusterEvalQ(cl,library(multiway))
system.time({scamod <- sca(X,nfac=nf,type="sca-ind",parallel=TRUE,cl=cl)})
scamod
stopCluster(cl)

# fit SCA-ECP model (10 random starts -- sequential computation)
set.seed(1)
system.time({scamod <- sca(X,nfac=nf,type="sca-ecp")})
scamod

# fit SCA-ECP model (10 random starts -- parallel computation)
set.seed(1)
cl <- makeCluster(detectCores())
ce <- clusterEvalQ(cl,library(multiway))
system.time({scamod <- sca(X,nfac=nf,type="sca-ecp",parallel=TRUE,cl=cl)})
scamod
stopCluster(cl)
}

}
