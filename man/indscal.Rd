\name{indscal}
\alias{indscal}
\title{
Individual Differences Scaling
}
\description{
Given a 3-way array \code{X = array(x,dim=c(J,J,K))} with \code{X[,,k]} denoting the \code{k}-th subject's dissimilarity matrix rating \code{J} objects, the INDSCAL model can be written as 
\tabular{c}{
\code{ Z[i,j,k] = sum B[i,r]*B[j,r]*C[k,r] + E[i,j,k] }
} 
where \code{Z} is the array of scalar products obtained from \code{X}, \code{B = matrix(b,J,R)} are the object weights, \code{C = matrix(c,K,R)} are the non-negative subject weights, and \code{E = array(e,dim=c(J,J,K))} is the 3-way residual array. The summation is for \code{r = seq(1,R)}.

Weight matrices are estimated using an alternating least squares algorithm with optional constraints.
}
\usage{
indscal(X,nfac,nstart=10,const=NULL,maxit=500,
        type=c("dissimilarity","similarity"),
        ctol=1e-4,parallel=FALSE,cl=NULL,
        output=c("best","all"))
}
\arguments{
  \item{X}{
  Three-way data array with \code{dim=c(J,J,K)} where \code{X[,,k]} is dissimilarity matrix. Can also input a list of (dis)similarity matrices or objects output by \code{\link{dist}}.
}
  \item{nfac}{
  Number of factors.
}
  \item{nstart}{
  Number of random starts.
}
  \item{const}{
  Constraints for Modes B and C. See Note.
}
  \item{maxit}{
  Maximum number of iterations.
}
  \item{type}{
  Character indicating if \code{X} contains dissimilarity data (default) or similarity data.
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
}
\value{
If \code{output="best"}, returns an object of class \code{"indscal"} with the following elements:
\item{B}{Mode B weight matrix.}
\item{C}{Mode C weight matrix.}
\item{SSE}{Sum of Squared Errors.}
\item{Rsq}{R-squared value.}
\item{GCV}{Generalized Cross-Validation.}
\item{edf}{Effective degrees of freedom.}
\item{iter}{Number of iterations.}
\item{cflag}{Convergence flag.}
\item{const}{See argument \code{const}.}

Otherwise returns a list of length \code{nstart} where each element is an object of class \code{"indscal"}.
}
\references{
Bro, R., & De Jong, S. (1997). A fast non-negativity-constrained least squares algorithm. \emph{Journal of Chemometrics, 11}, 393-401.

Carroll, J. D., & Chang, J-J. (1970). Analysis of individual differences in multidimensional scaling via an n-way generalization of "Eckart-Young" decomposition. \emph{Psychometrika, 35}, 283-319.

Harshman, R. A. (1970). Foundations of the PARAFAC procedure: Models and conditions for an "explanatory" multimodal factor analysis. \emph{UCLA Working Papers in Phonetics, 16}, 1-84.

Harshman, R. A., & Lundy, M. E. (1994). PARAFAC: Parallel factor analysis. \emph{Computational Statistics and Data Analysis, 18}, 39-72.
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Default use is 10 random strarts (\code{nstart=10}) with 500 maximum iterations of the ALS algorithm for each start (\code{maxit=500}) using a convergence tolerance of 1e-4 (\code{ctol=1e-4}). The algorithm is determined to have converged once the change in R^2 is less than or equal to \code{ctol}.

Input \code{const} should be a two element integer vector giving constraints for Modes B and C. There are four possible options:
\tabular{rl}{
\code{const=c(0,0)} \tab Unconstrained update for Modes B and C  \cr
\code{const=c(0,2)} \tab Unconstrained Mode B with non-negative Mode C \cr
\code{const=c(1,0)} \tab Orthogonal Mode B with unconstrained Mode C \cr
\code{const=c(1,2)} \tab Orthogonal Mode B with non-negative Mode C \cr
}
Default is unconstrained update for all modes, i.e., \code{const=c(0,0)}.

Output \code{cflag} gives convergence information: \code{cflag=0} if ALS algorithm converged normally, \code{cflag=1} if maximum iteration limit was reached before convergence, and \code{cflag=2} if ALS algorithm terminated abnormally due to problem with non-negativity constraints.
}
\section{Warnings }{
The ALS algorithm can perform poorly if the number of factors \code{nfac} is set too large.

Default is unconstrained ALS update, which may produce negative (invalid) Mode C weights. Use \code{const=c(0,2)} to force non-negativity via \code{\link{fnnls}}.
}

\examples{
##########   array example   ##########

# create random data array with INDSCAL structure
set.seed(3)
mydim <- c(50,5,10)
nf <- 2
X <- array(0,dim=c(rep(mydim[2],2),mydim[3]))
for(k in 1:mydim[3]) {
  X[,,k] <- as.matrix(dist(t(matrix(rnorm(prod(mydim[1:2])),mydim[1],mydim[2]))))
}

# fit INDSCAL model (unconstrained)
imod <- indscal(X,nfac=nf,nstart=1)
imod

# check solution
Xhat <- fitted(imod)
sum((array(apply(X,3,ed2sp),dim=dim(X))-Xhat)^2)
imod$SSE

# reorder and resign factors
imod$B[1:4,]
imod <- reorder(imod, 2:1)
imod$B[1:4,]
imod <- resign(imod, newsign=c(1,-1))
imod$B[1:4,]
sum((array(apply(X,3,ed2sp),dim=dim(X))-Xhat)^2)
imod$SSE

# rescale factors
colSums(imod$B^2)
colSums(imod$C^2)
imod <- rescale(imod, mode="C")
colSums(imod$B^2)
colSums(imod$C^2)
sum((array(apply(X,3,ed2sp),dim=dim(X))-Xhat)^2)
imod$SSE


##########   list example   ##########

# create random data array with INDSCAL structure
set.seed(4)
mydim <- c(100,8,20)
nf <- 3
X <- vector("list",mydim[3])
for(k in 1:mydim[3]) {
  X[[k]] <- dist(t(matrix(rnorm(prod(mydim[1:2])),mydim[1],mydim[2])))
}

# fit INDSCAL model (orthogonal B, non-negative C)
imod <- indscal(X,nfac=nf,nstart=1,const=c(1,2))
imod

# check solution
Xhat <- fitted(imod)
sum((array(unlist(lapply(X,ed2sp)),dim=mydim[c(2,2,3)])-Xhat)^2)
imod$SSE
crossprod(imod$B)


\dontrun{

##########   parallel computation   ##########

# create random data array with INDSCAL structure
set.seed(3)
mydim <- c(50,5,10)
nf <- 2
X <- array(0,dim=c(rep(mydim[2],2),mydim[3]))
for(k in 1:mydim[3]) {
  X[,,k] <- as.matrix(dist(t(matrix(rnorm(prod(mydim[1:2])),mydim[1],mydim[2]))))
}

# fit INDSCAL model (10 random starts -- sequential computation)
set.seed(1)
system.time({imod <- indscal(X,nfac=nf)})
imod

# fit INDSCAL model (10 random starts -- parallel computation)
set.seed(1)
cl <- makeCluster(detectCores())
ce <- clusterEvalQ(cl,library(multiway))
system.time({imod <- indscal(X,nfac=nf,parallel=TRUE,cl=cl)})
imod
stopCluster(cl)
}

}
