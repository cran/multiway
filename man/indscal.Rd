\name{indscal}
\alias{indscal}
\title{
Individual Differences Scaling
}
\description{
Fits Carroll and Chang's Individual Differences Scaling (INDSCAL) model to 3-way dissimilarity or similarity data. Parameters are estimated via alternating least squares with optional constraints.
}
\usage{
indscal(X, nfac, nstart = 10, const = NULL, control = NULL,
        type = c("dissimilarity", "similarity"),
        Bfixed = NULL, Bstart = NULL, Bstruc = NULL, Bmodes = NULL,
        Cfixed = NULL, Cstart = NULL, Cstruc = NULL, Cmodes = NULL,
        maxit = 500, ctol = 1e-4, parallel = FALSE, cl = NULL,
        output = c("best", "all"), verbose = TRUE, backfit = FALSE)
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
  Character vector of length 2 giving the constraints for modes B and C (defaults to unconstrained for B and non-negative for C). See \code{\link[CMLS]{const}} for the 24 available options. Constraints for Mode C weights are limited to one of the 8 possible non-negative options.
}
  \item{control}{
  List of parameters controlling options for smoothness constraints. This is passed to \code{\link{const.control}}, which describes the available options.
}
  \item{type}{
  Character indicating if \code{X} contains dissimilarity data (default) or similarity data.
}
  \item{Bfixed}{
  Used to fit model with fixed Mode B weights.
}
  \item{Bstart}{
  Starting Mode B weights. Default uses random weights.
}
  \item{Bstruc}{
  Structure constraints for Mode B weights. See Note.
}
  \item{Bmodes}{
  Mode ranges for Mode B weights (for unimodality constraints). See Note.
}
  \item{Cfixed}{
  Used to fit model with fixed Mode C weights.
}
  \item{Cstart}{
  Starting Mode C weights. Default uses random weights.
}
  \item{Cstruc}{
  Structure constraints for Mode C weights. See Note.
}
  \item{Cmodes}{
  Mode ranges for Mode C weights (for unimodality constraints). See Note.
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
  If \code{TRUE}, fitting progress is printed via \code{\link{txtProgressBar}}. Ignored if \code{parallel=TRUE}.
}
  \item{backfit}{
  Should backfitting algorithm be used for \code{\link[CMLS]{cmls}}? 
}
}

\details{
Given a 3-way array \code{X = array(x,dim=c(J,J,K))} with \code{X[,,k]} denoting the \code{k}-th subject's dissimilarity matrix rating \code{J} objects, the INDSCAL model can be written as 
\tabular{c}{
\code{ Z[i,j,k] = sum B[i,r]*B[j,r]*C[k,r] + E[i,j,k] }
} 
where \code{Z} is the array of scalar products obtained from \code{X}, \code{B = matrix(b,J,R)} are the object weights, \code{C = matrix(c,K,R)} are the non-negative subject weights, and \code{E = array(e,dim=c(J,J,K))} is the 3-way residual array. The summation is for \code{r = seq(1,R)}.

Weight matrices are estimated using an alternating least squares algorithm with optional constraints.
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
\item{cflag}{Convergence flag. See Note.}
\item{const}{See argument \code{const}.}
\item{control}{See argument \code{control}.}
\item{fixed}{Logical vector indicating whether 'fixed' weights were used for each mode.}
\item{struc}{Logical vector indicating whether 'struc' constraints were used for each mode.}

Otherwise returns a list of length \code{nstart} where each element is an object of class \code{"indscal"}.
}
\references{
Carroll, J. D., & Chang, J-J. (1970). Analysis of individual differences in multidimensional scaling via an n-way generalization of "Eckart-Young" decomposition. \emph{Psychometrika, 35}, 283-319. \doi{10.1007/BF02310791}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
Structure constraints should be specified with a matrix of logicals (TRUE/FALSE), such that FALSE elements indicate a weight should be constrained to be zero. Default uses unstructured weights, i.e., a matrix of all TRUE values.

When using unimodal constraints, the \code{*modes} inputs can be used to specify the mode search range for each factor. These inputs should be matrices with dimension \code{c(2,nfac)} where the first row gives the minimum mode value and the second row gives the maximum mode value (with respect to the indicies of the given corresponding matrix).

Output \code{cflag} gives convergence information: \code{cflag = 0} if algorithm converged normally, \code{cflag = 1} if maximum iteration limit was reached before convergence, and \code{cflag = 2} if algorithm terminated abnormally due to a problem with the constraints.
}
\section{Warnings }{
The algorithm can perform poorly if the number of factors \code{nfac} is set too large.
}

\seealso{
The \code{\link{fitted.indscal}} function creates the model-implied fitted values from a fit \code{"indscal"} object.

The \code{\link{resign.indscal}} function can be used to resign factors from a fit \code{"indscal"} object.

The \code{\link{rescale.indscal}} function can be used to rescale factors from a fit \code{"indscal"} object.

The \code{\link{reorder.indscal}} function can be used to reorder factors from a fit \code{"indscal"} object.

The \code{\link[CMLS]{cmls}} function (from \bold{CMLS} package) is called as a part of the alternating least squares algorithm.
}

\examples{
##########   array example   ##########

# create random data array with INDSCAL structure
set.seed(3)
mydim <- c(50,5,10)
nf <- 2
X <- array(0, dim = c(rep(mydim[2],2), mydim[3]))
for(k in 1:mydim[3]) {
  X[,,k] <- as.matrix(dist(t(matrix(rnorm(prod(mydim[1:2])), mydim[1], mydim[2]))))
}

# fit INDSCAL model
imod <- indscal(X, nfac = nf, nstart = 1)
imod

# check solution
Xhat <- fitted(imod)
sum((array(apply(X,3,ed2sp), dim = dim(X)) - Xhat)^2)
imod$SSE

# reorder and resign factors
imod$B[1:4,]
imod <- reorder(imod, 2:1)
imod$B[1:4,]
imod <- resign(imod, newsign = c(1,-1))
imod$B[1:4,]
sum((array(apply(X,3,ed2sp), dim = dim(X)) - Xhat)^2)
imod$SSE

# rescale factors
colSums(imod$B^2)
colSums(imod$C^2)
imod <- rescale(imod, mode = "C")
colSums(imod$B^2)
colSums(imod$C^2)
sum((array(apply(X,3,ed2sp), dim = dim(X)) - Xhat)^2)
imod$SSE


##########   list example   ##########

# create random data array with INDSCAL structure
set.seed(4)
mydim <- c(100, 8, 20)
nf <- 3
X <- vector("list", mydim[3])
for(k in 1:mydim[3]) {
  X[[k]] <- dist(t(matrix(rnorm(prod(mydim[1:2])), mydim[1], mydim[2])))
}

# fit INDSCAL model (orthogonal B, non-negative C)
imod <- indscal(X, nfac = nf, nstart = 1, const = c("orthog", "nonneg"))
imod

# check solution
Xhat <- fitted(imod)
sum((array(unlist(lapply(X,ed2sp)), dim = mydim[c(2,2,3)]) - Xhat)^2)
imod$SSE
crossprod(imod$B)


\dontrun{

##########   parallel computation   ##########

# create random data array with INDSCAL structure
set.seed(3)
mydim <- c(50,5,10)
nf <- 2
X <- array(0,dim=c(rep(mydim[2],2), mydim[3]))
for(k in 1:mydim[3]) {
  X[,,k] <- as.matrix(dist(t(matrix(rnorm(prod(mydim[1:2])), mydim[1], mydim[2]))))
}

# fit INDSCAL model (10 random starts -- sequential computation)
set.seed(1)
system.time({imod <- indscal(X, nfac = nf)})
imod

# fit INDSCAL model (10 random starts -- parallel computation)
cl <- makeCluster(detectCores())
ce <- clusterEvalQ(cl,library(multiway))
clusterSetRNGStream(cl, 1)
system.time({imod <- indscal(X, nfac = nf, parallel = TRUE, cl = cl)})
imod
stopCluster(cl)
}

}

\keyword{ models }% use one of  RShowDoc("KEYWORDS")
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")
