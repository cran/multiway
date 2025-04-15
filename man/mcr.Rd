\name{mcr}
\alias{mcr}
\title{
Multiway Covariates Regression
}
\description{
Fits Smilde and Kiers's Multiway Covariates Regression (MCR) model to connect a 3-way predictor array and a 2-way response array that share a common mode. Parameters are estimated via alternating least squares with optional constraints.
}
\usage{
mcr(X, Y, nfac = 1, alpha = 0.5, nstart = 10, 
    model = c("parafac", "parafac2", "tucker"),
    const = NULL, control = NULL, weights = NULL,
    Afixed = NULL, Bfixed = NULL, Cfixed = NULL, Dfixed = NULL,
    Astart = NULL, Bstart = NULL, Cstart = NULL, Dstart = NULL,
    Astruc = NULL, Bstruc = NULL, Cstruc = NULL, Dstruc = NULL,
    Amodes = NULL, Bmodes = NULL, Cmodes = NULL, Dmodes = NULL,
    maxit = 500, ctol = 1e-4, parallel = FALSE, cl = NULL, 
    output = c("best", "all"), verbose = TRUE, backfit = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  Three-way predictor array with \code{dim = c(I,J,K)}.
}
  \item{Y}{
  Two-way response array with \code{dim = c(K,L)}.
}
  \item{nfac}{
  Number of factors.
}
  \item{alpha}{
  Tuning parameter between 0 and 1. 
}
  \item{nstart}{
  Number of random starts.
}
  \item{model}{
  Model for \code{X}. Defaults to "parafac".
}
  \item{const}{
  Character vector of length 4 giving the constraints for \code{A}, \code{B}, \code{C}, and \code{D} (defaults to unconstrained). See \code{\link[CMLS]{const}} for the 24 available options. Ignored if \code{model = "tucker"}.
}
  \item{control}{
  List of parameters controlling options for smoothness constraints. This is passed to \code{\link{const.control}}, which describes the available options.
}
  \item{weights}{
  Vector of length \code{K} giving non-negative weights for fitting via weighted least squares. Defaults to vector of ones.
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
  Structure constraints for Mode C weights. Ignored.
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
  Mode ranges for Mode C weights (for unimodality constraints). Ignored.
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
Given a predictor array \code{X = array(x, dim=c(I,J,K))} and a response matrix \code{Y = matrix(y, nrow=K, ncol=L)}, the multiway covariates regression (MCR) model assumes a tensor model for \code{X} and a bilinear model for \code{Y}, which are linked through a common \code{C} weight matrix. For example, using the Parafac model for \code{X}, the MCR model has the form
\tabular{l}{
\code{ X[i,j,k] = sum A[i,r]*B[j,r]*C[k,r] + Ex[i,j,k] } \cr
and \cr
\code{   Y[k,l] = sum C[k,r]*D[l,r] + Ey[k,l] } \cr
}

Parameter matrices are estimated by minimizing the loss function
\tabular{c}{
\code{LOSS = alpha * (SSE.X / SSX) + (1 - alpha) * (SSE.Y / SSY)}
}
where 
\tabular{l}{
\code{SSE.X = sum((X - Xhat)^2)} \cr
\code{SSE.Y = sum((Y - Yhat)^2)} \cr
  \code{SSX = sum(X^2)} \cr
  \code{SSY = sum(Y^2)} \cr
}

When \code{weights} are input, \code{SSE.X}, \code{SSE.Y}, \code{SSX}, and \code{SSY} are replaced by the corresponding weighted versions.
}
\value{
\item{A}{Predictor A weight matrix.}
\item{B}{Predictor B weight matrix.}
\item{C}{Common C weight matrix.}
\item{D}{Response D weight matrix.}
\item{W}{Coefficients. See Note.}
\item{LOSS}{Value of \code{LOSS} function.}
\item{SSE}{Sum of Squared Errors for \code{X} and \code{Y}.}
\item{Rsq}{R-squared value for \code{X} and \code{Y}.}
\item{iter}{Number of iterations.}
\item{cflag}{Convergence flag. See Note.}
\item{model}{See argument \code{model}.}
\item{const}{See argument \code{const}.}
\item{control}{See argument \code{control}.}
\item{weights}{See argument \code{weights}.}
\item{alpha}{See argument \code{alpha}.}
\item{fixed}{Logical vector indicating whether 'fixed' weights were used for each matrix.}
\item{struc}{Logical vector indicating whether 'struc' constraints were used for each matrix.}
\item{Phi}{Mode A crossproduct matrix. Only if \code{model = "parafac2"}.}
\item{G}{Core array. Only if \code{model = "tucker"}.}
}
\references{
Smilde, A. K., & Kiers, H. A. L. (1999). Multiway covariates regression models, \emph{Journal of Chemometrics, 13}, 31-48. \doi{10.1002/(SICI)1099-128X(199901/02)13:1<31::AID-CEM528>3.0.CO;2-P}
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>
}
\note{
When \code{model = "parafac2"}, the arguments \code{Afixed}, \code{Astart}, and \code{Astruc} are treated as the arguments \code{Gfixed}, \code{Gstart}, and \code{Gstruc} from the \code{\link{parafac2}} function.

Structure constraints should be specified with a matrix of logicals (TRUE/FALSE), such that FALSE elements indicate a weight should be constrained to be zero. Default uses unstructured weights, i.e., a matrix of all TRUE values. Structure constraints are ignored if \code{model = "tucker"}.

When using unimodal constraints, the \code{*modes} inputs can be used to specify the mode search range for each factor. These inputs should be matrices with dimension \code{c(2,nfac)} where the first row gives the minimum mode value and the second row gives the maximum mode value (with respect to the indicies of the corresponding weight matrix).

\code{C = Xc \%*\% W} where \code{Xc = matrix(aperm(X,c(3,1,2)),K)}

Output \code{cflag} gives convergence information: \code{cflag = 0} if algorithm converged normally, \code{cflag = 1} if maximum iteration limit was reached before convergence, and \code{cflag = 2} if algorithm terminated abnormally due to a problem with the constraints.
}

\seealso{
The \code{\link{fitted.mcr}} function creates the model-implied fitted values from a fit \code{"mcr"} object.

The \code{\link{resign.mcr}} function can be used to resign factors from a fit \code{"mcr"} object.

The \code{\link{rescale.mcr}} function can be used to rescale factors from a fit \code{"mcr"} object.

The \code{\link{reorder.mcr}} function can be used to reorder factors from a fit \code{"mcr"} object.

The \code{\link[CMLS]{cmls}} function (from \bold{CMLS} package) is called as a part of the alternating least squares algorithm.

See \code{\link{parafac}}, \code{\link{parafac2}}, and \code{\link{tucker}} for more information about the Parafac, Parafac2, and Tucker models. 
}
\examples{

##########   multiway covariates regression   ##########

# create random data array with Parafac structure
set.seed(3)
mydim <- c(10, 20, 100)
nf <- 2
Amat <- matrix(rnorm(mydim[1]*nf), mydim[1], nf)
Bmat <- matrix(rnorm(mydim[2]*nf), mydim[2], nf)
Cmat <- matrix(rnorm(mydim[3]*nf), mydim[3], nf)
Xmat <- tcrossprod(Amat, krprod(Cmat, Bmat))
Xmat <- array(Xmat, dim = mydim)
EX <- array(rnorm(prod(mydim)), dim = mydim)
EX <- nscale(EX, 0, ssnew = sumsq(Xmat))   # SNR = 1
X <- Xmat + EX

# create response array
ydim <- c(mydim[3], 4)
Dmat <- matrix(rnorm(ydim[2]*nf), ydim[2], nf)
Ymat <- tcrossprod(Cmat, Dmat)
EY <- array(rnorm(prod(ydim)), dim = ydim)
EY <- nscale(EY, 0, ssnew = sumsq(Ymat))   # SNR = 1
Y <- Ymat + EY

# fit MCR model
mcr(X, Y, nfac = nf, nstart = 1)
mcr(X, Y, nfac = nf, nstart = 1, model = "parafac2")
mcr(X, Y, nfac = nf, nstart = 1, model = "tucker")



\dontrun{

##########   parallel computation   ##########

# create random data array with Parafac structure
set.seed(3)
mydim <- c(10, 20, 100)
nf <- 2
Amat <- matrix(rnorm(mydim[1]*nf), mydim[1], nf)
Bmat <- matrix(rnorm(mydim[2]*nf), mydim[2], nf)
Cmat <- matrix(rnorm(mydim[3]*nf), mydim[3], nf)
Xmat <- tcrossprod(Amat, krprod(Cmat, Bmat))
Xmat <- array(Xmat, dim = mydim)
EX <- array(rnorm(prod(mydim)), dim = mydim)
EX <- nscale(EX, 0, ssnew = sumsq(Xmat))   # SNR = 1
X <- Xmat + EX

# create response array
ydim <- c(mydim[3], 4)
Dmat <- matrix(rnorm(ydim[2]*nf), ydim[2], nf)
Ymat <- tcrossprod(Cmat, Dmat)
EY <- array(rnorm(prod(ydim)), dim = ydim)
EY <- nscale(EY, 0, ssnew = sumsq(Ymat))   # SNR = 1
Y <- Ymat + EY

# fit MCR-Parafac model (10 random starts -- sequential computation)
set.seed(1)
system.time({mod <- mcr(X, Y, nfac = nf)})
mod

# fit MCR-Parafac model (10 random starts -- parallel computation)
cl <- makeCluster(detectCores())
ce <- clusterEvalQ(cl, library(multiway))
clusterSetRNGStream(cl, 1)
system.time({mod <- mcr(X, Y, nfac = nf, parallel = TRUE, cl = cl)})
mod
stopCluster(cl)

}


}
\keyword{ models }
\keyword{ multivariate }
\keyword{ regression }
