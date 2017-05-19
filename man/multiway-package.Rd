\name{multiway-package}
\alias{multiway-package}
\alias{multiway}
\docType{package}
\title{Component Models for Multi-Way Data}
\description{
Fits multi-way component models via alternating least squares algorithms with optional constraints: orthogonal, non-negative, unimodal, monotonic, periodic, smooth, or structure. Fit models include Individual Differences Scaling, Parallel Factor Analysis (1 and 2), Simultaneous Component Analysis, and Tucker Factor Analysis.
}
\details{
%\tabular{ll}{
%Package: \tab multiway\cr
%Type: \tab Package\cr
%Version: \tab 1.0\cr
%Date: \tab 2014-01-31\cr
%License: \tab GPL (>=2)\cr
%}
\code{\link{indscal}} fits the Individual Differences Scaling model. \code{\link{parafac}} fits the 3-way and 4-way Parallel Factor Analysis-1 model. \code{\link{parafac2}} fits the 3-way and 4-way Parallel Factor Analysis-2 model. \code{\link{sca}} fits the four different Simultaneous Component Analysis models. \code{\link{tucker}} fits the 3-way and 4-way Tucker Factor Analysis model. 
}
\author{
Nathaniel E. Helwig <helwig@umn.edu>

Maintainer: Nathaniel E. Helwig <helwig@umn.edu>
}
\references{
Bro, R., & De Jong, S. (1997). A fast non-negativity-constrained least squares algorithm. \emph{Journal of Chemometrics, 11}, 393-401.

Bro, R., & Kiers, H.A.L. (2003). A new efficient method for determining the number of components in PARAFAC models. \emph{Journal of Chemometrics, 17}, 274-286.

Carroll, J. D., & Chang, J-J. (1970). Analysis of individual differences in multidimensional scaling via an n-way generalization of "Eckart-Young" decomposition. \emph{Psychometrika, 35}, 283-319.

Harshman, R. A. (1970). Foundations of the PARAFAC procedure: Models and conditions for an "explanatory" multimodal factor analysis. \emph{UCLA Working Papers in Phonetics, 16}, 1-84.

Harshman, R. A. (1972). PARAFAC2: Mathematical and technical notes. \emph{UCLA Working Papers in Phonetics, 22}, 30-44.

Harshman, R. A., & Lundy, M. E. (1994). PARAFAC: Parallel factor analysis. \emph{Computational Statistics and Data Analysis, 18}, 39-72.

Haughwout, S. P., LaVallee, R. A., & Castle, I-J. P. (2015). 
   Surveillance Report #102: Apparent Per Capita Alcohol Consumption: 
   National, State, and Regional Trends, 1977-2013. Bethesda, MD: NIAAA,
   Alcohol Epidemiologic Data System.

Helwig, N. E. (2013). The special sign indeterminacy of the direct-fitting Parafac2 model: Some implications, cautions, and recommendations, for Simultaneous Component Analysis. \emph{Psychometrika, 78}, 725-739.

Helwig, N. E. (2016). Estimating latent trends in multivariate longitudinal data via Parafac2 with functional and structural constraints. \emph{Biometrical Journal}.

Kiers, H. A. L., ten Berge, J. M. F., & Bro, R. (1999). PARAFAC2-part I: A direct-fitting algorithm for the PARAFAC2 model. \emph{Journal of Chemometrics, 13}, 275-294.

Kroonenberg, P. M., & de Leeuw, J. (1980). Principal component analysis of three-mode data by means of alternating least squares algorithms. \emph{Psychometrika, 45}, 69-97.

Moore, E. H. (1920). On the reciprocal of the general algebraic matrix. Bulletin of the American Mathematical Society 26, 394-395.

Nephew, T. M., Yi, H., Williams, G. D., Stinson, F. S., & Dufour, M.C., (2004). 
   U.S. Alcohol Epidemiologic Data Reference Manual, Vol. 1, 4th ed. 
   U.S. Apparent Consumption of Alcoholic Beverages Based on State Sales, 
   Taxation, or Receipt Data. Bethesda, MD: NIAAA, Alcohol Epidemiologic Data 
   System. NIH Publication No. 04-5563.

Penrose, R. (1950). A generalized inverse for matrices. Mathematical Proceedings of the Cambridge Philosophical Society 51, 406-413.

Ramsay, J. O. (1988). Monotone regression splines in action. \emph{Statistical Science, 3}, 425-441.

Timmerman, M. E., & Kiers, H. A. L. (2003). Four simultaneous component models for the analysis of multivariate time series from more than one subject to model intraindividual and interindividual differences. \emph{Psychometrika, 68}, 105-121.

Tucker, L. R. (1966). Some mathematical notes on three-mode factor analysis. \emph{Psychometrika, 31}, 279-311.
}
\keyword{ package }
\examples{
# See examples for indscal, parafac, parafac2, sca, and tucker
}
