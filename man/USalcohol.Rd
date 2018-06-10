\name{USalcohol}
\alias{USalcohol}
\docType{data}
\title{
United States Alcohol Consumption Data (1970-2013)
}
\description{
This dataset contains yearly (1970-2013) consumption data from the 50 United States and the District of Columbia for three types of alcoholic beverages: spirits, wine, and beer. The data were obtained from the National Institute on Alcohol Abuse and Alcoholism (NIAAA) Surveillance Report #102 (see below link).
}
\usage{data("USalcohol")}
\format{
  A data frame with 6732 observations on the following 8 variables.
  \describe{
    \item{\code{year}}{integer Year (1970-2013)}
    \item{\code{state}}{factor State Name (51 levels)}
    \item{\code{region}}{factor Region Name (4 levels)}
    \item{\code{type}}{factor Beverage Type (3 levels)}
    \item{\code{beverage}}{numeric Beverage Consumed (thousands of gallons)}
    \item{\code{ethanol}}{numeric Absolute Alcohol Consumed (thousands of gallons)}
    \item{\code{pop14}}{numeric Population Age 14 and Older (thousands of people)}
    \item{\code{pop21}}{numeric Population Age 21 and Older (thousands of people)}
  }
}
\details{
In the data source, the population age 21 and older for Mississippi in year 1989 is reported to be 3547.839 thousand, which is incorrect. In this dataset, the miscoded population value has been replaced with the average of the corresponding 1988 population (1709 thousand) and the 1990 population (1701.527 thousand). 
}
\source{
\url{https://pubs.niaaa.nih.gov/publications/surveillance102/pcyr19702013.txt}
}
\references{
Haughwout, S. P., LaVallee, R. A., & Castle, I-J. P. (2015). 
   Surveillance Report #102: Apparent Per Capita Alcohol Consumption: 
   National, State, and Regional Trends, 1977-2013. Bethesda, MD: NIAAA,
   Alcohol Epidemiologic Data System.
   
Helwig, N. E. (2017). Estimating latent trends in multivariate longitudinal data via Parafac2 with functional and structural constraints. \emph{Biometrical Journal, 59}(4), 783-803.

Nephew, T. M., Yi, H., Williams, G. D., Stinson, F. S., & Dufour, M.C., (2004). 
   U.S. Alcohol Epidemiologic Data Reference Manual, Vol. 1, 4th ed. 
   U.S. Apparent Consumption of Alcoholic Beverages Based on State Sales, 
   Taxation, or Receipt Data. Bethesda, MD: NIAAA, Alcohol Epidemiologic Data 
   System. NIH Publication No. 04-5563.
}
\examples{
# load data and print first six rows
data(USalcohol)
head(USalcohol)

# form tensor (time x variables x state)
Xbev <- with(USalcohol, tapply(beverage/pop21, list(year, type, state), c))
Xeth <- with(USalcohol, tapply(ethanol/pop21, list(year, type, state), c))
X <- array(0, dim=c(44, 6, 51))
X[, c(1,3,5) ,] <- Xbev
X[, c(2,4,6) ,] <- Xeth
dnames <- dimnames(Xbev)
dnames[[2]] <- c(paste0(dnames[[2]],".bev"), paste0(dnames[[2]],".eth"))[c(1,4,2,5,3,6)]
dimnames(X) <- dnames

# center each variable across time (within state)
Xc <- ncenter(X, mode = 1)

# scale each variable to have mean square of 1 (across time and states)
Xs <- nscale(Xc, mode = 2)

# fit parafac model with 3 factors
set.seed(1)
pfac <- parafac(Xs, nfac = 3, nstart = 1)

# fit parafac model with functional constraints
set.seed(1)
pfacF <- parafac(Xs, nfac = 3, nstart = 1, 
                 const = c("smooth", NA, NA))

# fit parafac model with functional and structural constraints
Bstruc <- matrix(c(rep(c(TRUE,FALSE), c(2,4)), 
                   rep(c(FALSE,TRUE,FALSE), c(2,2,2)),
                   rep(c(FALSE,TRUE), c(4,2))), nrow=6, ncol=3)
set.seed(1)
pfacFS <- parafac(Xs, nfac = 3, nstart = 1, 
                  const = c("smooth", NA, NA), Bstruc = Bstruc)

}
\keyword{datasets}
