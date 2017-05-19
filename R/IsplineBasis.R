IsplineBasis <-
  function(x, df = NULL, knots = NULL, degree = 3, Boundary.knots = range(x)){
    
    # get needed M-spline basis
    if(!is.null(df)) df <- df + 1L
    Xm <- MsplineBasis(x = x, df = df, knots = knots, degree = degree,
                       intercept = TRUE, Boundary.knots = Boundary.knots)
    
    # form knots and get info
    Aknots <- sort(c(rep(Xm$Boundary.knots, degree+1L), Xm$knots))
    nx <- length(x)
    nk <- length(Aknots)
    df <- ncol(Xm$X)
    
    # make design matrix
    X <- matrix(0, nrow=nx, ncol=df-1L)
    for(h in 1:nx){
      j <- sum(x[h] >= Aknots)
      for(i in 2:df){
        if(i > j){
          X[h,i-1L] <- 0
        } else if(i < (j - degree + 1)){
          X[h,i-1L] <- 1 
        } else {
          for(m in i:j){
            X[h,i-1L] <- X[h,i-1L] + (Aknots[m+degree+1] - Aknots[m]) * Xm$X[h,m] / (degree + 1)
          }
        }
      } # end for(i in 2:df)
    } # end for(h in 1:nx)
    
    return(list(X=X, knots=Xm$knots, degree=degree, Boundary.knots=Xm$Boundary.knots))
    
  }