fsmls <- 
  function(x, y, df = NULL, knots = NULL, degree = 3, 
           Boundary.knots = range(x), nonneg=FALSE, periodic=FALSE){
    # Fast Smooth Least Squares via M-Splines
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 12, 2017
    
    
    ###*###   INITIAL CHECKS   ###*###
    
    # check x and y
    nx <- length(x)
    if(length(y) != nx) stop("lengths of 'x' and 'y' must match")
    x <- as.numeric(x)
    y <- as.numeric(y)
    
    
    ###*###   ESTIMATE COEFFICIENTS   ###*###
    
    # make design matrix
    if(periodic){
      mbasis <- MsplineBasis(x, df = df, knots = knots, degree = degree, 
                             intercept = FALSE, Boundary.knots = Boundary.knots)
      mbasis$X <- cbind(1, mbasis$X[,-ncol(mbasis$X)])
    } else {
      mbasis <- MsplineBasis(x, df = df, knots = knots, degree = degree, 
                             intercept = TRUE, Boundary.knots = Boundary.knots)
    }
    ncoef <- ncol(mbasis$X)
    xtx <- crossprod(mbasis$X)
    xty <- crossprod(mbasis$X, y)
    
    # solve for coefficients
    if(nonneg){
      Gmat <- diag(ncoef)
      bvec <- rep(0, ncoef)
      qpfit <- solve.QP(Dmat = xtx, dvec = xty, Amat = t(Gmat), bvec = bvec)
      beta <- qpfit$solution
    } else {
      beta <- solve(xtx) %*% xty
    }
    
    # output results
    edf <- sum(abs(beta) > .Machine$double.eps)
    list(fitted.values=as.numeric(mbasis$X %*% beta), coef=as.numeric(beta),
         edf=edf, knots=mbasis$knots, Boundary.knots=mbasis$Boundary.knots, 
         degree=degree, nonneg=nonneg, periodic=periodic)
    
  }