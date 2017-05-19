fmnls <- 
  function(x, y, df = NULL, knots = NULL, degree = 3,
           Boundary.knots = range(x), lower = -Inf, upper = Inf){
    # Fast Monotonic Least Squares via I-Splines
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 11, 2017
    
    
    ###*###   INITIAL CHECKS   ###*###
    
    # check x and y
    nx <- length(x)
    if(length(y) != nx) stop("lengths of 'x' and 'y' must match")
    x <- as.numeric(x)
    y <- as.numeric(y)
    
    # check lower and upper
    lower <- lower[1]
    upper <- upper[1]
    if(!is.infinite(lower) && !is.infinite(upper)){
      if(lower >= upper) stop("Inputs 'lower' and 'upper' must satisfy: lower < upper")
    }
    
    
    ###*###   ESTIMATE COEFFICIENTS   ###*###
    
    # make design matrix
    if(!is.null(df)) df <- df - 1L
    ibasis <- IsplineBasis(x, df = df, knots = knots, degree = degree, Boundary.knots = Boundary.knots)
    X <- cbind(1, ibasis$X)
    ncoef <- ncol(X)
    xtx <- crossprod(X)
    xty <- crossprod(X, y)
    
    # make constraint matrix
    Gmat <- cbind(0, diag(ncoef - 1L))
    bvec <- rep(0, ncoef - 1L)
    
    # check for lower
    if(!is.infinite(lower)){
      Gmat <- rbind(Gmat, X[which.min(x),])
      bvec <- c(bvec, lower)
    }
    
    # check for upper
    if(!is.infinite(upper)){
      Gmat <- rbind(Gmat, -X[which.max(x),])
      bvec <- c(bvec, -upper)
    }
    
    # call solve.QP to get coefficients
    qpfit <- solve.QP(Dmat = xtx, dvec = xty, Amat = t(Gmat), bvec = bvec)
    beta <- qpfit$solution
    edf <- sum(abs(beta) > .Machine$double.eps)
    
    # output results
    list(fitted.values=as.numeric(X %*% beta), coef=as.numeric(beta),
         edf=edf, knots=ibasis$knots, Boundary.knots=ibasis$Boundary.knots, 
         degree=degree, lower=lower, upper=upper)
    
  }