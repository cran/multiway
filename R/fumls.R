fumls <- 
  function(x, y, mode = NULL, df = NULL, knots = NULL, degree = 3,
           Boundary.knots = range(x), lower = -Inf, upper = Inf){
    # Fast Unimodal Least Squares via I-Splines
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
    iknots <- sort(c(ibasis$knots, ibasis$Boundary.knots))
    nk <- length(iknots)
    X <- cbind(1, ibasis$X)
    ncoef <- ncol(X)
    xtx <- crossprod(X)
    xty <- crossprod(X, y)
    
    # make constraint matrix
    Gmat <- cbind(0, diag(ncoef - 1L))
    bvec <- rep(0, ncoef - 1L)
    
    # check if mode is provided
    if(is.null(mode)){
      
      betas <- matrix(0, ncoef, nk)
      fits <- matrix(0, nx, nk)
      mse <- rep(0, nk)
      for(modeIndex in 1:nk){
        
        # make resigning vector
        signvec <- c(1,rep(c(1,-1), c(modeIndex, ncoef-modeIndex-1)))
        Sxtx <- diag(signvec) %*% xtx %*% diag(signvec)
        Sxty <- xty*signvec
        
        # check for lower
        Gtemp <- Gmat
        btemp <- bvec
        if(!is.infinite(lower)){
          Gtemp <- rbind(Gtemp, X[which.min(x),] * signvec, X[which.max(x),] * signvec)
          btemp <- c(btemp, rep(lower,2))
        }
        
        #browser(1>0)
        
        # solve QP problem
        qpfit <- solve.QP(Dmat = Sxtx, dvec = Sxty, Amat = t(Gtemp), bvec = btemp)
        betas[,modeIndex] <- qpfit$solution
        fits[,modeIndex] <- X %*% (betas[,modeIndex] * signvec)
        mse[modeIndex] <- mean((y - fits[,modeIndex])^2)
        
      }
      
      modeIndex <- which.min(mse)
      beta <- betas[,modeIndex]
      edf <- sum(abs(beta) > .Machine$double.eps)
      fit <- fits[,modeIndex]
      
    } else {
      
      # make resigning vector
      xmode <- IsplineBasis(mode, knots = ibasis$knots, degree = ibasis$degree, Boundary.knots = ibasis$Boundary.knots)$X
      modeIndex <- sum(xmode >= 0.5)
      signvec <- c(1,rep(c(1,-1), c(modeIndex, ncoef-modeIndex-1)))
      Sxtx <- diag(signvec) %*% xtx %*% diag(signvec)
      Sxty <- xty*signvec
      
      # check for lower
      if(!is.infinite(lower)){
        Gmat <- rbind(Gmat, X[which.min(x),] * signvec, X[which.max(x),] * signvec)
        bvec <- c(bvec, rep(lower,2))
      }
      
      # check for upper
      if(!is.infinite(upper)){
        Gmat <- rbind(Gmat, -cbind(1,xmode)*signvec)
        bvec <- c(bvec, -upper)
      }
      
      # solve QP problem
      qpfit <- solve.QP(Dmat = Sxtx, dvec = Sxty, Amat = t(Gmat), bvec = bvec)
      beta <- qpfit$solution
      edf <- sum(abs(beta) > .Machine$double.eps)
      fit <- X %*% (beta * signvec)
      
    } # end if(is.null(mode))
    
    # output results
    list(fitted.values=as.numeric(fit), coef=as.numeric(beta),
         edf=edf, knots=ibasis$knots, Boundary.knots=ibasis$Boundary.knots, 
         degree=degree, lower=lower, upper=upper, index=modeIndex)
    
  }