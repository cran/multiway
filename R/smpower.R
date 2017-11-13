smpower <- 
  function(X, power = 0.5, tol = NULL){
    # Symmetric Matrix Power
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: September 5, 2017
    
    X <- as.matrix(X)
    xdim <- dim(X)
    if(xdim[1] != xdim[2]) stop("Input 'X' must be a symmetric matrix.")
    if(is.null(tol)){
      tol <- .Machine$double.eps * xdim[1]
    } else {
      tol <- as.numeric(tol[1])
      if(tol <= 0) stop("Input 'tol' must be greater than zero.")
    }
    Xeig <- eigen(X, symmetric = TRUE)
    nze <- sum( Xeig$values > (tol*Xeig$values[1]) )
    if(nze > 1L){
      return( tcrossprod(Xeig$vectors[,1:nze] %*% diag(sqrt(Xeig$values[1:nze])^power)) )
    } else if(nze == 1L){
      return( outer(Xeig$vectors[,1], Xeig$vectors[,1]) * (Xeig$values[1]^power) )
    } else {
      return( matrix(0, xdim[1], xdim[2]) )
    }
    
  }