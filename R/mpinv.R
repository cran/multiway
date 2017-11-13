mpinv <- 
  function(X, tol = NULL){
    # Moore-Penrose Pseudoinverse
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: September 5, 2017
    
    X <- as.matrix(X)
    if(is.null(tol)){
      tol <- .Machine$double.eps * max(dim(X))
    } else {
      tol <- as.numeric(tol[1])
      if(tol <= 0) stop("Input 'tol' must be greater than zero.")
    }
    xsvd <- svd(X)
    nze <- sum( xsvd$d > (tol*xsvd$d[1]) )
    if(nze > 1L){
      return( xsvd$v[,1:nze] %*% diag(1/xsvd$d[1:nze]) %*% t(xsvd$u[,1:nze]) )
    } else if (nze == 1L) {
      return( outer(xsvd$v[,1],xsvd$u[,1]) / xsvd$d[1] )
    } else {
      return( matrix(0, nrow(X), ncol(X)) )
    }
    
}