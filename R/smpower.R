smpower <- 
  function(X,power=0.5,tol=.Machine$double.eps){
    # Symmetric Matrix Power
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 15, 2017
    
    X <- as.matrix(X)
    Xeig <- eigen(X, symmetric=TRUE)
    nze <- sum( Xeig$values > (tol*Xeig$values[1]) )
    if(nze == 1L){
      return( outer(Xeig$vectors[,1], Xeig$vectors[,1]) * (Xeig$values[1]^power) )
    } else {
      return( tcrossprod(Xeig$vectors[,1:nze] %*% diag(sqrt(Xeig$values[1:nze])^power)) )
    }
    
  }