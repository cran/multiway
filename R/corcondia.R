corcondia <- 
  function(X, object, divisor = c("nfac", "core")){
    # Core Consistency Diagnostic
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: November 12, 2017
    
    clob <- class(object)
    if(clob=="parafac"){
      
      nfac <- ncol(object$A)
      Ai <- mpinv(object$A)
      imat <- kronecker(mpinv(object$C), mpinv(object$B))
      if(is.null(object$D)){
        g <- tcrossprod(Ai %*% matrix(X, nrow = nrow(object$A)), imat)
        G <- array(g, dim=rep(nfac,3))
        super <- array(0, dim=rep(nfac,3))
        for(k in 1:nfac) super[k,k,k] <- 1
      } else {
        imat <- kronecker(mpinv(object$D), imat)
        g <- tcrossprod(Ai %*% matrix(X, nrow = nrow(object$A)), imat)
        G <- array(g, dim=rep(nfac,4))
        super <- array(0, dim=rep(nfac,4))
        for(k in 1:nfac) super[k,k,k,k] <- 1
      }
      
    } else if(clob=="parafac2"){
      
      if(is.array(X)){
        xdim <- dim(X)
        lxdim <- length(xdim)
        if(lxdim==3L){
          Xlist <- vector("list",nrow(object$C))
          for(k in 1:nrow(object$C)) Xlist[[k]] <- X[,,k]
        } else {
          Xlist <- vector("list",nrow(object$D))
          for(k in 1:nrow(object$D)) Xlist[[k]] <- X[,,,k]
        }
        X <- Xlist
        rm(Xlist)
      }
      
      nfac <- ncol(object$B)
      Ai <- smpower(object$Phi, power = -1)
      imat <- kronecker(mpinv(object$C), mpinv(object$B))
      if(is.null(object$D)){
        for(k in 1:nrow(object$C)) X[[k]] <- crossprod(object$A[[k]], X[[k]])
        g <- tcrossprod(Ai %*% matrix(unlist(X), nrow = nfac), imat)
        G <- array(g, dim=rep(nfac,3))
        super <- array(0, dim=rep(nfac,3))
        for(k in 1:nfac) super[k,k,k] <- 1
      } else {
        imat <- kronecker(mpinv(object$D), imat)
        for(k in 1:nrow(object$D)) X[[k]] <- crossprod(object$A[[k]], matrix(X[[k]], nrow = nrow(object$A[[k]])))
        g <- tcrossprod(Ai %*% matrix(unlist(X), nrow = nfac), imat)
        G <- array(g, dim=rep(nfac,4))
        super <- array(0, dim=rep(nfac,4))
        for(k in 1:nfac) super[k,k,k,k] <- 1
      }
      
    } else {
      stop("Input 'object' must be object of class 'parafac' or 'parafac2'")
    }
    
    if(divisor[1]=="nfac"){
      corcon <- 100 * (1 - sum((G-super)^2)/nfac)
    } else {
      corcon <- 100 * (1 - sum((G-super)^2)/sum(G^2))
    }
    
    return(corcon)
  
}