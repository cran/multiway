fitted.parafac2 <-
  function(object,simplify=TRUE,...){
    # Calculates Fitted Values (lists) for fit Parafac2 Models
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 26, 2018
    
    mydim <- c(NA,nrow(object$B),nrow(object$C))
    nf <- ncol(object$B)
    nk <- sapply(object$A,nrow)
    if(is.null(object$D)){
      fit <- vector("list",mydim[3])
      for(k in 1:mydim[3]){
        fit[[k]] <- tcrossprod(object$A[[k]]%*%(diag(nf)*object$C[k,]),object$B)
      }
    } else {
      mydim <- c(mydim,nrow(object$D))
      fit <- vector("list",mydim[4])
      CkrB <- krprod(object$C, object$B)
      for(k in 1:mydim[4]){
        fit[[k]] <- array(tcrossprod(object$A[[k]]%*%(diag(nf)*object$D[k,]),CkrB),
                          dim=c(nk[k],mydim[2],mydim[3]))
      }
    }
    
    if(min(nk) == max(nk) && simplify) {
      mydim[1] <- nk[1]
      fit <- array(unlist(fit), dim=mydim)
    }
    fit
    
  }