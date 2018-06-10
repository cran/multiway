fitted.mcr <- 
  function(object, type = c("X", "Y"), ...){
    # Calculates Fitted Values for fit MCR Models
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: February 2, 2018
    
    # check 'type'
    type <- as.character(type[1])
    if(!any(type == c("X", "Y"))) stop("Input 'type' must be either 'X' or 'Y'.")
    
    # get fitted values (type == Y)
    if(type == "Y") return(object$C %*% t(object$D))
    
    # get fitted values (type == X)
    if(object$model == "parafac"){
      xdim <- c(nrow(object$A), nrow(object$B), nrow(object$C))
      Xhat <- array(object$A %*% t(krprod(object$C, object$B)), dim = xdim)
      return(Xhat)
    } else if(object$model == "parafac2"){
      nx <- sapply(object$A, nrow)
      nfac <- ncol(object$B)
      if(min(nx) == max(nx)){
        xdim <- c(min(nx), nrow(object$B), nrow(object$C))
        Xhat <- array(0, dim = xdim)
        for(k in 1:length(object$A)) Xhat[,,k] <- object$A[[k]] %*% (diag(nfac) * object$C[k,]) %*% t(object$B)
        return(Xhat)
      } else {
        Xhat <- vector("list", length(object$A))
        for(k in 1:length(object$A)) Xhat[[k]] <- object$A[[k]] %*% (diag(nfac) * object$C[k,]) %*% t(object$B)
      } # end if(min(nx) == max(nx))
      return(Xhat)
    } else {
      xdim <- c(nrow(object$A), nrow(object$B), nrow(object$C))
      Ga <- matrix(object$G, nrow = ncol(object$A))
      Xhat <- array(object$A %*% Ga %*% t(kronecker(object$C, object$B)), dim = xdim)
      return(Xhat)
    } # end if(object$model == "parafac")
    
  }