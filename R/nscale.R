nscale <- 
  function(X, mode = 1, ssnew = NULL, newscale = 1){
    # Scale n-th Dimension of Array
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 25, 2018
    
    mode <- as.integer(mode[1])
    if(mode < 0L)(stop("mode must be nonnegative integer."))
    
    # old method: sum-of-squares
    if(!is.null(ssnew)){
      #warning("You may want to consider using the 'newscale' input instead of 'ssnew' input.")
      
      ssnew <- as.numeric(ssnew[1])
      if(ssnew<=0L)(stop("ssnew must be positive value."))
      if(mode==0L){
        
        if(is.list(X)){
          X <- mapply("*",X,sqrt(ssnew/sumsq(X)),SIMPLIFY=FALSE)
        } else{
          X <- X*sqrt(ssnew/sumsq(X))
        }
        
      } else {
        
        if(is.array(X)){
          xdim <- dim(X)
          if(mode>length(xdim)){stop("mode must be less than maximum mode.")}
          if(length(xdim)==2L){
            if(mode==1L){
              xsqs <- rowSums(X^2)
              X <- (diag(xdim[1])*sqrt(ssnew/xsqs))%*%X
            } else {
              xsqs <- colSums(X^2)
              X <- X%*%(diag(xdim[2])*sqrt(ssnew/xsqs))
            }
          } else if(length(xdim)==3L){
            if(mode==1L){
              X <- matrix(X,xdim[1],xdim[2]*xdim[3])
              X <- nscale(X,ssnew=ssnew)
              X <- array(X,dim=xdim)
            } else if(mode==2L){
              X <- matrix(aperm(X,perm=c(2,1,3)),xdim[2],xdim[1]*xdim[3])
              X <- nscale(X,ssnew=ssnew)
              X <- aperm(array(X,dim=xdim[c(2,1,3)]),perm=c(2,1,3))
            } else {
              X <- matrix(aperm(X,perm=c(3,1,2)),xdim[3],xdim[1]*xdim[2])
              X <- nscale(X,ssnew=ssnew)
              X <- aperm(array(X,dim=xdim[c(3,1,2)]),perm=c(2,3,1))
            }
          } else if(length(xdim)==4L){
            if(mode==1L){
              X <- matrix(X,xdim[1],prod(xdim[2:4]))
              X <- nscale(X,ssnew=ssnew)
              X <- array(X,dim=xdim)
            } else if(mode==2L){
              X <- matrix(aperm(X,perm=c(2,1,3,4)),xdim[2],prod(xdim[c(1,3,4)]))
              X <- nscale(X,ssnew=ssnew)
              X <- aperm(array(X,dim=xdim[c(2,1,3,4)]),perm=c(2,1,3,4))
            } else if(mode==3L){
              X <- matrix(aperm(X,perm=c(3,1,2,4)),xdim[3],prod(xdim[c(1,2,4)]))
              X <- nscale(X,ssnew=ssnew)
              X <- aperm(array(X,dim=xdim[c(3,1,2,4)]),perm=c(2,3,1,4))
            } else {
              X <- matrix(aperm(X,perm=c(4,1,2,3)),xdim[4],prod(xdim[1:3]))
              X <- nscale(X,ssnew=ssnew)
              X <- aperm(array(X,dim=xdim[c(4,1,2,3)]),perm=c(2,3,4,1))
            }
          } else{stop("nscale only rescales 2-, 3-, and 4-way arrays.")}
          
        } else if(is.list(X)){
          
          lxdim <- length(dim(X[[1]]))
          if(mode>(lxdim+1L)){stop("mode must be less than maximum mode.")}
          dimchk <- sapply(X, function(x) length(dim(x))==lxdim)
          if(any(!dimchk)) stop("When 'X' is a list, each element of X \nmust be an array with the same number of dimensions.")
          
          if(lxdim==2L){
            if(mode==1L){
              rowSS <- sapply(X, function(x) rowSums(x^2))
              if(is.list(rowSS)) stop("Incompatible dimensions for scaling mode across lists")
              rowSS <- rowSums(rowSS)
              rowMat <- diag( sqrt(ssnew/rowSS) )
              for(k in 1:length(X)) X[[k]] <- rowMat %*% X[[k]]
            } else if(mode==2L){
              colSS <- sapply(X, function(x) colSums(x^2))
              if(is.list(colSS)) stop("Incompatible dimensions for scaling mode across lists")
              colSS <- rowSums(colSS)
              colMat <- diag( sqrt(ssnew/colSS) )
              for(k in 1:length(X)) X[[k]] <- X[[k]] %*% colMat
            } else {
              X <- mapply(nscale,X,MoreArgs=list(mode=0,ssnew=ssnew),SIMPLIFY=FALSE)
            }
          } else if(lxdim==3L){
            if(mode==1L){
              rowSS <- sapply(X, function(x) apply(x,1,sumsq))
              if(is.list(rowSS)) stop("Incompatible dimensions for scaling mode across lists")
              rowSS <- rowSums(rowSS)
              rowMat <- diag( sqrt(ssnew/rowSS) )
              for(k in 1:length(X)){
                xdim <- dim(X[[k]])
                Xmat <- matrix(X[[k]],nrow=nrow(rowMat))
                X[[k]] <- array(rowMat %*% Xmat, dim=xdim)
              }
            } else if(mode==2L){
              colSS <- sapply(X, function(x) apply(x,2,sumsq))
              if(is.list(colSS)) stop("Incompatible dimensions for scaling mode across lists")
              colSS <- rowSums(colSS)
              colMat <- diag( sqrt(ssnew/colSS) )
              for(k in 1:length(X)){
                xdim <- dim(X[[k]])
                Xmat <- matrix(aperm(X[[k]],c(2,1,3)),nrow=nrow(colMat))
                X[[k]] <- aperm(array(colMat %*% Xmat, dim=xdim[c(2,1,3)]),c(2,1,3))
              }
            } else if(mode==3L){
              slabSS <- sapply(X, function(x) apply(x,3,sumsq))
              if(is.list(slabSS)) stop("Incompatible dimensions for scaling mode across lists")
              slabSS <- rowSums(slabSS)
              slabMat <- diag( sqrt(ssnew/slabSS) )
              for(k in 1:length(X)){
                xdim <- dim(X[[k]])
                Xmat <- matrix(aperm(X[[k]],c(3,1,2)),nrow=nrow(slabMat))
                X[[k]] <- aperm(array(slabMat %*% Xmat, dim=xdim[c(3,1,2)]),c(2,3,1))
              }
            } else{
              X <- mapply(nscale,X,MoreArgs=list(mode=0,ssnew=ssnew),SIMPLIFY=FALSE)
            }
          } else {
            stop("When 'X' is a list, each element of X \nmust be a matrix or 3-way array.")
          }
        } else{stop("Input X must be array or list with array elements.")}
        
      } # end if(mode==0L)
      
      return(X)
      
    } # end if(!is.null(ssnew))
    
    # new method: root-mean-square
    newscale <- as.numeric(newscale[1])
    if(newscale <= 0L)(stop("newscale must be positive value."))
    if(mode == 0L){
      
      if(is.list(X)){
        X <- mapply("*", X, newscale / sqrt(meansq(X)), SIMPLIFY = FALSE)
      } else{
        X <- X * newscale / sqrt(meansq(X))
      }
      
    } else {
      
      if(is.array(X)){
        xdim <- dim(X)
        if(mode > length(xdim)){stop("mode must be less than maximum mode.")}
        if(length(xdim) == 2L){
          if(mode == 1L){
            xsqs <- rowMeans(X^2)
            X <- X * matrix(newscale / sqrt(xsqs), nrow = xdim[1], ncol = xdim[2])
          } else {
            xsqs <- colMeans(X^2)
            X <- X * matrix(newscale / sqrt(xsqs), nrow = xdim[1], ncol = xdim[2], byrow = TRUE)
          }
        } else if(length(xdim) == 3L){
          if(mode == 1L){
            X <- matrix(X, nrow = xdim[1], ncol = xdim[2]*xdim[3])
            X <- nscale(X, newscale = newscale)
            X <- array(X, dim = xdim)
          } else if(mode == 2L){
            X <- matrix(aperm(X, perm = c(2,1,3)), nrow = xdim[2], ncol = xdim[1]*xdim[3])
            X <- nscale(X, newscale = newscale)
            X <- aperm(array(X, dim = xdim[c(2,1,3)]), perm = c(2,1,3))
          } else {
            X <- matrix(aperm(X, perm = c(3,1,2)), nrow = xdim[3], ncol = xdim[1]*xdim[2])
            X <- nscale(X, newscale = newscale)
            X <- aperm(array(X, dim = xdim[c(3,1,2)]), perm = c(2,3,1))
          }
        } else if(length(xdim) == 4L){
          if(mode == 1L){
            X <- matrix(X, nrow = xdim[1], ncol = prod(xdim[2:4]))
            X <- nscale(X, newscale = newscale)
            X <- array(X, dim = xdim)
          } else if(mode == 2L){
            X <- matrix(aperm(X, perm = c(2,1,3,4)), nrow = xdim[2], ncol = prod(xdim[c(1,3,4)]))
            X <- nscale(X, newscale = newscale)
            X <- aperm(array(X, dim = xdim[c(2,1,3,4)]), perm = c(2,1,3,4))
          } else if(mode == 3L){
            X <- matrix(aperm(X, perm = c(3,1,2,4)), nrow = xdim[3], ncol = prod(xdim[c(1,2,4)]))
            X <- nscale(X, newscale = newscale)
            X <- aperm(array(X, dim = xdim[c(3,1,2,4)]), perm = c(2,3,1,4))
          } else {
            X <- matrix(aperm(X, perm = c(4,1,2,3)), nrow = xdim[4], ncol = prod(xdim[1:3]))
            X <- nscale(X, newscale = newscale)
            X <- aperm(array(X, dim = xdim[c(4,1,2,3)]), perm = c(2,3,4,1))
          }
        } else{stop("nscale only scales 2-, 3-, and 4-way arrays.")}
        
      } else if(is.list(X)){
        
        lengthX <- length(X)
        xdim <- dim(X[[1]])
        lxdim <- length(xdim)
        if(mode > (lxdim + 1L)){stop("mode must be less than maximum mode.")}
        dimchk <- sapply(X, function(x) length(dim(x)) == lxdim)
        if(any(!dimchk)) stop("When 'X' is a list, each element of X \nmust be an array with the same number of dimensions.")
        
        if(lxdim == 2L){
          if(mode == 1L){
            rowSS <- sapply(X, function(x) rowSums(x^2))
            rowNS <- sapply(X, ncol)
            if(is.list(rowSS)) stop("Incompatible dimensions for scaling mode across lists")
            rowSS <- rowSums(rowSS) / sum(rowNS)
            theta <- newscale / sqrt(rowSS)
            for(k in 1:length(X)) X[[k]] <- X[[k]] * matrix(theta, nrow = xdim[1], ncol = rowNS[k])
          } else if(mode == 2L){
            colSS <- sapply(X, function(x) colSums(x^2))
            colNS <- sapply(X, nrow)
            if(is.list(colSS)) stop("Incompatible dimensions for scaling mode across lists")
            colSS <- rowSums(colSS) / sum(colNS)
            theta <- newscale / sqrt(colSS)
            for(k in 1:length(X)) X[[k]] <- X[[k]] * matrix(theta, nrow = colNS[k], ncol = xdim[2])
          } else {
            X <- mapply(nscale, X, MoreArgs = list(mode = 0, newscale = newscale), SIMPLIFY = FALSE)
          }
        } else if(lxdim == 3L){
          if(mode == 1L){
            rowSS <- sapply(X, function(x) apply(x, 1, sumsq))
            rowNS <- sapply(X, function(x) prod(dim(x)[-1]))
            if(is.list(rowSS)) stop("Incompatible dimensions for scaling mode across lists")
            rowSS <- rowSums(rowSS) / sum(rowNS)
            rowMat <- diag( newscale / sqrt(rowSS) )
            for(k in 1:length(X)){
              xdim <- dim(X[[k]])
              Xmat <- matrix(X[[k]], nrow = nrow(rowMat))
              X[[k]] <- array(rowMat %*% Xmat, dim = xdim)
            }
          } else if(mode == 2L){
            colSS <- sapply(X, function(x) apply(x, 2, sumsq))
            colNS <- sapply(X, function(x) prod(dim(x)[-2]))
            if(is.list(colSS)) stop("Incompatible dimensions for scaling mode across lists")
            colSS <- rowSums(colSS) / sum(colNS)
            colMat <- diag( newscale / sqrt(colSS) )
            for(k in 1:length(X)){
              xdim <- dim(X[[k]])
              Xmat <- matrix(aperm(X[[k]], perm = c(2,1,3)), nrow = nrow(colMat))
              X[[k]] <- aperm(array(colMat %*% Xmat, dim = xdim[c(2,1,3)]), perm = c(2,1,3))
            }
          } else if(mode == 3L){
            slabSS <- sapply(X, function(x) apply(x, 3, sumsq))
            slabNS <- sapply(X, function(x) prod(dim(x)[-3]))
            if(is.list(slabSS)) stop("Incompatible dimensions for scaling mode across lists")
            slabSS <- rowSums(slabSS) / sum(slabNS)
            slabMat <- diag( newscale / sqrt(slabSS) )
            for(k in 1:length(X)){
              xdim <- dim(X[[k]])
              Xmat <- matrix(aperm(X[[k]], perm = c(3,1,2)), nrow = nrow(slabMat))
              X[[k]] <- aperm(array(slabMat %*% Xmat, dim = xdim[c(3,1,2)]), perm = c(2,3,1))
            }
          } else{
            X <- mapply(nscale, X, MoreArgs = list(mode = 0, newscale = newscale), SIMPLIFY = FALSE)
          }
        } else {
          stop("When 'X' is a list, each element of X \nmust be a matrix or 3-way array.")
        }
      } else{stop("Input X must be array or list with array elements.")}
      
    } # end if(mode==0L)
    
    return(X)
    
  }