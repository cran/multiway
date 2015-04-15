nrescale <- 
  function(X,nmode=1,ssnew=1){
    # Rescale n-th Dimension of Array
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: March 25, 2015
    
    nmode <- as.integer(nmode[1])
    if(nmode<0L)(stop("nmode must be nonnegative integer."))
    ssnew <- as.numeric(ssnew[1])
    if(ssnew<1L)(stop("ssnew must be positive value."))
    if(nmode==0L){
      
      if(is.list(X)){
        X <- mapply("*",X,sqrt(ssnew/sumsq(X)),SIMPLIFY=FALSE)
      } else{
        X <- X*sqrt(ssnew/sumsq(X))
      }
      
    } else {
      
      if(is.array(X)){
        xdim <- dim(X)
        if(nmode>length(xdim)){stop("nmode must be less than maximum mode.")}
        if(length(xdim)==2L){
          if(nmode==1L){
            xsqs <- rowSums(X^2)
            X <- (diag(xdim[1])*sqrt(ssnew/xsqs))%*%X
          } else {
            xsqs <- colSums(X^2)
            X <- X%*%(diag(xdim[2])*sqrt(ssnew/xsqs))
          }
        } else if(length(xdim)==3L){
          if(nmode==1L){
            X <- matrix(X,xdim[1],xdim[2]*xdim[3])
            X <- nrescale(X,ssnew=ssnew)
            X <- array(X,dim=xdim)
          } else if(nmode==2L){
            X <- matrix(aperm(X,perm=c(2,1,3)),xdim[2],xdim[1]*xdim[3])
            X <- nrescale(X,ssnew=ssnew)
            X <- aperm(array(X,dim=xdim[c(2,1,3)]),perm=c(2,1,3))
          } else {
            X <- matrix(aperm(X,perm=c(3,1,2)),xdim[3],xdim[1]*xdim[2])
            X <- nrescale(X,ssnew=ssnew)
            X <- aperm(array(X,dim=xdim[c(3,1,2)]),perm=c(2,3,1))
          }
        } else if(length(xdim)==4L){
          if(nmode==1L){
            X <- matrix(X,xdim[1],prod(xdim[2:4]))
            X <- nrescale(X,ssnew=ssnew)
            X <- array(X,dim=xdim)
          } else if(nmode==2L){
            X <- matrix(aperm(X,perm=c(2,1,3,4)),xdim[2],prod(xdim[c(1,3,4)]))
            X <- nrescale(X,ssnew=ssnew)
            X <- aperm(array(X,dim=xdim[c(2,1,3,4)]),perm=c(2,1,3,4))
          } else if(nmode==3L){
            X <- matrix(aperm(X,perm=c(3,1,2,4)),xdim[3],prod(xdim[c(1,2,4)]))
            X <- nrescale(X,ssnew=ssnew)
            X <- aperm(array(X,dim=xdim[c(3,1,2,4)]),perm=c(2,3,1,4))
          } else {
            X <- matrix(aperm(X,perm=c(4,1,2,3)),xdim[4],prod(xdim[1:3]))
            X <- nrescale(X,ssnew=ssnew)
            X <- aperm(array(X,dim=xdim[c(4,1,2,3)]),perm=c(2,3,4,1))
          }
        } else{stop("nrescale only rescales 2-, 3-, and 4-way arrays.")}
        
      } else if(is.list(X)){
        X <- mapply(nrescale,X,MoreArgs=list(nmode=nmode,ssnew=ssnew),SIMPLIFY=FALSE)
      } else{stop("Input X must be array or list with array elements.")}
      
    } # end if(nmode==0L)
    
    X
    
  }