ncenter <- 
  function(X,nmode=1){
    # Center n-th Dimension of Array
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: April 9, 2015
    
    nmode <- as.integer(nmode[1])
    if(nmode<1L)(stop("nmode must be positive integer."))
    if(is.array(X)){
      xdim <- dim(X)
      if(nmode>length(xdim)){stop("nmode must be less than maximum mode.")}
      if(length(xdim)==2L){
        if(nmode==1L){
          xmns <- colMeans(X)
          X <- X - matrix(xmns,xdim[1],xdim[2],byrow=TRUE)
        } else {
          xmns <- rowMeans(X)
          X <- X - matrix(xmns,xdim[1],xdim[2])
        }
      } else if(length(xdim)==3L){
        if(nmode==1L){
          X <- matrix(X,xdim[1],xdim[2]*xdim[3])
          X <- ncenter(X)
          X <- array(X,dim=xdim)
        } else if(nmode==2L){
          X <- matrix(aperm(X,perm=c(2,1,3)),xdim[2],xdim[1]*xdim[3])
          X <- ncenter(X)
          X <- aperm(array(X,dim=xdim[c(2,1,3)]),perm=c(2,1,3))
        } else {
          X <- matrix(aperm(X,perm=c(3,1,2)),xdim[3],xdim[1]*xdim[2])
          X <- ncenter(X)
          X <- aperm(array(X,dim=xdim[c(3,1,2)]),perm=c(2,3,1))
        }
      } else if(length(xdim)==4L){
        if(nmode==1L){
          X <- matrix(X,xdim[1],prod(xdim[2:4]))
          X <- ncenter(X)
          X <- array(X,dim=xdim)
        } else if(nmode==2L){
          X <- matrix(aperm(X,perm=c(2,1,3,4)),xdim[2],prod(xdim[c(1,3,4)]))
          X <- ncenter(X)
          X <- aperm(array(X,dim=xdim[c(2,1,3,4)]),perm=c(2,1,3,4))
        } else if(nmode==3L){
          X <- matrix(aperm(X,perm=c(3,1,2,4)),xdim[3],prod(xdim[c(1,2,4)]))
          X <- ncenter(X)
          X <- aperm(array(X,dim=xdim[c(3,1,2,4)]),perm=c(2,3,1,4))
        } else {
          X <- matrix(aperm(X,perm=c(4,1,2,3)),xdim[4],prod(xdim[1:3]))
          X <- ncenter(X)
          X <- aperm(array(X,dim=xdim[c(4,1,2,3)]),perm=c(2,3,4,1))
        }
      } else {stop("ncenter only centers 2-, 3-, and 4-way arrays.")}
    } else if(is.list(X)){
      X <- mapply(ncenter,X,MoreArgs=list(nmode=nmode),SIMPLIFY=FALSE)
    } else{stop("Input X must be array or list with array elements.")}
    
    X
    
  }