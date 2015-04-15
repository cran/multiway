parafac <- 
  function(X,nfac,nstart=10,const=NULL,maxit=500,
           Bfixed=NULL,Cfixed=NULL,Dfixed=NULL,
           Bstart=NULL,Cstart=NULL,Dstart=NULL,
           ctol=10^-7,parallel=FALSE,cl=NULL){
    # 3-way or 4-way Parallel Factor Analysis (Parafac)
    # via alternating least squares (ALS) with optional constraints
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: April 9, 2015
    
    # check 'X' input
    xdim <- dim(X)
    lxdim <- length(xdim)
    if(lxdim<3L | lxdim>4L){stop("Input 'X' must be 3-way or 4-way array")}
    if(any(is.na(X)) | any(is.nan(X)) | any(is.infinite(X))){stop("Input 'X' cannot contain NA, NaN, or Inf values")}
    xcx <- sumsq(X)
    
    # check 'nfac' and 'nstart' inputs
    nfac <- as.integer(nfac[1])
    if(nfac<1L){stop("Input 'nfac' must be positive integer")}
    nstart <- as.integer(nstart[1])
    if(nstart<1L){stop("Input 'nstart' must be positive integer")}
    
    # check 'maxit' and 'ctol' inputs
    maxit <- as.integer(maxit[1])
    if(maxit<1L){stop("Input 'maxit' must be positive integer")}
    ctol <- as.numeric(ctol[1])
    if(ctol<0L){stop("Input 'ctol' must be positive numeric")}
    
    # check 'const' input
    if(is.null(const)){
      const <- rep(0L,lxdim)
    } else {
      const <- as.integer(const)
      if(length(const)!=lxdim){stop(paste("Input 'const' must be ",lxdim," element vector specifying constraint for each mode"))}
      if(min(const)<0L | max(const)>2L){stop("'const[j]' must be 0 (unconstrained), 1 (orthogonal), or 2 (non-negative)")}
    }
    
    # check 'Bfixed' and 'Cfixed' inputs
    if(!is.null(Bfixed)){
      Bfixed <- as.matrix(Bfixed)
      if(nrow(Bfixed)!=xdim[2]){stop("Input 'Bfixed' must have the same number of rows as dim(X)[2]")}
      if(ncol(Bfixed)!=nfac){stop("Input 'Bfixed' must have 'nfac' columns")}
    }
    if(!is.null(Cfixed)){
      Cfixed <- as.matrix(Cfixed)
      if(nrow(Cfixed)!=xdim[3]){stop("Input 'Cfixed' must have the same number of rows as dim(X)[3]")}
      if(ncol(Cfixed)!=nfac){stop("Input 'Cfixed' must have 'nfac' columns")}
    }
    
    # check 'Bstart' and 'Cstart' inputs
    if(!is.null(Bstart)){
      Bstart <- as.matrix(Bstart)
      if(nrow(Bstart)!=xdim[2]){stop("Input 'Bstart' must have the same number of rows as dim(X)[2]")}
      if(ncol(Bstart)!=nfac){stop("Input 'Bstart' must have 'nfac' columns")}
    }
    if(!is.null(Cstart)){
      Cstart <- as.matrix(Cstart)
      if(nrow(Cstart)!=xdim[3]){stop("Input 'Cstart' must have the same number of rows as dim(X)[3]")}
      if(ncol(Cstart)!=nfac){stop("Input 'Cstart' must have 'nfac' columns")}
    }
    
    # check 'parallel' and 'cl' inputs
    if(parallel && !any(class(cl)=="cluster")){
      stop("Input 'cl' must be cluster (created by makeCluster) when parallel=TRUE \n  See examples in documentation:  ?parafac")
    }
    
    # parafac fitting
    if(lxdim==3L){
      if(parallel){
        nstartlist <- vector("list",nstart)
        nstartlist[1:nstart] <- nfac
        pfaclist <- parLapply(cl=cl,X=nstartlist,fun="parafac_3way",data=X,xcx=xcx,
                              const=const,maxit=maxit,ctol=ctol,Bfixed=Bfixed,Cfixed=Cfixed,
                              Bstart=Bstart,Cstart=Cstart)
      } else {
        pfaclist <- vector("list",nstart)
        for(j in 1:nstart){
          pfaclist[[j]] <- parafac_3way(X,nfac,xcx,const,maxit,ctol,Bfixed,Cfixed,Bstart,Cstart)
        }
      }
      widx <- which.max(sapply(pfaclist,function(x) x$Rsq))
      pfac <- c(pfaclist[[widx]],list(const=const))
      class(pfac) <- "parafac"
      return(pfac)
    } else if(lxdim==4L){
      # check 'Dfixed' and 'Dstart' inputs
      if(!is.null(Dfixed)){
        Dfixed <- as.matrix(Dfixed)
        if(nrow(Dfixed)!=xdim[4]){stop("Input 'Dfixed' must have the same number of rows as dim(X)[4]")}
        if(ncol(Dfixed)!=nfac){stop("Input 'Dfixed' must have 'nfac' columns")}
      }
      if(!is.null(Dstart)){
        Dstart <- as.matrix(Dstart)
        if(nrow(Dstart)!=xdim[4]){stop("Input 'Dstart' must have the same number of rows as dim(X)[4]")}
        if(ncol(Dstart)!=nfac){stop("Input 'Dstart' must have 'nfac' columns")}
      }
      if(parallel){
        nstartlist <- vector("list",nstart)
        nstartlist[1:nstart] <- nfac
        pfaclist <- parLapply(cl=cl,X=nstartlist,fun="parafac_4way",data=X,xcx=xcx,
                              const=const,maxit=maxit,ctol=ctol,Bfixed=Bfixed,Cfixed=Cfixed,
                              Dfixed=Dfixed,Bstart=Bstart,Cstart=Cstart,Dstart=Dstart)
      } else {
        pfaclist <- vector("list",nstart)
        for(j in 1:nstart){
          pfaclist[[j]] <- parafac_4way(X,nfac,xcx,const,maxit,ctol,Bfixed,Cfixed,Dfixed,Bstart,Cstart,Dstart)
        }
      }
      widx <- which.max(sapply(pfaclist,function(x) x$Rsq))
      pfac <- c(pfaclist[[widx]],list(const=const))
      class(pfac) <- "parafac"
      return(pfac)
    } 
    
  }