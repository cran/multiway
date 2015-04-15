indscal <- 
  function(X,nfac,nstart=10,const=NULL,maxit=500,
           type=c("dissimilarity","similarity"),
           ctol=10^-7,parallel=FALSE,cl=NULL){
    # Individual Differences Scaling (INDSCAL)
    # via alternating least squares (ALS) with optional constraints
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: April 9, 2015
    
    # check 'X' input
    if(is.list(X)){
      X <- lapply(X,as.matrix)
      px <- ncol(X[[1]])
      nx <- length(X)
      if(sum((sapply(X,dim)-matrix(px,2,nx))^2)>0L){
        stop("Input 'X' must be list of square matrices with same dimension.")
      }
      xdim <- c(px,px,nx)
      X <- array(unlist(X),dim=xdim)
    } else{
      xdim <- dim(X)
      if(length(xdim)!=3L){stop("Input 'X' must be 3-way array")}
      if(xdim[1]!=xdim[2]){stop("Input 'X' must be 3-way array with dim(X)[1]==dim(X)[2].")}
    }
    if(any(is.na(X)) | any(is.nan(X)) | any(is.infinite(X))){stop("Input 'X' cannot contain NA, NaN, or Inf values")}
    if(type[1]=="dissimilarity") { X <- array(apply(X,3,ed2sp),dim=xdim) }
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
      const <- c(0L,0L)
    } else {
      const <- as.integer(const[1:2])
      if(!any(const[1]==c(0L,1L))){stop("'const[1]' must be 0 (unconstrained) or 1 (orthogonal)")}
      if(!any(const[2]==c(0L,2L))){stop("'const[2]' must be 0 (unconstrained) or 2 (nonnegative)")}
    }
    
    # check 'parallel' and 'cl' inputs
    if(parallel && !any(class(cl)=="cluster")){
      stop("Input 'cl' must be cluster (created by makeCluster) when parallel=TRUE \n  See examples in documentation:  ?parafac")
    }
    
    # indscal fitting
    if(parallel){
      nstartlist <- vector("list",nstart)
      nstartlist[1:nstart] <- nfac
      pfaclist <- parLapply(cl=cl,X=nstartlist,fun="parafac_3way",data=X,xcx=xcx,
                            const=c(rep(const[1],2),const[2]),maxit=maxit,ctol=ctol)
    } else {
      pfaclist <- vector("list",nstart)
      for(j in 1:nstart){
        pfaclist[[j]] <- parafac_3way(X,nfac,xcx,c(rep(const[1],2),const[2]),maxit,ctol)
      }
    }
    widx <- which.max(sapply(pfaclist,function(x) x$Rsq))
    pfac <- pfaclist[[widx]]
    adg <- colMeans(pfac$A^2)
    pfac$C <- pfac$C%*%(diag(nfac)*(adg^0.5))
    if(const[2]==0L){
      csg <- sign(colSums(pfac$C^3))
      pfac$C <- pfac$C%*%(diag(nfac)*(csg))
      if(any(sign(pfac$C)<0L)){warning("Unconstrained ALS produced negative C weights. \nTry refitting with nonnegativity for Mode C, e.g., const=c(0,2)")}
    }
    bsg <- sign(colSums(pfac$B^3))
    pfac$B <- pfac$B%*%(diag(nfac)*(bsg))
    inds <- list(B=pfac$B,C=pfac$C,Rsq=pfac$Rsq,iter=pfac$iter,
                 cflag=pfac$cflag,const=const,strain=(1-pfac$Rsq)*sum(X^2))
    class(inds) <- "indscal"
    return(inds)
    
  }