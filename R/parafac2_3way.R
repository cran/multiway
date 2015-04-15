parafac2_3way <- 
  function(data,nfac,xcx=sumsq(data),const=rep(0L,3),
           maxit=500,ctol=10^-7,Bfixed=NULL,Cfixed=NULL,
           Bstart=NULL,Cstart=NULL){
    # 3-way Parallel Factor Analysis 2 (Parafac2)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: April 11, 2015
    
    ### initialize Khatri-Rao product matrices
    xdims <- rep(NA,3)
    xdims[2] <- ncol(data[[1]])
    xdims[3] <- length(data)
    if(is.null(Cfixed)){BkrA <- matrix(0,nfac*xdims[2],nfac)}
    if(is.null(Bfixed)){CkrA <- matrix(0,nfac*xdims[3],nfac)}
    CkrB <- matrix(0,xdims[2]*xdims[3],nfac)
    
    ### initialize stuff for Mode A update
    Rknew <- vector("list",xdims[3])
    Xtilde <- array(0,dim=c(nfac,xdims[2],xdims[3]))
    
    ### initialize parameter matrices
    if(const[1]==0L){
      Gold <- crossprod(matrix(rnorm(nfac^2),nfac,nfac))
    } else if(const[1]==1L){
      Gold <- diag(nfac)
    } 
    if(is.null(Bfixed)){
      if(!is.null(Bstart)){
        Bold <- Bstart
      } else if(const[2]==0L){
        Bold <- matrix(rnorm(xdims[2]*nfac),xdims[2],nfac)
      } else if(const[2]==1L){
        Bold <- svd(matrix(rnorm(xdims[2]*nfac),xdims[2],nfac),nu=nfac,nv=0)$u
      } else if(const[2]==2L){
        Bold <- Bnew <- matrix(runif(xdims[2]*nfac),xdims[2],nfac)
      }
    } else {Bold <- Bnew <- Bfixed}
    if(is.null(Cfixed)){
      if(!is.null(Cstart)){
        Cold <- Cstart
      } else if(const[3]==0L){
        Cold <- matrix(rnorm(xdims[3]*nfac),xdims[3],nfac)
      } else if(const[3]==1L){
        Cold <- svd(matrix(rnorm(xdims[3]*nfac),xdims[3],nfac),nu=nfac,nv=0)$u
      } else if(const[3]==2L){
        Cold <- Cnew <- matrix(runif(xdims[3]*nfac),xdims[3],nfac)
      }
    } else {Cold <- Cnew <- Cfixed}
    
    ### iterative update of matrices
    vtol <- sseold <- xcx
    iter <- 0
    cflag <- NA
    while(vtol>ctol && iter<maxit) {
      
      ## Step 1: update mode A weights
      # 1a: update orthogonal projections
      for(kk in 1:xdims[3]){
        xsvd <- svd(data[[kk]]%*%Bold%*%tcrossprod((diag(nfac)*Cold[kk,]),Gold))
        Rknew[[kk]] <- tcrossprod(xsvd$u,xsvd$v)
        Xtilde[,,kk] <- crossprod(Rknew[[kk]],data[[kk]])
      }
      # 1b: update correlation matrix
      if(const[1]==0L){
        Xa <- matrix(Xtilde,nfac,xdims[2]*xdims[3])
        for(u in 1:nfac){CkrB[,u] <- kronecker(Cold[,u],Bold[,u])}
        Gnew <- Xa%*%CkrB%*%smpower(crossprod(CkrB),-1)
      } else{Gnew <- diag(nfac)}
      
      ## Step 2: update mode B weights
      if(is.null(Bfixed)){
        Xb <- matrix(aperm(Xtilde,perm=c(2,1,3)),xdims[2],nfac*xdims[3])
        for(u in 1:nfac){CkrA[,u] <- kronecker(Cold[,u],Gnew[,u])}
        if(const[2]==0L){
          Bnew <- Xb%*%CkrA%*%smpower(crossprod(CkrA),-1)
        } else if(const[2]==1L) {
          Zmat <- Xb%*%CkrA
          Bnew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
        } else if(const[2]==2L) {
          cpmat <- crossprod(CkrA)
          for(ii in 1:xdims[2]){Bnew[ii,] <- fnnls(cpmat,crossprod(CkrA,Xb[ii,]))}
          if(any(colSums(Bnew)==0)){Bnew <- Bold;  cflag <- 2}
        }
      }
      
      ## Step 3: update mode C weights
      if(is.null(Cfixed)){
        Xc <- matrix(aperm(Xtilde,perm=c(3,1,2)),xdims[3],nfac*xdims[2])
        for(u in 1:nfac){BkrA[,u] <- kronecker(Bnew[,u],Gnew[,u])}
        if(const[3]==0L){
          Cnew <- Xc%*%BkrA%*%smpower(crossprod(BkrA),-1)
        } else if(const[3]==1L) {
          Zmat <- Xc%*%BkrA
          Cnew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
        } else if(const[3]==2L) {
          cpmat <- crossprod(BkrA)
          for(ii in 1:xdims[3]){Cnew[ii,] <- fnnls(cpmat,crossprod(BkrA,Xc[ii,]))}
          if(any(colSums(Cnew)==0)){Cnew <- Cold;  cflag <- 2}
        }
      }
      
      ## Step 4: check for convergence
      ssenew <- 0
      for(kk in 1:xdims[3]){
        ssenew <- ssenew + sum((data[[kk]]-tcrossprod(Rknew[[kk]]%*%Gnew%*%(diag(nfac)*Cnew[kk,]),Bnew))^2)
      }
      vtol <- (sseold-ssenew)/sseold
      Gold <- Gnew
      Bold <- Bnew
      Cold <- Cnew
      sseold <- ssenew
      iter <- iter + 1
      
    }  # end while(vtol>ctol && iter<maxit)
    
    ### put the scale in Mode A      
    bdg <- colMeans(Bnew^2)
    Bnew <- Bnew%*%(diag(nfac)*(bdg^-0.5))
    cdg <- colMeans(Cnew^2)
    Cnew <- Cnew%*%(diag(nfac)*(cdg^-0.5))
    Gnew <- Gnew%*%(diag(nfac)*((bdg*cdg)^0.5))
    
    ### order the solution
    if(any(c(!is.null(Bfixed),!is.null(Cfixed)))){
      if(!is.null(Bfixed)){
        bdgfx <- colMeans(Bfixed^2)
        bdg <- colMeans(Bnew^2)
        Bnew <- Bnew%*%(diag(nfac)*((bdgfx/bdg)^0.5))
        Gnew <- Gnew%*%(diag(nfac)*((bdg/bdgfx)^0.5))
        bsgfx <- sign(colSums(Bfixed^3))
        bsg <- sign(colSums(Bnew^3))
        Bnew <- Bnew%*%(diag(nfac)*(bsg*bsgfx))
        Gnew <- Gnew%*%(diag(nfac)*(bsg*bsgfx))
      }
      if(!is.null(Cfixed)){
        cdgfx <- colMeans(Cfixed^2)
        cdg <- colMeans(Cnew^2)
        Cnew <- Cnew%*%(diag(nfac)*((cdgfx/cdg)^0.5))
        Gnew <- Gnew%*%(diag(nfac)*((cdg/cdgfx)^0.5))
        csgfx <- sign(colSums(Cfixed^3))
        csg <- sign(colSums(Cnew^3))
        Cnew <- Cnew%*%(diag(nfac)*(csg*csgfx))
        Gnew <- Gnew%*%(diag(nfac)*(csg*csgfx))
      } 
    } else {
      fordr <- order(colSums(Gnew^2),decreasing=TRUE)
      Gnew <- as.matrix(Gnew[,fordr])
      Bnew <- as.matrix(Bnew[,fordr])
      Cnew <- as.matrix(Cnew[,fordr])
    }
    
    ### collect results
    Rsq <- 1 - ssenew/xcx
    if(is.na(cflag)){
      if(vtol<=ctol){cflag <- 0} else {cflag <- 1}
    }
    pfac <- list(A=list(H=Rknew,G=Gnew),B=Bnew,C=Cnew,Rsq=Rsq,iter=iter,cflag=cflag)
    return(pfac)
    
  }