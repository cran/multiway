parafac2_4way <- 
  function(data,nfac,xcx=sumsq(data),const=rep(0L,4),
           maxit=500,ctol=10^-7,Bfixed=NULL,Cfixed=NULL,
           Dfixed=NULL,Bstart=NULL,Cstart=NULL,Dstart=NULL){
    # 4-way Parallel Factor Analysis 2 (Parafac2)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: April 11, 2015
    
    ### initialize Khatri-Rao product matrices
    xdims <- rep(NA,4)
    xdims[2] <- dim(data[[1]])[2]
    xdims[3] <- dim(data[[1]])[3]
    xdims[4] <- length(data)
    if(is.null(Dfixed)){CBkrA <- matrix(0,nfac*xdims[2]*xdims[3],nfac)}
    if(is.null(Cfixed)){DBkrA <- matrix(0,nfac*xdims[2]*xdims[4],nfac)}
    if(is.null(Bfixed)){DCkrA <- matrix(0,nfac*xdims[3]*xdims[4],nfac)}
    DCkrB <- matrix(0,xdims[2]*xdims[3]*xdims[4],nfac)
    CkrB <- matrix(0,xdims[2]*xdims[3],nfac)
    
    ### initialize stuff for Mode A update
    Rknew <- vector("list",xdims[4])
    Xtilde <- array(0,dim=c(nfac,xdims[2],xdims[3],xdims[4]))
    
    ### reshape raw data
    for(kk in 1:xdims[4]){
      mdim <- dim(data[[kk]])
      data[[kk]] <- matrix(data[[kk]],mdim[1],xdims[2]*xdims[3])
    }
    
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
    if(is.null(Dfixed)){
      if(!is.null(Dstart)){
        Dold <- Dstart
      } else if(const[4]==0L){
        Dold <- matrix(rnorm(xdims[4]*nfac),xdims[4],nfac)
      } else if(const[4]==1L){
        Dold <- svd(matrix(rnorm(xdims[4]*nfac),xdims[4],nfac),nu=nfac,nv=0)$u
      } else if(const[4]==2L){
        Dold <- Dnew <- matrix(runif(xdims[4]*nfac),xdims[4],nfac)
      }
    } else {Dold <- Dnew <- Dfixed}
    
    ### iterative update of matrices
    vtol <- sseold <- xcx
    iter <- 0
    cflag <- NA
    while(vtol>ctol && iter<maxit) {
      
      ## Step 1: update mode A weights
      # 1a: update orthogonal projections
      for(u in 1:nfac){CkrB[,u] <- kronecker(Cold[,u],Bold[,u])}
      for(kk in 1:xdims[4]){
        xsvd <- svd(data[[kk]]%*%CkrB%*%tcrossprod((diag(nfac)*Dold[kk,]),Gold))
        Rknew[[kk]] <- tcrossprod(xsvd$u,xsvd$v)
        Xtilde[,,,kk] <- array(crossprod(Rknew[[kk]],data[[kk]]),dim=c(nfac,xdims[2],xdims[3]))
      }
      # 1b: update correlation matrix
      if(const[1]==0L){
        Xa <- matrix(Xtilde,nfac,xdims[2]*xdims[3]*xdims[4])
        for(u in 1:nfac){DCkrB[,u] <- kronecker(Dold[,u],kronecker(Cold[,u],Bold[,u]))}
        Gnew <- Xa%*%DCkrB%*%smpower(crossprod(DCkrB),-1)
      } else{Gnew <- diag(nfac)}
      
      ## Step 2: update mode B weights
      if(is.null(Bfixed)){
        Xb <- matrix(aperm(Xtilde,perm=c(2,1,3,4)),xdims[2],nfac*xdims[3]*xdims[4])
        for(u in 1:nfac){DCkrA[,u] <- kronecker(Dold[,u],kronecker(Cold[,u],Gnew[,u]))}
        if(const[2]==0L){
          Bnew <- Xb%*%DCkrA%*%smpower(crossprod(DCkrA),-1)
        } else if(const[2]==1L) {
          Zmat <- Xb%*%DCkrA
          Bnew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
        } else if(const[2]==2L) {
          cpmat <- crossprod(DCkrA)
          for(ii in 1:xdims[2]){Bnew[ii,] <- fnnls(cpmat,crossprod(DCkrA,Xb[ii,]))}
          if(any(colSums(Bnew)==0)){Bnew <- Bold;  cflag <- 2}
        }
      }
      
      ## Step 3: update mode C weights
      if(is.null(Cfixed)){
        Xc <- matrix(aperm(Xtilde,perm=c(3,1,2,4)),xdims[3],nfac*xdims[2]*xdims[4])
        for(u in 1:nfac){DBkrA[,u] <- kronecker(Dold[,u],kronecker(Bnew[,u],Gnew[,u]))}
        if(const[3]==0L){
          Cnew <- Xc%*%DBkrA%*%smpower(crossprod(DBkrA),-1)
        } else if(const[3]==1L) {
          Zmat <- Xc%*%DBkrA
          Cnew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
        } else if(const[3]==2L) {
          cpmat <- crossprod(DBkrA)
          for(ii in 1:xdims[3]){Cnew[ii,] <- fnnls(cpmat,crossprod(DBkrA,Xc[ii,]))}
          if(any(colSums(Cnew)==0)){Cnew <- Cold;  cflag <- 2}
        }
      }
      
      ## Step 4: update mode D weights
      if(is.null(Dfixed)){
        Xd <- matrix(aperm(Xtilde,perm=c(4,1,2,3)),xdims[4],nfac*xdims[2]*xdims[3])
        for(u in 1:nfac){CBkrA[,u] <- kronecker(Cnew[,u],kronecker(Bnew[,u],Gnew[,u]))}
        if(const[4]==0L){
          Dnew <- Xd%*%CBkrA%*%smpower(crossprod(CBkrA),-1)
        } else if(const[4]==1L) {
          Zmat <- Xd%*%CBkrA
          Dnew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
        } else if(const[4]==2L) {
          cpmat <- crossprod(CBkrA)
          for(ii in 1:xdims[4]){Dnew[ii,] <- fnnls(cpmat,crossprod(CBkrA,Xd[ii,]))}
          if(any(colSums(Dnew)==0)){Dnew <- Dold;  cflag <- 2}
        }
      }
      
      ## Step 5: check for convergence
      for(u in 1:nfac){CkrB[,u] <- kronecker(Cnew[,u],Bnew[,u])}
      ssenew <- 0
      for(kk in 1:xdims[4]){
        ssenew <- ssenew + sum((data[[kk]]-tcrossprod(Rknew[[kk]]%*%Gnew%*%(diag(nfac)*Dnew[kk,]),CkrB))^2)
      }
      vtol <- (sseold-ssenew)/sseold
      Gold <- Gnew
      Bold <- Bnew
      Cold <- Cnew
      Dold <- Dnew
      sseold <- ssenew
      iter <- iter + 1
      
    }  # end while(vtol>ctol && iter<maxit)
    
    ### put the scale in Mode A      
    bdg <- colMeans(Bnew^2)
    Bnew <- Bnew%*%(diag(nfac)*(bdg^-0.5))
    cdg <- colMeans(Cnew^2)
    Cnew <- Cnew%*%(diag(nfac)*(cdg^-0.5))
    ddg <- colMeans(Dnew^2)
    Dnew <- Dnew%*%(diag(nfac)*(ddg^-0.5))
    Gnew <- Gnew%*%(diag(nfac)*((bdg*cdg*ddg)^0.5))
    
    ### order the solution
    if(any(c(!is.null(Bfixed),!is.null(Cfixed),!is.null(Dfixed)))){
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
      if(!is.null(Dfixed)){
        ddgfx <- colMeans(Dfixed^2)
        ddg <- colMeans(Dnew^2)
        Dnew <- Dnew%*%(diag(nfac)*((ddgfx/ddg)^0.5))
        Gnew <- Gnew%*%(diag(nfac)*((ddg/ddgfx)^0.5))
        dsgfx <- sign(colSums(Dfixed^3))
        dsg <- sign(colSums(Dnew^3))
        Dnew <- Dnew%*%(diag(nfac)*(dsg*dsgfx))
        Gnew <- Gnew%*%(diag(nfac)*(dsg*dsgfx))
      }
    } else {
      fordr <- order(colSums(Gnew^2),decreasing=TRUE)
      Gnew <- as.matrix(Gnew[,fordr])
      Bnew <- as.matrix(Bnew[,fordr])
      Cnew <- as.matrix(Cnew[,fordr])
      Dnew <- as.matrix(Dnew[,fordr])
    }
    
    ### collect results
    Rsq <- 1 - ssenew/xcx
    if(is.na(cflag)){
      if(vtol<=ctol){cflag <- 0} else {cflag <- 1}
    }
    pfac <- list(A=list(H=Rknew,G=Gnew),B=Bnew,C=Cnew,D=Dnew,Rsq=Rsq,iter=iter,cflag=cflag)
    return(pfac)
    
  }