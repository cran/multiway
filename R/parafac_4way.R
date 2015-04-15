parafac_4way <-
  function(data,nfac,xcx=sumsq(data),const=rep(0L,4),
           maxit=500,ctol=10^-7,Bfixed=NULL,Cfixed=NULL,
           Dfixed=NULL,Bstart=NULL,Cstart=NULL,Dstart=NULL){
    # 4-way Parallel Factor Analysis (Parafac)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: April 9, 2015
    
    ### initialize Khatri-Rao product matrices
    xdims <- dim(data)
    if(is.null(Dfixed)){CBkrA <- matrix(0,xdims[1]*xdims[2]*xdims[3],nfac)}
    if(is.null(Cfixed)){DBkrA <- matrix(0,xdims[1]*xdims[2]*xdims[4],nfac)}
    if(is.null(Bfixed)){DCkrA <- matrix(0,xdims[1]*xdims[3]*xdims[4],nfac)}
    DCkrB <- matrix(0,xdims[2]*xdims[3]*xdims[4],nfac)
    
    ### initialize reshaped data matrices
    Xa <- matrix(data,xdims[1],xdims[2]*xdims[3]*xdims[4])
    if(is.null(Bfixed)){Xb <- matrix(aperm(data,perm=c(2,1,3,4)),xdims[2],xdims[1]*xdims[3]*xdims[4])}
    if(is.null(Cfixed)){Xc <- matrix(aperm(data,perm=c(3,1,2,4)),xdims[3],xdims[1]*xdims[2]*xdims[4])}
    if(is.null(Dfixed)){Xd <- matrix(aperm(data,perm=c(4,1,2,3)),xdims[4],xdims[1]*xdims[2]*xdims[3])}
    rm(data)
    
    ### initialize parameter matrices
    if(const[1]==0L){
      Aold <- matrix(rnorm(xdims[1]*nfac),xdims[1],nfac)
    } else if(const[1]==1L){
      Aold <- svd(matrix(rnorm(xdims[1]*nfac),xdims[1],nfac),nu=nfac,nv=0)$u
    } else if(const[1]==2L){
      Aold <- Anew <- matrix(runif(xdims[1]*nfac),xdims[1],nfac)
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
    } else {Dold=Dnew=Dfixed}
    
    ### iterative update of matrices
    vtol <- sseold <- xcx
    iter <- 0
    cflag <- NA
    while(vtol>ctol && iter<maxit) {
      
      # Step 1: update mode A weights
      for(u in 1:nfac){DCkrB[,u] <- kronecker(Dold[,u],kronecker(Cold[,u],Bold[,u]))}
      if(const[1]==0L){
        Anew <- Xa%*%DCkrB%*%smpower(crossprod(DCkrB),-1)
      } else if(const[1]==1L) {
        Zmat <- Xa%*%DCkrB
        Anew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
      } else if(const[1]==2L) {
        cpmat <- crossprod(DCkrB)
        for(ii in 1:xdims[1]){Anew[ii,] <- fnnls(cpmat,crossprod(DCkrB,Xa[ii,]))}
        if(any(colSums(Anew)==0)){Anew <- Aold;  cflag <- 2}
      }
      
      # Step 2: update mode B weights
      if(is.null(Bfixed)){
        for(u in 1:nfac){DCkrA[,u] <- kronecker(Dold[,u],kronecker(Cold[,u],Anew[,u]))}
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
      
      # Step 3: update mode C weights
      if(is.null(Cfixed)){
        for(u in 1:nfac){DBkrA[,u] <- kronecker(Dold[,u],kronecker(Bnew[,u],Anew[,u]))}
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
      
      # Step 4: update mode D weights
      if(is.null(Dfixed)){
        for(u in 1:nfac){CBkrA[,u] <- kronecker(Cnew[,u],kronecker(Bnew[,u],Anew[,u]))}
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
      
      # Step 5: check for convergence
      for(u in 1:nfac){DCkrB[,u] <- kronecker(Dnew[,u],kronecker(Cnew[,u],Bnew[,u]))}
      ssenew <- sum((Xa-tcrossprod(Anew,DCkrB))^2)
      vtol <- (sseold-ssenew)/sseold
      Aold <- Anew
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
    Anew <- Anew%*%(diag(nfac)*((bdg*cdg*ddg)^0.5))
    
    ### order the solution
    if(any(c(!is.null(Bfixed),!is.null(Cfixed),!is.null(Dfixed)))){
      if(!is.null(Bfixed)){
        bdgfx <- colMeans(Bfixed^2)
        bdg <- colMeans(Bnew^2)
        Bnew <- Bnew%*%(diag(nfac)*((bdgfx/bdg)^0.5))
        Anew <- Anew%*%(diag(nfac)*((bdg/bdgfx)^0.5))
        bsgfx <- sign(colSums(Bfixed^3))
        bsg <- sign(colSums(Bnew^3))
        Bnew <- Bnew%*%(diag(nfac)*(bsg*bsgfx))
        Anew <- Anew%*%(diag(nfac)*(bsg*bsgfx))
      }
      if(!is.null(Cfixed)){
        cdgfx <- colMeans(Cfixed^2)
        cdg <- colMeans(Cnew^2)
        Cnew <- Cnew%*%(diag(nfac)*((cdgfx/cdg)^0.5))
        Anew <- Anew%*%(diag(nfac)*((cdg/cdgfx)^0.5))
        csgfx <- sign(colSums(Cfixed^3))
        csg <- sign(colSums(Cnew^3))
        Cnew <- Cnew%*%(diag(nfac)*(csg*csgfx))
        Anew <- Anew%*%(diag(nfac)*(csg*csgfx))
      }
      if(!is.null(Dfixed)){
        ddgfx <- colMeans(Dfixed^2)
        ddg <- colMeans(Dnew^2)
        Dnew <- Dnew%*%(diag(nfac)*((ddgfx/ddg)^0.5))
        Anew <- Anew%*%(diag(nfac)*((ddg/ddgfx)^0.5))
        dsgfx <- sign(colSums(Dfixed^3))
        dsg <- sign(colSums(Dnew^3))
        Dnew <- Dnew%*%(diag(nfac)*(dsg*dsgfx))
        Anew <- Anew%*%(diag(nfac)*(dsg*dsgfx))
      }
    } else {
      fordr <- order(colSums(Anew^2),decreasing=TRUE)
      Anew <- as.matrix(Anew[,fordr])
      Bnew <- as.matrix(Bnew[,fordr])
      Cnew <- as.matrix(Cnew[,fordr])
      Dnew <- as.matrix(Dnew[,fordr])
    }
    
    ### collect results
    Rsq <- 1 - ssenew/xcx
    if(is.na(cflag)){
      if(vtol<=ctol){cflag <- 0} else {cflag <- 1}
    }
    pfac <- list(A=Anew,B=Bnew,C=Cnew,D=Dnew,Rsq=Rsq,iter=iter,cflag=cflag)
    return(pfac)
    
  }