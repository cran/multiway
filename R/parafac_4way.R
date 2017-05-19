parafac_4way <-
  function(data,nfac,xcx=sumsq(data),const=rep(0L,4),
           maxit=500,ctol=1e-4,Bfixed=NULL,Cfixed=NULL,
           Dfixed=NULL,Bstart=NULL,Cstart=NULL,Dstart=NULL,
           Bstruc=NULL,Cstruc=NULL,Dstruc=NULL,
           control=const.control(const)){
    # 4-way Parallel Factor Analysis (Parafac)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 16, 2017
    
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
    } else if(const[1]==3L){
      Aold <- Anew <- replicate(nfac, rfunc(1:xdims[1], df=control$df[1], degree=control$degree[1], type="uni", nonneg=control$nonneg[1]))
      Adf <- rep(0, nfac)
    } else if(const[1]==4L){
      Aold <- Anew <- replicate(nfac, rfunc(1:xdims[1], df=control$df[1], degree=control$degree[1], type="mon", nonneg=control$nonneg[1]))
      Adf <- rep(0, nfac)
    } else if(const[1]==5L){
      Aold <- Anew <- replicate(nfac, rfunc(1:xdims[1], df=control$df[1], degree=control$degree[1], type="per", nonneg=control$nonneg[1]))
      Adf <- rep(0, nfac)
      if(!control$nonneg[1]) {
        Ra <- tcrossprod(svd(cbind(1,MsplineBasis(1:xdims[1], df=control$df[1], degree=control$degree[1], intercept=FALSE)$X[,-control$df[1]]))$u)
        Adf <- rep(control$df[1], nfac)
      }
    } else if(const[1]==6L){
      Aold <- Anew <- replicate(nfac, rfunc(1:xdims[1], df=control$df[1], degree=control$degree[1], type="smo", nonneg=control$nonneg[1]))
      Adf <- rep(0, nfac)
      if(!control$nonneg[1]) {
        Ra <- tcrossprod(svd(MsplineBasis(1:xdims[1], df=control$df[1], degree=control$degree[1], intercept=TRUE)$X)$u)
        Adf <- rep(control$df[1], nfac)
      }
    }
    if(is.null(Bfixed)){
      if(!is.null(Bstart)){
        Bold <- Bnew <- Bstart
      } else if(const[2]==0L){
        Bold <- matrix(rnorm(xdims[2]*nfac),xdims[2],nfac)
      } else if(const[2]==1L){
        Bold <- svd(matrix(rnorm(xdims[2]*nfac),xdims[2],nfac),nu=nfac,nv=0)$u
      } else if(const[2]==2L){
        Bold <- Bnew <- matrix(runif(xdims[2]*nfac),xdims[2],nfac)
      } else if(const[2]==3L){
        Bold <- Bnew <- replicate(nfac, rfunc(1:xdims[2], df=control$df[2], degree=control$degree[2], type="uni", nonneg=control$nonneg[2]))
        Bdf <- rep(0, nfac)
      } else if(const[2]==4L){
        Bold <- Bnew <- replicate(nfac, rfunc(1:xdims[2], df=control$df[2], degree=control$degree[2], type="mon", nonneg=control$nonneg[2]))
        Bdf <- rep(0, nfac)
      } else if(const[2]==5L){
        Bold <- Bnew <- replicate(nfac, rfunc(1:xdims[2], df=control$df[2], degree=control$degree[2], type="per", nonneg=control$nonneg[2]))
        Bdf <- rep(0, nfac)
        if(!control$nonneg[2]) {
          Rb <- tcrossprod(svd(cbind(1,MsplineBasis(1:xdims[2], df=control$df[2], degree=control$degree[2], intercept=FALSE)$X[,-control$df[2]]))$u)
          Bdf <- rep(control$df[2], nfac)
        }
      } else if(const[2]==6L){
        Bold <- Bnew <- replicate(nfac, rfunc(1:xdims[2], df=control$df[2], degree=control$degree[2], type="smo", nonneg=control$nonneg[2]))
        Bdf <- rep(0, nfac)
        if(!control$nonneg[2]) {
          Rb <- tcrossprod(svd(MsplineBasis(1:xdims[2], df=control$df[2], degree=control$degree[2], intercept=TRUE)$X)$u)
          Bdf <- rep(control$df[2], nfac)
        }
      }
      if(!is.null(Bstruc)) Bold <- Bnew <- Bold * Bstruc
    } else {
      Bold <- Bnew <- Bfixed
    }
    if(is.null(Cfixed)){
      if(!is.null(Cstart)){
        Cold <- Cnew <- Cstart
      } else if(const[3]==0L){
        Cold <- matrix(rnorm(xdims[3]*nfac),xdims[3],nfac)
      } else if(const[3]==1L){
        Cold <- svd(matrix(rnorm(xdims[3]*nfac),xdims[3],nfac),nu=nfac,nv=0)$u
      } else if(const[3]==2L){
        Cold <- Cnew <- matrix(runif(xdims[3]*nfac),xdims[3],nfac)
      } else if(const[3]==3L){
        Cold <- Cnew <- replicate(nfac, rfunc(1:xdims[3], df=control$df[3], degree=control$degree[3], type="uni", nonneg=control$nonneg[3]))
        Cdf <- rep(0, nfac)
      } else if(const[3]==4L){
        Cold <- Cnew <- replicate(nfac, rfunc(1:xdims[3], df=control$df[3], degree=control$degree[3], type="mon", nonneg=control$nonneg[3]))
        Cdf <- rep(0, nfac)
      } else if(const[3]==5L){
        Cold <- Cnew <- replicate(nfac, rfunc(1:xdims[3], df=control$df[3], degree=control$degree[3], type="per", nonneg=control$nonneg[3]))
        Cdf <- rep(0, nfac)
        if(!control$nonneg[3]) {
          Rc <- tcrossprod(svd(cbind(1,MsplineBasis(1:xdims[3], df=control$df[3], degree=control$degree[3], intercept=FALSE)$X[,-control$df[3]]))$u)
          Cdf <- rep(control$df[3], nfac)
        }
      } else if(const[3]==6L){
        Cold <- Cnew <- replicate(nfac, rfunc(1:xdims[3], df=control$df[3], degree=control$degree[3], type="smo", nonneg=control$nonneg[3]))
        Cdf <- rep(0, nfac)
        if(!control$nonneg[3]) {
          Rc <- tcrossprod(svd(MsplineBasis(1:xdims[3], df=control$df[3], degree=control$degree[3], intercept=TRUE)$X)$u)
          Cdf <- rep(control$df[3], nfac)
        }
      }
      if(!is.null(Cstruc)) Cold <- Cnew <- Cold * Cstruc
    } else {
      Cold <- Cnew <- Cfixed
    }
    if(is.null(Dfixed)){
      if(!is.null(Dstart)){
        Dold <- Dnew <- Dstart
      } else if(const[4]==0L){
        Dold <- matrix(rnorm(xdims[4]*nfac),xdims[4],nfac)
      } else if(const[4]==1L){
        Dold <- svd(matrix(rnorm(xdims[4]*nfac),xdims[4],nfac),nu=nfac,nv=0)$u
      } else if(const[4]==2L){
        Dold <- Dnew <- matrix(runif(xdims[4]*nfac),xdims[4],nfac)
      } else if(const[4]==3L){
        Dold <- Dnew <- replicate(nfac, rfunc(1:xdims[4], df=control$df[4], degree=control$degree[4], type="uni", nonneg=control$nonneg[4]))
        Ddf <- rep(0, nfac)
      } else if(const[4]==4L){
        Dold <- Dnew <- replicate(nfac, rfunc(1:xdims[4], df=control$df[4], degree=control$degree[4], type="mon", nonneg=control$nonneg[4]))
        Ddf <- rep(0, nfac)
      } else if(const[4]==5L){
        Dold <- Dnew <- replicate(nfac, rfunc(1:xdims[4], df=control$df[4], degree=control$degree[4], type="per", nonneg=control$nonneg[4]))
        Ddf <- rep(0, nfac)
        if(!control$nonneg[4]) {
          Rd <- tcrossprod(svd(cbind(1,MsplineBasis(1:xdims[4], df=control$df[4], degree=control$degree[4], intercept=FALSE)$X[,-control$df[4]]))$u)
          Ddf <- rep(control$df[4], nfac)
        }
      } else if(const[4]==6L){
        Dold <- Dnew <- replicate(nfac, rfunc(1:xdims[4], df=control$df[4], degree=control$degree[4], type="smo", nonneg=control$nonneg[4]))
        Ddf <- rep(0, nfac)
        if(!control$nonneg[4]) {
          Rd <- tcrossprod(svd(MsplineBasis(1:xdims[4], df=control$df[4], degree=control$degree[4], intercept=TRUE)$X)$u)
          Ddf <- rep(control$df[4], nfac)
        }
      }
      if(!is.null(Dstruc)) Dold <- Dnew <- Dold * Dstruc
    } else {
      Dold <- Dnew <- Dfixed
    }
    
    ### iterative update of matrices
    vtol <- sseold <- xcx
    iter <- 0
    cflag <- NA
    while(vtol>ctol && iter<maxit) {
      
      # Step 1: update mode A weights
      for(u in 1:nfac){DCkrB[,u] <- kronecker(Dold[,u],kronecker(Cold[,u],Bold[,u]))}
      if(const[1]==0L){
        #Anew <- Xa%*%DCkrB%*%smpower(crossprod(DCkrB),-1)
        Anew <- Xa%*%DCkrB%*%smpower(crossprod(Dold)*crossprod(Cold)*crossprod(Bold),-1)
      } else if(const[1]==1L) {
        Zmat <- Xa%*%DCkrB
        Anew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
      } else if(const[1]==2L) {
        cpmat <- crossprod(DCkrB)
        for(ii in 1:xdims[1]){Anew[ii,] <- fnnls(cpmat,crossprod(DCkrB,Xa[ii,]))}
        zix <- which(Anew <= 1e-9)
        if(length(zix) > 0) Anew[zix] <- 0
        if(any(colSums(Anew)==0)){
          Anew <- Aold
          vtol <- 0
          cflag <- 2
        }
      } else if(const[1]==3L) {
        for(u in 1:nfac){
          Zhat <- Xa - tcrossprod(Anew[,-u],DCkrB[,-u])
          fit <- fumls(1:xdims[1], Zhat %*% DCkrB[,u] / sum(DCkrB[,u]^2), df=control$df[1], degree=control$degree[1], lower=ifelse(control$nonneg[1],0,-Inf))
          Adf[u] <- fit$edf
          Anew[,u] <- fit$fitted.values
        }
        zix <- which(Anew <= 1e-9)
        if(length(zix) > 0) Anew[zix] <- 0
        if(any(colSums(Anew)==0)){
          Anew <- Aold
          vtol <- 0
          cflag <- 2
        }
      } else if(const[1]==4L){
        for(u in 1:nfac){
          Zhat <- Xa - tcrossprod(Anew[,-u],DCkrB[,-u])
          fit <- fmnls(1:xdims[1], Zhat %*% DCkrB[,u] / sum(DCkrB[,u]^2), df=control$df[1], degree=control$degree[1], lower=ifelse(control$nonneg[1],0,-Inf))
          Adf[u] <- fit$edf
          Anew[,u] <- fit$fitted.values
        }
        zix <- which(Anew <= 1e-9)
        if(length(zix) > 0) Anew[zix] <- 0
        if(any(colSums(Anew)==0)){
          Anew <- Aold
          vtol <- 0
          cflag <- 2
        }
      } else if(const[1]==5L | const[1]==6L){
        for(u in 1:nfac){
          Zhat <- Xa - tcrossprod(Anew[,-u],DCkrB[,-u])
          if(control$nonneg[1]){
            fit <- fsmls(1:xdims[1], Zhat %*% DCkrB[,u] / sum(DCkrB[,u]^2), df=control$df[1], degree=control$degree[1], nonneg=TRUE, periodic=ifelse(const[1]==5L,TRUE,FALSE))
            Adf[u] <- fit$edf
            Anew[,u] <- fit$fitted.values
          } else {
            Anew[,u] <- Ra %*% Zhat %*% DCkrB[,u] / sum(DCkrB[,u]^2)
          }
        }
        if(control$nonneg[1]){
          zix <- which(Anew <= 1e-9)
          if(length(zix) > 0) Anew[zix] <- 0
          if(any(colSums(Anew)==0)){
            Anew <- Aold
            vtol <- 0
            cflag <- 2
          }
        }
      } # end if(const[1]==0L)
      
      # Step 2: update mode B weights
      if(is.null(Bfixed)){
        for(u in 1:nfac){DCkrA[,u] <- kronecker(Dold[,u],kronecker(Cold[,u],Anew[,u]))}
        if(is.null(Bstruc)){
          if(const[2]==0L){
            #Bnew <- Xb%*%DCkrA%*%smpower(crossprod(DCkrA),-1)
            Bnew <- Xb%*%DCkrA%*%smpower(crossprod(Dold)*crossprod(Cold)*crossprod(Anew),-1)
          } else if(const[2]==1L) {
            Zmat <- Xb%*%DCkrA
            Bnew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
          } else if(const[2]==2L) {
            cpmat <- crossprod(DCkrA)
            for(ii in 1:xdims[2]){Bnew[ii,] <- fnnls(cpmat,crossprod(DCkrA,Xb[ii,]))}
            zix <- which(Bnew <= 1e-9)
            if(length(zix) > 0) Bnew[zix] <- 0
            if(any(colSums(Bnew)==0)){
              Bnew <- Bold
              vtol <- 0
              cflag <- 2
            }
          } else if(const[2]==3L){
            for(u in 1:nfac){
              Zhat <- Xb - tcrossprod(Bnew[,-u],DCkrA[,-u])
              fit <- fumls(1:xdims[2], Zhat %*% DCkrA[,u] / sum(DCkrA[,u]^2), df=control$df[2], degree=control$degree[2], lower=ifelse(control$nonneg[2],0,-Inf))
              Bdf[u] <- fit$edf
              Bnew[,u] <- fit$fitted.values
            }
            zix <- which(Bnew <= 1e-9)
            if(length(zix) > 0) Bnew[zix] <- 0
            if(any(colSums(Bnew)==0)){
              Bnew <- Bold
              vtol <- 0
              cflag <- 2
            }
          } else if(const[2]==4L){
            for(u in 1:nfac){
              Zhat <- Xb - tcrossprod(Bnew[,-u],DCkrA[,-u])
              fit <- fmnls(1:xdims[2], Zhat %*% DCkrA[,u] / sum(DCkrA[,u]^2), df=control$df[2], degree=control$degree[2], lower=ifelse(control$nonneg[2],0,-Inf))
              Bdf[u] <- fit$edf
              Bnew[,u] <- fit$fitted.values
            }
            zix <- which(Bnew <= 1e-9)
            if(length(zix) > 0) Bnew[zix] <- 0
            if(any(colSums(Bnew)==0)){
              Bnew <- Bold
              vtol <- 0
              cflag <- 2
            }
          } else if(const[2]==5L | const[2]==6L){
            for(u in 1:nfac){
              Zhat <- Xb - tcrossprod(Bnew[,-u],DCkrA[,-u])
              if(control$nonneg[2]){
                fit <- fsmls(1:xdims[2], Zhat %*% DCkrA[,u] / sum(DCkrA[,u]^2), df=control$df[2], degree=control$degree[2], nonneg=TRUE, periodic=ifelse(const[2]==5L,TRUE,FALSE))
                Bdf[u] <- fit$edf
                Bnew[,u] <- fit$fitted.values
              } else {
                Bnew[,u] <- Rb %*% Zhat %*% DCkrA[,u] / sum(DCkrA[,u]^2)
              }
            }
            if(control$nonneg[2]){
              zix <- which(Bnew <= 1e-9)
              if(length(zix) > 0) Bnew[zix] <- 0
              if(any(colSums(Bnew)==0)){
                Bnew <- Bold
                vtol <- 0
                cflag <- 2
              }
            }
          } # end if(const[2]==0L)
        } else {
          for(u in 1:nfac){
            Zhat <- Xb - tcrossprod(Bnew[,-u],DCkrA[,-u])
            Bnew[,u] <- ( (Zhat %*% DCkrA[,u]) / sum(DCkrA[,u]^2) ) * Bstruc[,u]
          }
        } # end if(is.null(Bstruc))
      } # end if(is.null(Bfixed))
      
      # Step 3: update mode C weights
      if(is.null(Cfixed)){
        for(u in 1:nfac){DBkrA[,u] <- kronecker(Dold[,u],kronecker(Bnew[,u],Anew[,u]))}
        if(is.null(Cstruc)){
          if(const[3]==0L){
            #Cnew <- Xc%*%DBkrA%*%smpower(crossprod(DBkrA),-1)
            Cnew <- Xc%*%DBkrA%*%smpower(crossprod(Dold)*crossprod(Bnew)*crossprod(Anew),-1)
          } else if(const[3]==1L) {
            Zmat <- Xc%*%DBkrA
            Cnew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
          } else if(const[3]==2L) {
            cpmat <- crossprod(DBkrA)
            for(ii in 1:xdims[3]){Cnew[ii,] <- fnnls(cpmat,crossprod(DBkrA,Xc[ii,]))}
            zix <- which(Cnew <= 1e-9)
            if(length(zix) > 0) Cnew[zix] <- 0
            if(any(colSums(Cnew)==0)){
              Cnew <- Cold
              vtol <- 0
              cflag <- 2
            }
          } else if(const[3]==3L){
            for(u in 1:nfac){
              Zhat <- Xc - tcrossprod(Cnew[,-u],DBkrA[,-u])
              fit <- fumls(1:xdims[3], Zhat %*% DBkrA[,u] / sum(DBkrA[,u]^2), df=control$df[3], degree=control$degree[3], lower=ifelse(control$nonneg[3],0,-Inf))
              Cdf[u] <- fit$edf
              Cnew[,u] <- fit$fitted.values
            }
            zix <- which(Cnew <= 1e-9)
            if(length(zix) > 0) Cnew[zix] <- 0
            if(any(colSums(Cnew)==0)){
              Cnew <- Cold
              vtol <- 0
              cflag <- 2
            }
          } else if(const[3]==4L){
            for(u in 1:nfac){
              Zhat <- Xc - tcrossprod(Cnew[,-u],DBkrA[,-u])
              fit <- fmnls(1:xdims[3], Zhat %*% DBkrA[,u] / sum(DBkrA[,u]^2), df=control$df[3], degree=control$degree[3], lower=ifelse(control$nonneg[3],0,-Inf))
              Cdf[u] <- fit$edf
              Cnew[,u] <- fit$fitted.values
            }
            zix <- which(Cnew <= 1e-9)
            if(length(zix) > 0) Cnew[zix] <- 0
            if(any(colSums(Cnew)==0)){
              Cnew <- Cold
              vtol <- 0
              cflag <- 2
            }
          } else if(const[3]==5L | const[3]==6L){
            for(u in 1:nfac){
              Zhat <- Xc - tcrossprod(Cnew[,-u],DBkrA[,-u])
              if(control$nonneg[3]){
                fit <- fsmls(1:xdims[3], Zhat %*% DBkrA[,u] / sum(DBkrA[,u]^2), df=control$df[3], degree=control$degree[3], nonneg=TRUE, periodic=ifelse(const[3]==5L,TRUE,FALSE))
                Cdf[u] <- fit$edf
                Cnew[,u] <- fit$fitted.values
              } else {
                Cnew[,u] <- Rc %*% Zhat %*% DBkrA[,u] / sum(DBkrA[,u]^2)
              }
            }
            if(control$nonneg[3]){
              zix <- which(Cnew <= 1e-9)
              if(length(zix) > 0) Cnew[zix] <- 0
              if(any(colSums(Cnew)==0)){
                Cnew <- Cold
                vtol <- 0
                cflag <- 2
              }
            }
          } # end if(const[3]==0L)
        } else {
          for(u in 1:nfac){
            Zhat <- Xc - tcrossprod(Cnew[,-u],DBkrA[,-u])
            Cnew[,u] <- ( (Zhat %*% DBkrA[,u]) / sum(DBkrA[,u]^2) ) * Cstruc[,u]
          }
        } # end if(is.null(Cstruc))
      } # end if(is.null(Cfixed))
      
      # Step 4: update mode D weights
      if(is.null(Dfixed)){
        for(u in 1:nfac){CBkrA[,u] <- kronecker(Cnew[,u],kronecker(Bnew[,u],Anew[,u]))}
        if(is.null(Dstruc)){
          if(const[4]==0L){
            #Dnew <- Xd%*%CBkrA%*%smpower(crossprod(CBkrA),-1)
            Dnew <- Xd%*%CBkrA%*%smpower(crossprod(Cnew)*crossprod(Bnew)*crossprod(Anew),-1)
          } else if(const[4]==1L) {
            Zmat <- Xd%*%CBkrA
            Dnew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
          } else if(const[4]==2L) {
            cpmat <- crossprod(CBkrA)
            for(ii in 1:xdims[4]){Dnew[ii,] <- fnnls(cpmat,crossprod(CBkrA,Xd[ii,]))}
            zix <- which(Dnew <= 1e-9)
            if(length(zix) > 0) Dnew[zix] <- 0
            if(any(colSums(Dnew)==0)){
              Dnew <- Dold
              vtol <- 0
              cflag <- 2
            }
          } else if(const[4]==3L){
            for(u in 1:nfac){
              Zhat <- Xd - tcrossprod(Dnew[,-u],CBkrA[,-u])
              fit <- fumls(1:xdims[4], Zhat %*% CBkrA[,u] / sum(CBkrA[,u]^2), df=control$df[4], degree=control$degree[4], lower=ifelse(control$nonneg[4],0,-Inf))
              Ddf[u] <- fit$edf
              Dnew[,u] <- fit$fitted.values
            }
            zix <- which(Dnew <= 1e-9)
            if(length(zix) > 0) Dnew[zix] <- 0
            if(any(colSums(Dnew)==0)){
              Dnew <- Dold
              vtol <- 0
              cflag <- 2
            }
          } else if(const[4]==4L){
            for(u in 1:nfac){
              Zhat <- Xd - tcrossprod(Dnew[,-u],CBkrA[,-u])
              fit <- fmnls(1:xdims[4], Zhat %*% CBkrA[,u] / sum(CBkrA[,u]^2), df=control$df[4], degree=control$degree[4], lower=ifelse(control$nonneg[4],0,-Inf))
              Ddf[u] <- fit$edf
              Dnew[,u] <- fit$fitted.values
            }
            zix <- which(Dnew <= 1e-9)
            if(length(zix) > 0) Dnew[zix] <- 0
            if(any(colSums(Dnew)==0)){
              Dnew <- Dold
              vtol <- 0
              cflag <- 2
            }
          } else if(const[4]==5L | const[4]==6L){
            for(u in 1:nfac){
              Zhat <- Xd - tcrossprod(Dnew[,-u],CBkrA[,-u])
              if(control$nonneg[4]){
                fit <- fsmls(1:xdims[4], Zhat %*% CBkrA[,u] / sum(CBkrA[,u]^2), df=control$df[4], degree=control$degree[4], nonneg=TRUE, periodic=ifelse(const[4]==5L,TRUE,FALSE))
                Ddf[u] <- fit$edf
                Dnew[,u] <- fit$fitted.values
              } else {
                Dnew[,u] <- Rd %*% Zhat %*% CBkrA[,u] / sum(CBkrA[,u]^2)
              }
            }
            if(control$nonneg[4]){
              zix <- which(Dnew <= 1e-9)
              if(length(zix) > 0) Dnew[zix] <- 0
              if(any(colSums(Dnew)==0)){
                Dnew <- Dold
                vtol <- 0
                cflag <- 2
              }
            }
          } # end if(const[3]==0L)
        } else {
          for(u in 1:nfac){
            Zhat <- Xd - tcrossprod(Dnew[,-u],CBkrA[,-u])
            Dnew[,u] <- ( (Zhat %*% CBkrA[,u]) / sum(CBkrA[,u]^2) ) * Dstruc[,u]
          }
        } # end if(is.null(Dstruc))
      } # end if(is.null(Dfixed))
      
      # Step 5: check for convergence
      for(u in 1:nfac){DCkrB[,u] <- kronecker(Dnew[,u],kronecker(Cnew[,u],Bnew[,u]))}
      ssenew <- sum((Xa-tcrossprod(Anew,DCkrB))^2)
      #vtol <- (sseold-ssenew)/sseold
      vtol <- (sseold - ssenew) / xcx
      Aold <- Anew
      Bold <- Bnew
      Cold <- Cnew
      Dold <- Dnew
      sseold <- ssenew
      iter <- iter + 1
      
    }  # end while(vtol>ctol && iter<maxit)
    
    ### scale and order solution
    if(is.null(Bfixed) & is.null(Cfixed) & is.null(Dfixed)){
      
      # put the scale in Mode D
      adg <- colMeans(Anew^2)
      Anew <- Anew%*%(diag(nfac)*(adg^-0.5))
      bdg <- colMeans(Bnew^2)
      Bnew <- Bnew%*%(diag(nfac)*(bdg^-0.5))
      cdg <- colMeans(Cnew^2)
      Cnew <- Cnew%*%(diag(nfac)*(cdg^-0.5))
      Dnew <- Dnew%*%(diag(nfac)*((adg*bdg*cdg)^0.5))
      
      # order according to sum-of-squares
      if(is.null(Bstruc) & is.null(Cstruc) & is.null(Dstruc)){
        fordr <- order(colSums(Dnew^2),decreasing=TRUE)
        Anew <- as.matrix(Anew[,fordr])
        Bnew <- as.matrix(Bnew[,fordr])
        Cnew <- as.matrix(Cnew[,fordr])
        Dnew <- as.matrix(Dnew[,fordr])
      }
      
    }
    
    
    # NEED TO FIX EDF!!
    
    ### effective degrees of freedom (Mode A)
    if(const[1]==0L){
      Adf <- nfac * (xdims[1] - 1L)
    } else if(const[1]==1L){
      Adf <- nfac * (xdims[1] - (nfac+1)/2)
    } else if(const[1]==2L){
      Adf <- nfac * (xdims[1] - 1L) - sum(Anew==0)
    } else if(const[1]>=3L & const[1]<=6L){
      Adf <- sum(Adf) - nfac
    }
    
    ### effective degrees of freedom (Mode B)
    if(is.null(Bfixed)){
      if(is.null(Bstruc)){
        if(const[2]==0L){
          Bdf <- nfac * (xdims[2] - 1L)
        } else if(const[2]==1L){
          Bdf <- nfac * (xdims[2] - (nfac+1)/2)
        } else if(const[2]==2L){
          Bdf <- nfac * (xdims[2] - 1L) - sum(Bnew==0)
        } else if(const[2]>=3L & const[2]<=6L){
          Bdf <- sum(Bdf) - nfac
        } 
      } else {
        Bdf <- sum(Bstruc) - nfac
      } # end if(is.null(Bstruc))
    } else {
      Bdf <- 0
    } # end if(is.null(Bfixed))
    
    ### effective degrees of freedom (Mode C)
    if(is.null(Cfixed)){
      if(is.null(Cstruc)){
        if(const[3]==0L){
          Cdf <- nfac * (xdims[3] - 1L)
        } else if(const[3]==1L){
          Cdf <- nfac * (xdims[3] - (nfac+1)/2)
        } else if(const[3]==2L){
          Cdf <- nfac * (xdims[3] - 1L) - sum(Cnew==0)
        } else if(const[3]>=3L & const[3]<=6L){
          Cdf <- sum(Cdf) - nfac
        } 
      } else {
        Cdf <- sum(Cstruc) - nfac
      } # end if(is.null(Cstruc))
    } else {
      Cdf <- 0
    } # end if(is.null(Cfixed))
    
    ### effective degrees of freedom (Mode D)
    if(is.null(Dfixed)){
      if(is.null(Dstruc)){
        if(const[4]==0L){
          Ddf <- nfac * xdims[4]
        } else if(const[4]==1L){
          Ddf <- nfac * (xdims[4] - (nfac-1)/2)
        } else if(const[4]==2L){
          Ddf <- nfac * xdims[4] - sum(Dnew==0)
        } else if(const[4]>=3L & const[4]<=6L){
          Ddf <- sum(Ddf)
        } 
      } else {
        Ddf <- sum(Dstruc)
      } # end if(is.null(Dstruc))
    } else {
      Ddf <- 0
      Adf <- Adf + nfac
    } # end if(is.null(Dfixed))
    
    ### GCV criterion
    edf <- c(Adf,Bdf,Cdf,Ddf)
    pxdim <- prod(xdims)
    GCV <- (ssenew/pxdim) / (1 - sum(edf)/pxdim)^2
    
    ### collect results
    Rsq <- 1 - ssenew/xcx
    if(is.na(cflag)){
      if(vtol<=ctol){cflag <- 0} else {cflag <- 1}
    }
    names(edf) <- c("A","B","C","D")
    fixed <- c(FALSE, ifelse(is.null(Bfixed), FALSE, TRUE), 
               ifelse(is.null(Cfixed), FALSE, TRUE), ifelse(is.null(Dfixed), FALSE, TRUE))
    struc <- c(FALSE, ifelse(is.null(Bstruc), FALSE, TRUE), 
               ifelse(is.null(Cstruc), FALSE, TRUE), ifelse(is.null(Dstruc), FALSE, TRUE))
    pfac <- list(A=Anew,B=Bnew,C=Cnew,D=Dnew,SSE=ssenew,Rsq=Rsq,GCV=GCV,edf=edf,
                 iter=iter,cflag=cflag,const=const,control=control,fixed=fixed,struc=struc)
    return(pfac)
    
  }