parafac2_4wayna <- 
  function(data,nfac,naid=NULL,const=rep(0L,4),
           maxit=500,ctol=1e-4,Gfixed=NULL,Bfixed=NULL,
           Cfixed=NULL,Dfixed=NULL,Gstart=NULL,Bstart=NULL,
           Cstart=NULL,Dstart=NULL,Gstruc=NULL,Bstruc=NULL,
           Cstruc=NULL,Dstruc=NULL,control=const.control(const)){
    # 4-way Parallel Factor Analysis 2 (Parafac2)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: November 9, 2017
    
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
    
    ### initialize missing data
    if(is.null(naid)) naid <- lapply(data, function(x) which(is.na(x)))
    nmiss <- sapply(naid, length)
    for(k in 1:xdims[4]){
      if(nmiss[k] > 0) data[[k]][naid[[k]]] <- rnorm(nmiss[k])
    }
    xcx <- sumsq(data)
    Xhat <- vector("list", xdims[4])
    
    ### initialize stuff for Mode A update
    Rknew <- vector("list",xdims[4])
    Xtilde <- array(0,dim=c(nfac,xdims[2],xdims[3],xdims[4]))
    
    ### reshape raw data
    for(kk in 1:xdims[4]){
      mdim <- dim(data[[kk]])
      data[[kk]] <- matrix(data[[kk]],mdim[1],xdims[2]*xdims[3])
    }
    
    ### check for smoothness on Mode A
    nx <- sapply(data, nrow)
    if(const[1]==5L){
      if(min(nx) == max(nx)){
        SAMEMODEA <- TRUE
        Ra <- tcrossprod(svd(cbind(1,MsplineBasis(1:nx[1], df=control$df[1], degree=control$degree[1], intercept=FALSE)$X[,-control$df[1]]))$u)
      } else {
        SAMEMODEA <- FALSE
        Ra <- vector("list",xdims[4])
        for(k in 1:xdims[4]) Ra[[k]] <- tcrossprod(svd(cbind(1,MsplineBasis(1:nx[k], df=control$df[1], degree=control$degree[1], intercept=FALSE)$X[,-control$df[1]]))$u)
      }
    }
    if(const[1]==6L){
      if(min(nx) == max(nx)){
        SAMEMODEA <- TRUE
        Ra <- tcrossprod(svd(MsplineBasis(1:nx[1], df=control$df[1], degree=control$degree[1], intercept=TRUE)$X)$u)
      } else {
        SAMEMODEA <- FALSE
        Ra <- vector("list",xdims[4])
        for(k in 1:xdims[4]) Ra[[k]] <- tcrossprod(svd(MsplineBasis(1:nx[k], df=control$df[1], degree=control$degree[1], intercept=TRUE)$X)$u)
      }
    }
    
    ### initialize parameter matrices
    if(is.null(Gfixed)){
      if(!is.null(Gstart)){
        Gold <- Gnew <- Gstart
      } else if(const[1]==0L | const[1]==5L | const[1]==6L){
        Gold <- Gnew <- matrix(rnorm(nfac^2),nfac,nfac)
        #Gold <- Gold %*% (diag(nfac)/sqrt(colSums(Gold^2)))
      } else if(const[1]==1L){
        Gold <- Gnew <- diag(nfac)
      }
      if(!is.null(Gstruc)) Gold <- Gnew <- Gold * Gstruc
    } else {
      Gold <- Gnew <- Gfixed
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
      
      ## Step 1: update mode A weights
      # 1a: update orthogonal projections
      for(u in 1:nfac){CkrB[,u] <- kronecker(Cold[,u],Bold[,u])}
      if(any(const[1]==c(5L,6L))){
        if(SAMEMODEA){
          for(kk in 1:xdims[4]){
            xsvd <- svd(Ra%*%data[[kk]]%*%CkrB%*%tcrossprod((diag(nfac)*Dold[kk,]),Gold))
            Rknew[[kk]] <- tcrossprod(xsvd$u,xsvd$v)
            Xtilde[,,,kk] <- array(crossprod(Rknew[[kk]],data[[kk]]),dim=c(nfac,xdims[2],xdims[3]))
          }
        } else {
          for(kk in 1:xdims[4]){
            xsvd <- svd(Ra[[kk]]%*%data[[kk]]%*%CkrB%*%tcrossprod((diag(nfac)*Dold[kk,]),Gold))
            Rknew[[kk]] <- tcrossprod(xsvd$u,xsvd$v)
            Xtilde[,,,kk] <- array(crossprod(Rknew[[kk]],data[[kk]]),dim=c(nfac,xdims[2],xdims[3]))
          }
        }
      } else {
        for(kk in 1:xdims[4]){
          xsvd <- svd(data[[kk]]%*%CkrB%*%tcrossprod((diag(nfac)*Dold[kk,]),Gold))
          Rknew[[kk]] <- tcrossprod(xsvd$u,xsvd$v)
          Xtilde[,,,kk] <- array(crossprod(Rknew[[kk]],data[[kk]]),dim=c(nfac,xdims[2],xdims[3]))
        }
      }
      
      # 1b: update correlation matrix
      if(is.null(Gfixed)){
        if(is.null(Gstruc)){
          if(const[1]!=1L){
            Xa <- matrix(Xtilde,nfac,xdims[2]*xdims[3]*xdims[4])
            for(u in 1:nfac){DCkrB[,u] <- kronecker(Dold[,u],kronecker(Cold[,u],Bold[,u]))}
            #Gnew <- Xa%*%DCkrB%*%smpower(crossprod(DCkrB),-1)
            Gnew <- Xa%*%DCkrB%*%smpower(crossprod(Dold)*crossprod(Cold)*crossprod(Bold),-1)
          }
        } else {
          Xa <- matrix(Xtilde,nfac,xdims[2]*xdims[3]*xdims[4])
          for(u in 1:nfac){DCkrB[,u] <- kronecker(Dold[,u],kronecker(Cold[,u],Bold[,u]))}
          for(u in 1:nfac){
            Zhat <- Xa - tcrossprod(Gnew[,-u],DCkrB[,-u])
            Gnew[,u] <- ( (Zhat %*% DCkrB[,u]) / sum(DCkrB[,u]^2) ) * Gstruc[,u]
          }
        } # end if(is.null(Gstruc))
      } # end if(is.null(Gfixed))
      
      ## Step 2: update mode B weights
      if(is.null(Bfixed)){
        Xb <- matrix(aperm(Xtilde,perm=c(2,1,3,4)),xdims[2],nfac*xdims[3]*xdims[4])
        for(u in 1:nfac){DCkrA[,u] <- kronecker(Dold[,u],kronecker(Cold[,u],Gnew[,u]))}
        if(is.null(Bstruc)){
          if(const[2]==0L){
            #Bnew <- Xb%*%DCkrA%*%smpower(crossprod(DCkrA),-1)
            Bnew <- Xb%*%DCkrA%*%smpower(crossprod(Dold)*crossprod(Cold)*crossprod(Gnew),-1)
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
      
      ## Step 3: update mode C weights
      if(is.null(Cfixed)){
        Xc <- matrix(aperm(Xtilde,perm=c(3,1,2,4)),xdims[3],nfac*xdims[2]*xdims[4])
        for(u in 1:nfac){DBkrA[,u] <- kronecker(Dold[,u],kronecker(Bnew[,u],Gnew[,u]))}
        if(is.null(Cstruc)){
          if(const[3]==0L){
            #Cnew <- Xc%*%DBkrA%*%smpower(crossprod(DBkrA),-1)
            Cnew <- Xc%*%DBkrA%*%smpower(crossprod(Dold)*crossprod(Bnew)*crossprod(Gnew),-1)
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
      
      ## Step 4: update mode D weights
      if(is.null(Dfixed)){
        Xd <- matrix(aperm(Xtilde,perm=c(4,1,2,3)),xdims[4],nfac*xdims[2]*xdims[3])
        for(u in 1:nfac){CBkrA[,u] <- kronecker(Cnew[,u],kronecker(Bnew[,u],Gnew[,u]))}
        if(is.null(Dstruc)){
          if(const[4]==0L){
            #Dnew <- Xd%*%CBkrA%*%smpower(crossprod(CBkrA),-1)
            Dnew <- Xd%*%CBkrA%*%smpower(crossprod(Cnew)*crossprod(Bnew)*crossprod(Gnew),-1)
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
          } # end if(const[3]==0L)
        } else {
          for(u in 1:nfac){
            Zhat <- Xd - tcrossprod(Dnew[,-u],CBkrA[,-u])
            Dnew[,u] <- ( (Zhat %*% CBkrA[,u]) / sum(CBkrA[,u]^2) ) * Dstruc[,u]
          }
        } # end if(is.null(Dstruc))
      } # end if(is.null(Dfixed))
      
      ## Step 5: check for convergence
      for(u in 1:nfac){CkrB[,u] <- kronecker(Cnew[,u],Bnew[,u])}
      ssenew <- 0
      for(kk in 1:xdims[4]){
        Xhat[[kk]] <- tcrossprod(Rknew[[kk]]%*%Gnew%*%(diag(nfac)*Dnew[kk,]),CkrB)
        ssenew <- ssenew + sum((data[[kk]]-Xhat[[kk]])^2)
      }
      #vtol <- (sseold-ssenew)/sseold
      vtol <- (sseold - ssenew) / xcx
      Gold <- Gnew
      Bold <- Bnew
      Cold <- Cnew
      Dold <- Dnew
      sseold <- ssenew
      iter <- iter + 1
      
      # impute missing data
      for(k in 1:xdims[4]){
        if(nmiss[k] > 0) data[[k]][naid[[k]]] <- Xhat[[k]][naid[[k]]]
      }
      xcx <- sumsq(data)
      
    }  # end while(vtol>ctol && iter<maxit)
    
    ### update SSE
    ssenew <- 0
    for(kk in 1:xdims[4]) ssenew <- ssenew + sum((data[[kk]]-Xhat[[kk]])^2)
    
    ### scale and order solution
    if(is.null(Gfixed) & is.null(Bfixed) & is.null(Cfixed) & is.null(Dfixed)){
      
      # put the scale in Mode D
      adg <- colSums(Gnew^2)
      Gnew <- Gnew%*%(diag(nfac)*(adg^-0.5))
      bdg <- colMeans(Bnew^2)
      Bnew <- Bnew%*%(diag(nfac)*(bdg^-0.5))
      cdg <- colMeans(Cnew^2)
      Cnew <- Cnew%*%(diag(nfac)*(cdg^-0.5))
      Dnew <- Dnew%*%(diag(nfac)*((adg*bdg*cdg)^0.5))
      
      # order according to sum-of-squares
      if(is.null(Gstruc) & is.null(Bstruc) & is.null(Cstruc) & is.null(Dstruc)){
        fordr <- order(colSums(Dnew^2),decreasing=TRUE)
        Gnew <- as.matrix(Gnew[,fordr])
        Bnew <- as.matrix(Bnew[,fordr])
        Cnew <- as.matrix(Cnew[,fordr])
        Dnew <- as.matrix(Dnew[,fordr])
      }
      
    }
    
    ### effective degrees of freedom (Mode A; orthogonal basis)
    if(any(const[1]==c(5L,6L))){
      ntotal <- xdims[3] * control$df[1]
    } else {
      ntotal <- sum(nx)
    } # end if(any(const[1]==c(5L,6L)))
    Adf <- nfac * (ntotal - xdims[4]*(nfac+1)/2)
    
    ### effective degrees of freedom (Mode A; crossproducts)
    if(is.null(Gfixed)){
      if(is.null(Gstruc)){
        if(const[1]==1L){
          Gdf <- 0
        } else {
          Gdf <- nfac * (nfac - 1L) / 2
        }
      } else {
        GtG <- crossprod(Gstruc)
        Gdf <- sum(GtG[lower.tri(GtG)]>0L)
      } # end if(is.null(Gstruc))
    } else {
      Gdf <- 0
    } # end if(is.null(Gfixed))
    
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
      } # end if(is.null(Cstruc))
    } else {
      Ddf <- 0
      if(const[1]==1L){
        Bdf <- Bdf + nfac
      } else {
        Gdf <- Gdf + nfac
      }
    } # end if(is.null(Cfixed))
    
    ### GCV criterion
    edf <- c(Adf+Gdf,Bdf,Cdf,Ddf)
    pxdim <- sum(nx) * prod(xdims[2:3])
    GCV <- (ssenew/pxdim) / (1 - sum(edf)/pxdim)^2
    
    ### collect results
    Rsq <- 1 - ssenew/xcx
    if(is.na(cflag)){
      if(vtol<=ctol){cflag <- 0} else {cflag <- 1}
    }
    names(edf) <- c("A","B","C","D")
    fixed <- c(ifelse(is.null(Gfixed), FALSE, TRUE), ifelse(is.null(Bfixed), FALSE, TRUE), 
               ifelse(is.null(Cfixed), FALSE, TRUE), ifelse(is.null(Dfixed), FALSE, TRUE))
    struc <- c(ifelse(is.null(Gstruc), FALSE, TRUE), ifelse(is.null(Bstruc), FALSE, TRUE), 
               ifelse(is.null(Cstruc), FALSE, TRUE), ifelse(is.null(Dstruc), FALSE, TRUE))
    Ak <- vector("list", xdims[4])
    for(k in 1:xdims[4]) Ak[[k]] <- Rknew[[k]] %*% Gnew
    pfac <- list(A=Ak,B=Bnew,C=Cnew,D=Dnew,Phi=crossprod(Gnew),SSE=ssenew,Rsq=Rsq,
                 GCV=GCV,edf=edf,iter=iter,cflag=cflag,const=const,control=control,
                 fixed=fixed,struc=struc)
    return(pfac)
    
  }