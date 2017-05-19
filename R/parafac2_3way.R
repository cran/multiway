parafac2_3way <- 
  function(data,nfac,xcx=sumsq(data),const=rep(0L,3),
           maxit=500,ctol=1e-4,Gfixed=NULL,Bfixed=NULL,
           Cfixed=NULL,Gstart=NULL,Bstart=NULL,Cstart=NULL,
           Gstruc=NULL,Bstruc=NULL,Cstruc=NULL,
           control=const.control(const)){
    # 3-way Parallel Factor Analysis 2 (Parafac2)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 16, 2017
    
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
    
    ### check for periodicity and smoothness on Mode A
    nx <- sapply(data, nrow)
    if(const[1]==5L){
      if(min(nx) == max(nx)){
        SAMEMODEA <- TRUE
        Ra <- tcrossprod(svd(cbind(1,MsplineBasis(1:nx[1], df=control$df[1], degree=control$degree[1], intercept=FALSE)$X[,-control$df[1]]))$u)
      } else {
        SAMEMODEA <- FALSE
        Ra <- vector("list",xdims[3])
        for(k in 1:xdims[3]) Ra[[k]] <- tcrossprod(svd(cbind(1,MsplineBasis(1:nx[k], df=control$df[1], degree=control$degree[1], intercept=FALSE)$X[,-control$df[1]]))$u)
      }
    } else if(const[1]==6L){
      if(min(nx) == max(nx)){
        SAMEMODEA <- TRUE
        Ra <- tcrossprod(svd(MsplineBasis(1:nx[1], df=control$df[1], degree=control$degree[1], intercept=TRUE)$X)$u)
      } else {
        SAMEMODEA <- FALSE
        Ra <- vector("list",xdims[3])
        for(k in 1:xdims[3]) Ra[[k]] <- tcrossprod(svd(MsplineBasis(1:nx[k], df=control$df[1], degree=control$degree[1], intercept=TRUE)$X)$u)
      }
    }
    
    ### initialize parameter matrices
    if(is.null(Gfixed)){
      if(!is.null(Gstart)){
        Gold <- Gnew <- Gstart
      } else if(const[1]==0L | const[1]==5L | const[1]==6L){
        GG <- matrix(rnorm(nfac^2),nfac,nfac)
        Gold <- Gnew <- GG %*% diag(sign(colSums(GG^3)))
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
      } else if (const[2]==6L){
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
        #Cold <- matrix(rnorm(xdims[3]*nfac),xdims[3],nfac)
        Cold <- matrix(runif(xdims[3]*nfac),xdims[3],nfac)
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
    
    ### iterative update of matrices
    vtol <- sseold <- xcx
    iter <- 0
    cflag <- NA
    while(vtol>ctol && iter<maxit) {
      
      ## Step 1: update mode A weights
      # 1a: update orthogonal projections
      if(any(const[1]==c(5L,6L))){
        if(SAMEMODEA){
          for(kk in 1:xdims[3]){
            xsvd <- svd(Ra%*%data[[kk]]%*%Bold%*%tcrossprod((diag(nfac)*Cold[kk,]),Gold))
            Rknew[[kk]] <- tcrossprod(xsvd$u,xsvd$v)
            Xtilde[,,kk] <- crossprod(Rknew[[kk]],data[[kk]])
          }
        } else {
          for(kk in 1:xdims[3]){
            xsvd <- svd(Ra[[kk]]%*%data[[kk]]%*%Bold%*%tcrossprod((diag(nfac)*Cold[kk,]),Gold))
            Rknew[[kk]] <- tcrossprod(xsvd$u,xsvd$v)
            Xtilde[,,kk] <- crossprod(Rknew[[kk]],data[[kk]])
          }
        }
      } else {
        for(kk in 1:xdims[3]){
          xsvd <- svd(data[[kk]]%*%Bold%*%tcrossprod((diag(nfac)*Cold[kk,]),Gold))
          Rknew[[kk]] <- tcrossprod(xsvd$u,xsvd$v)
          Xtilde[,,kk] <- crossprod(Rknew[[kk]],data[[kk]])
        }
      }
      # 1b: update correlation matrix
      if(is.null(Gfixed)){
        if(is.null(Gstruc)){
          if(const[1]!=1L){
            Xa <- matrix(Xtilde,nfac,xdims[2]*xdims[3])
            for(u in 1:nfac){CkrB[,u] <- kronecker(Cold[,u],Bold[,u])}
            #Gnew <- Xa%*%CkrB%*%smpower(crossprod(CkrB),-1)
            Gnew <- Xa%*%CkrB%*%smpower(crossprod(Cold)*crossprod(Bold),-1)
          }
        } else {
          Xa <- matrix(Xtilde,nfac,xdims[2]*xdims[3])
          for(u in 1:nfac){CkrB[,u] <- kronecker(Cold[,u],Bold[,u])}
          for(u in 1:nfac){
            Zhat <- Xa - tcrossprod(Gnew[,-u],CkrB[,-u])
            Gnew[,u] <- ( (Zhat %*% CkrB[,u]) / sum(CkrB[,u]^2) ) * Gstruc[,u]
          }
        } # end if(is.null(Gstruc))
      } # end if(is.null(Gfixed))
      
      ## Step 2: update mode B weights
      if(is.null(Bfixed)){
        Xb <- matrix(aperm(Xtilde,perm=c(2,1,3)),xdims[2],nfac*xdims[3])
        for(u in 1:nfac){CkrA[,u] <- kronecker(Cold[,u],Gnew[,u])}
        if(is.null(Bstruc)){
          if(const[2]==0L){
            #Bnew <- Xb%*%CkrA%*%smpower(crossprod(CkrA),-1)
            Bnew <- Xb%*%CkrA%*%smpower(crossprod(Cold)*crossprod(Gnew),-1)
          } else if(const[2]==1L) {
            Zmat <- Xb%*%CkrA
            Bnew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
          } else if(const[2]==2L) {
            cpmat <- crossprod(CkrA)
            for(ii in 1:xdims[2]){Bnew[ii,] <- fnnls(cpmat,crossprod(CkrA,Xb[ii,]))}
            zix <- which(Bnew <= 1e-9)
            if(length(zix) > 0) Bnew[zix] <- 0
            if(any(colSums(Bnew)==0)){
              Bnew <- Bold
              vtol <- 0
              cflag <- 2
            }
          } else if(const[2]==3L){
            for(u in 1:nfac){
              Zhat <- Xb - tcrossprod(Bnew[,-u],CkrA[,-u])
              fit <- fumls(1:xdims[2], Zhat %*% CkrA[,u] / sum(CkrA[,u]^2), df=control$df[2], degree=control$degree[2], lower=ifelse(control$nonneg[2],0,-Inf))
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
              Zhat <- Xb - tcrossprod(Bnew[,-u],CkrA[,-u])
              fit <- fmnls(1:xdims[2], Zhat %*% CkrA[,u] / sum(CkrA[,u]^2), df=control$df[2], degree=control$degree[2], lower=ifelse(control$nonneg[2],0,-Inf))
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
              Zhat <- Xb - tcrossprod(Bnew[,-u],CkrA[,-u])
              if(control$nonneg[2]){
                fit <- fsmls(1:xdims[2], Zhat %*% CkrA[,u] / sum(CkrA[,u]^2), df=control$df[2], degree=control$degree[2], nonneg=TRUE, periodic=ifelse(const[2]==5L,TRUE,FALSE))
                Bdf[u] <- fit$edf
                Bnew[,u] <- fit$fitted.values
              } else {
                Bnew[,u] <- Rb %*% Zhat %*% CkrA[,u] / sum(CkrA[,u]^2)
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
            Zhat <- Xb - tcrossprod(Bnew[,-u],CkrA[,-u])
            Bnew[,u] <- ( (Zhat %*% CkrA[,u]) / sum(CkrA[,u]^2) ) * Bstruc[,u]
          }
        } # end if(is.null(Bstruc))
      } # end if(is.null(Bfixed))
      
      ## Step 3: update mode C weights
      if(is.null(Cfixed)){
        Xc <- matrix(aperm(Xtilde,perm=c(3,1,2)),xdims[3],nfac*xdims[2])
        for(u in 1:nfac){BkrA[,u] <- kronecker(Bnew[,u],Gnew[,u])}
        if(is.null(Cstruc)){
          if(const[3]==0L){
            #Cnew <- Xc%*%BkrA%*%smpower(crossprod(BkrA),-1)
            Cnew <- Xc%*%BkrA%*%smpower(crossprod(Bnew)*crossprod(Gnew),-1)
          } else if(const[3]==1L) {
            Zmat <- Xc%*%BkrA
            Cnew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
          } else if(const[3]==2L) {
            cpmat <- crossprod(BkrA)
            for(ii in 1:xdims[3]){Cnew[ii,] <- fnnls(cpmat,crossprod(BkrA,Xc[ii,]))}
            zix <- which(Cnew <= 1e-9)
            if(length(zix) > 0) Cnew[zix] <- 0
            if(any(colSums(Cnew)==0)){
              Cnew <- Cold
              vtol <- 0
              cflag <- 2
            }
          } else if(const[3]==3L){
            for(u in 1:nfac){
              Zhat <- Xc - tcrossprod(Cnew[,-u],BkrA[,-u])
              fit <- fumls(1:xdims[3], Zhat %*% BkrA[,u] / sum(BkrA[,u]^2), df=control$df[3], degree=control$degree[3], lower=ifelse(control$nonneg[3],0,-Inf))
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
              Zhat <- Xc - tcrossprod(Cnew[,-u],BkrA[,-u])
              fit <- fmnls(1:xdims[3], Zhat %*% BkrA[,u] / sum(BkrA[,u]^2), df=control$df[3], degree=control$degree[3], lower=ifelse(control$nonneg[3],0,-Inf))
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
              Zhat <- Xc - tcrossprod(Cnew[,-u],BkrA[,-u])
              if(control$nonneg[3]){
                fit <- fsmls(1:xdims[3], Zhat %*% BkrA[,u] / sum(BkrA[,u]^2), df=control$df[3], degree=control$degree[3], nonneg=TRUE, periodic=ifelse(const[3]==5L,TRUE,FALSE))
                Cdf[u] <- fit$edf
                Cnew[,u] <- fit$fitted.values
              } else {
                Cnew[,u] <- Rc %*% Zhat %*% BkrA[,u] / sum(BkrA[,u]^2)
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
            Zhat <- Xc - tcrossprod(Cnew[,-u],BkrA[,-u])
            Cnew[,u] <- ( (Zhat %*% BkrA[,u]) / sum(BkrA[,u]^2) ) * Cstruc[,u]
          }
        } # end if(is.null(Cstruc))
      } # end if(is.null(Cfixed))
      
      ## Step 4: check for convergence
      ssenew <- 0
      for(kk in 1:xdims[3]){
        ssenew <- ssenew + sum((data[[kk]]-tcrossprod(Rknew[[kk]]%*%Gnew%*%(diag(nfac)*Cnew[kk,]),Bnew))^2)
      }
      #vtol <- (sseold-ssenew)/sseold
      vtol <- (sseold - ssenew) / xcx
      Gold <- Gnew
      Bold <- Bnew
      Cold <- Cnew
      sseold <- ssenew
      iter <- iter + 1
      
    }  # end while(vtol>ctol && iter<maxit)
    
    ### scale and order solution
    if(is.null(Gfixed) & is.null(Bfixed) & is.null(Cfixed)){
      
      # put the scale in Mode C
      adg <- colSums(Gnew^2)
      Gnew <- Gnew%*%(diag(nfac)*(adg^-0.5))
      bdg <- colMeans(Bnew^2)
      Bnew <- Bnew%*%(diag(nfac)*(bdg^-0.5))
      Cnew <- Cnew%*%(diag(nfac)*((adg*bdg)^0.5))
      
      # order according to sum-of-squares
      if(is.null(Gstruc) & is.null(Bstruc) & is.null(Cstruc)){
        fordr <- order(colSums(Cnew^2),decreasing=TRUE)
        Gnew <- as.matrix(Gnew[,fordr])
        Bnew <- as.matrix(Bnew[,fordr])
        Cnew <- as.matrix(Cnew[,fordr])
      }
      
    }
    
    ### effective degrees of freedom (Mode A; orthogonal basis)
    if(any(const[1]==c(5L,6L))){
      ntotal <- xdims[3] * control$df[1]
    } else {
      ntotal <- sum(nx)
    } # end if(const[1]==5L)
    Adf <- nfac * (ntotal - xdims[3]*(nfac+1)/2)
    
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
          Cdf <- nfac * xdims[3]
        } else if(const[3]==1L){
          Cdf <- nfac * (xdims[3] - (nfac-1)/2)
        } else if(const[3]==2L){
          Cdf <- nfac * xdims[3] - sum(Cnew==0)
        } else if(const[3]>=3L & const[3]<=6L){
          Cdf <- sum(Cdf)
        } 
      } else {
        Cdf <- sum(Cstruc)
      } # end if(is.null(Cstruc))
    } else {
      Cdf <- 0
      if(const[1]==1L) {
        Bdf <- Bdf + nfac
      } else {
        Gdf <- Gdf + nfac
      }
    } # end if(is.null(Cfixed))
    
    ### GCV criterion
    edf <- c(Adf+Gdf,Bdf,Cdf)
    pxdim <- sum(nx) * xdims[2]
    GCV <- (ssenew/pxdim) / (1 - sum(edf)/pxdim)^2
    
    ### collect results
    Rsq <- 1 - ssenew/xcx
    if(is.na(cflag)){
      if(vtol<=ctol){cflag <- 0} else {cflag <- 1}
    }
    names(edf) <- c("A","B","C")
    fixed <- c(ifelse(is.null(Gfixed), FALSE, TRUE), ifelse(is.null(Bfixed), FALSE, TRUE), ifelse(is.null(Cfixed), FALSE, TRUE))
    struc <- c(ifelse(is.null(Gstruc), FALSE, TRUE), ifelse(is.null(Bstruc), FALSE, TRUE), ifelse(is.null(Cstruc), FALSE, TRUE))
    Ak <- vector("list", xdims[3])
    for(k in 1:xdims[3]) Ak[[k]] <- Rknew[[k]] %*% Gnew
    pfac <- list(A=Ak,B=Bnew,C=Cnew,Phi=crossprod(Gnew),SSE=ssenew,Rsq=Rsq,
                 GCV=GCV,edf=edf,iter=iter,cflag=cflag,const=const,control=control,
                 fixed=fixed,struc=struc)
    return(pfac)
    
  }