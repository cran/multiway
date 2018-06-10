mcr_tucker <- 
  function(dataX, dataY, nfac, alpha = 0.5, const = NULL,
           ssx = NULL, ssy = NULL, maxit = 500, ctol = 1e-4,
           Afixed = NULL, Bfixed = NULL, Cfixed = NULL, Dfixed = NULL,
           Astart = NULL, Bstart = NULL, Cstart = NULL, Dstart = NULL,
           projY = NULL, mode = 3){
    # Multi-way Covariates Regression (Tucker 3-way)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: Feb 1, 2018
    
    ### get dimensions
    xdims <- dim(dataX)
    ydims <- dim(dataY)
    
    ### check ssx and ssy
    if(is.null(ssx)) ssx <- sumsq(dataX)
    if(is.null(ssy)) ssy <- sumsq(dataY)
    
    ### initialize reshaped data matrices
    Xa <- matrix(dataX, nrow = xdims[1], ncol = xdims[2] * xdims[3])
    Xb <- matrix(aperm(dataX, perm = c(2,1,3)), nrow = xdims[2], ncol = xdims[1] * xdims[3])
    Xc <- matrix(aperm(dataX, perm = c(3,1,2)), nrow = xdims[3], ncol = xdims[1] * xdims[2])
    rm(dataX)
    
    ### project response
    if(xdims[mode] > prod(xdims[-mode])){
      if(is.null(projY)){
        if(mode == 1L){
          xsvd <- svd(Xa, nv = 0)
        } else if(mode == 2L){
          xsvd <- svd(Xb, nv = 0)
        } else if(mode == 3L){
          xsvd <- svd(Xc, nv = 0)
        }
        projY <- xsvd$u %*% crossprod(xsvd$u, dataY)
      }
    } else {
      projY <- dataY
    }
    
    ### initialize matrices
    Amat <- Afixed
    Bmat <- Bfixed
    if(is.null(Bmat)) Bmat <- Bstart
    if(is.null(Bmat)){
      Bmat <- svd(matrix(rnorm(xdims[2]*nfac[2]), nrow = xdims[2], ncol = nfac[2]), nu = nfac[2], nv = 0)$u
    }
    Cmat <- Cfixed
    if(is.null(Cmat)) Cmat <- Cstart
    if(is.null(Cmat)){
      Cmat <- svd(matrix(rnorm(xdims[3]*nfac[3]), nrow = xdims[3], ncol = nfac[3]), nu = nfac[3], nv = 0)$u
    }
    Dmat <- Dfixed
    if(is.null(Dmat)) Dmat <- Dstart
    if(is.null(Dmat)){
      Dmat <- matrix(rnorm(ydims[2]*nfac[mode]), nrow = ydims[2], ncol = nfac[mode])
    }
    Gmat <- array(rnorm(prod(nfac)), dim = nfac)
    
    ### iterative update of matrices
    wts <- c(alpha / ssx, (1-alpha) / ssy)
    ssxssy <- 1 + ctol
    vtol <- mseold <- ssxssy
    iter <- 0
    cflag <- NA
    while((vtol > ctol) && (iter < maxit)) {
      
      # Step 1: update mode A weights
      if(is.null(Afixed)){
        CkrB <- kronecker(Cmat, Bmat)
        if(mode == 1L){
          Ga <- matrix(Gmat, nrow = nfac[1], ncol = prod(nfac[-1]))
          Qmat <- wts[1] * Xa %*% CkrB %*% t(Ga) + wts[2] * projY %*% Dmat
          Qsvd <- svd(Qmat, nu = nfac[1], nv = nfac[1])
          Amat <- tcrossprod(Qsvd$u, Qsvd$v)
        } else {
          Amat <- svd(Xa %*% CkrB, nu=nfac[1], nv=0)$u
        }
      }
      
      # Step 2: update mode B weights
      if(is.null(Bfixed)){
        CkrA <- kronecker(Cmat,Amat)
        if(mode == 2L){
          Gb <- matrix(aperm(Gmat, perm = c(2,1,3)), nrow = nfac[2], ncol = prod(nfac[-2]))
          Qmat <- wts[1] * Xb %*% CkrA %*% t(Gb) + wts[2] * projY %*% Dmat
          Qsvd <- svd(Qmat, nu = nfac[2], nv = nfac[2])
          Bmat <- tcrossprod(Qsvd$u, Qsvd$v)
        } else {
          Bmat <- svd(Xb %*% CkrA, nu=nfac[2], nv=0)$u
        }
      }
      
      # Step 3: update mode C weights
      if(is.null(Cfixed)){
        BkrA <- kronecker(Bmat,Amat)
        if(mode == 3L){
          Gc <- matrix(aperm(Gmat, perm = c(3,1,2)), nrow = nfac[3], ncol = prod(nfac[-3]))
          Qmat <- wts[1] * Xc %*% BkrA %*% t(Gc) + wts[2] * projY %*% Dmat
          Qsvd <- svd(Qmat, nu = nfac[3], nv = nfac[3])
          Cmat <- tcrossprod(Qsvd$u, Qsvd$v)
        } else {
          Cmat <- svd(Xc %*% BkrA, nu=nfac[3], nv=0)$u
        }
      }
      
      # Step 4: update regression coefficients
      if(is.null(Dfixed)){
        if(mode == 1L){
          Dmat <- crossprod(dataY, Amat)
        } else if(mode == 2L){
          Dmat <- crossprod(dataY, Bmat)
        } else if(mode == 3L){
          Dmat <- crossprod(dataY, Cmat)
        }
      } # end if(is.null(Dfixed))
      if(mode == 1L){
        ssenew2 <- sum((dataY - tcrossprod(Amat, Dmat))^2)
      } else if(mode == 2L){
        ssenew2 <- sum((dataY - tcrossprod(Bmat, Dmat))^2)
      } else {
        ssenew2 <- sum((dataY - tcrossprod(Cmat, Dmat))^2)
      }
      
      # Step 5: check for convergence
      Ga <- crossprod(Amat, Xa %*% kronecker(Cmat, Bmat))
      Gmat <- array(Ga, dim = nfac)
      ssenew1 <- ssx - sum(Ga^2)
      msenew <- wts[1] * ssenew1 + wts[2] * ssenew2
      vtol <- (mseold - msenew) / (mseold + ctol)
      mseold <- msenew
      iter <- iter + 1
      
    } # end while(vtol>ctol && iter<maxit)
    
    ### get coefficients
    if(mode == 1L){
      Wmat <- mpinv(Xa) %*% Amat
    } else if(mode == 2L){
      Wmat <- mpinv(Xb) %*% Bmat
    } else if(mode == 3L){
      Wmat <- mpinv(Xc) %*% Cmat
    }
    
    ### collect results
    SSE <- c(X = ssenew1, Y = ssenew2)
    RsqX <- 1 - ssenew1 / ssx
    RsqY <- 1 - ssenew2 / ssy
    Rsq <- c(X = RsqX, Y = RsqY)
    cflag <- ifelse(vtol <= ctol, 0, 1)
    fixed <- c(ifelse(is.null(Afixed), FALSE, TRUE), ifelse(is.null(Bfixed), FALSE, TRUE), 
               ifelse(is.null(Cfixed), FALSE, TRUE), ifelse(is.null(Dfixed), FALSE, TRUE))
    struc <- rep(FALSE, 4L)
    tfac <- list(A = Amat, B = Bmat, C = Cmat, D = Dmat, W = Wmat,
                 LOSS = msenew, SSE = SSE, Rsq = Rsq, iter = iter,
                 cflag = cflag, model = "tucker", const = rep("uncons", 4L), 
                 control = NULL, weights = NULL, alpha = alpha, 
                 fixed = fixed, struc = struc, G = Gmat)
    return(tfac)
    
  } # end mcr_tucker