mcr_parafac <- 
  function(dataX, dataY, nfac, alpha = 0.5, ssx = NULL, ssy = NULL, 
           const = rep("uncons", 4L), control = const.control(const), 
           Afixed = NULL, Bfixed = NULL, Cfixed = NULL, Dfixed = NULL,
           Astart = NULL, Bstart = NULL, Cstart = NULL, Dstart = NULL,
           Astruc = NULL, Bstruc = NULL, Cstruc = NULL, Dstruc = NULL,
           Amodes = NULL, Bmodes = NULL, Cmodes = NULL, Dmodes = NULL,
           maxit = 500, ctol = 1e-4, projY = NULL, mode = 3, backfit = FALSE){
    # Multi-way Covariates Regression (Parafac 3-way)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: April 19, 2018
    
    ### get dimensions
    xdims <- dim(dataX)
    ydims <- dim(dataY)
    
    ### check ssx and ssy
    if(is.null(ssx)) ssx <- sumsq(dataX)
    if(is.null(ssy)) ssy <- sumsq(dataY)
    
    ### initialize Khatri-Rao product matrices
    BkrA <- matrix(0, nrow = xdims[1] * xdims[2], ncol = nfac)
    CkrA <- matrix(0, nrow = xdims[1] * xdims[3], ncol = nfac)
    CkrB <- matrix(0, nrow = xdims[2] * xdims[3], ncol = nfac)
    
    ### initialize reshaped data matrices
    Xa <- matrix(dataX, nrow = xdims[1], ncol = xdims[2]*xdims[3])
    Xb <- matrix(aperm(dataX, perm = c(2,1,3)), nrow = xdims[2], ncol = xdims[1] * xdims[3])
    Xc <- matrix(aperm(dataX, perm = c(3,1,2)), nrow = xdims[3], ncol = xdims[1] * xdims[2])
    rm(dataX)
    
    ### check mode
    mode <- as.integer(mode[1])
    if(!any(mode == c(1L, 2L, 3L))) mode <- 3L
    
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
      Bmat <- initcmls(nobs = xdims[2], nfac = nfac, const = const[2],
                       struc = Bstruc, df = control$df[2], 
                       degree = control$degree[2], 
                       intercept = control$intercept[2],
                       mode.range = Bmodes)
    }
    Cmat <- Cfixed
    if(is.null(Cmat)) Cmat <- Cstart
    if(is.null(Cmat)){
      Cmat <- initcmls(nobs = xdims[3], nfac = nfac, const = const[3],
                       struc = Cstruc, df = control$df[3], 
                       degree = control$degree[3], 
                       intercept = control$intercept[3],
                       mode.range = Cmodes)
    }
    Dmat <- Dfixed
    if(is.null(Dmat)) Dmat <- Dstart
    if(is.null(Dmat)){
      Dmat <- initcmls(nobs = ydims[2], nfac = nfac, const = const[4],
                       struc = Dstruc, df = control$df[4], 
                       degree = control$degree[4], 
                       intercept = control$intercept[4],
                       mode.range = Dmodes)
    }
    
    ### transpose struc arguments (if given)
    if(!is.null(Astruc)) Astruc <- t(Astruc)
    if(!is.null(Bstruc)) Bstruc <- t(Bstruc)
    if(!is.null(Cstruc)) Cstruc <- t(Cstruc)
    if(!is.null(Dstruc)) Dstruc <- t(Dstruc)
    
    ### iterative update of matrices
    wts <- c(alpha / ssx, (1-alpha) / ssy)
    ssxssy <- 1 + ctol
    xtol <- c(xdims, ydims[2]) * .Machine$double.eps
    vtol <- mseold <- ssxssy
    iter <- 0
    cflag <- NA
    while((vtol > ctol) && (iter < maxit)) {
      
      # Step 1: update mode A weights
      if(is.null(Afixed) && is.na(cflag)){
        
        for(u in 1:nfac) CkrB[,u] <- kronecker(Cmat[,u], Bmat[,u])
        if(mode == 1L){
          Qmat <- wts[1] * Xa %*% CkrB + wts[2] * projY %*% Dmat
          if(const[1] == "orthog"){
            Qsvd <- svd(Qmat)
            Amat <- tcrossprod(Qsvd$u, Qsvd$v)
          } else {
            Tmat <- wts[1] * crossprod(Bmat) * crossprod(Cmat) + wts[2] * crossprod(Dmat)
            Amat <- Qmat %*% smpower(Tmat, power = -1)
          }
        } else {
          Amat <- t(cmls(X = CkrB, Y = t(Xa), const = const[1],  struc = Astruc,
                         df = control$df[1], degree = control$degree[1],
                         intercept = control$intercept[1],
                         XtX = crossprod(Cmat) * crossprod(Bmat),
                         backfit = backfit, mode.range = Amodes))
          if(any(colSums(abs(Amat)) <= xtol[1])) cflag <- 2
        }
        
      } # end if(is.na(cflag))
      
      # Step 2: update mode B weights
      if(is.null(Bfixed) && is.na(cflag)){
        
        for(u in 1:nfac) CkrA[,u] <- kronecker(Cmat[,u], Amat[,u])
        if(mode == 2L){
          Qmat <- wts[1] * Xb %*% CkrA + wts[2] * projY %*% Dmat
          if(const[2] == "orthog"){
            Qsvd <- svd(Qmat)
            Bmat <- tcrossprod(Qsvd$u, Qsvd$v)
          } else {
            Tmat <- wts[1] * crossprod(Amat) * crossprod(Cmat) + wts[2] * crossprod(Dmat)
            Bmat <- Qmat %*% smpower(Tmat, power = -1)
          }
        } else {
          Bmat <- t(cmls(X = CkrA, Y = t(Xb), const = const[2], struc = Bstruc,
                         df = control$df[2], degree = control$degree[2],
                         intercept = control$intercept[2], 
                         XtX = crossprod(Cmat) * crossprod(Amat),
                         backfit = backfit, mode.range = Bmodes))
          if(any(colSums(abs(Bmat)) <= xtol[2])) cflag <- 2
        }
        
      } # end if(is.na(cflag))
      
      # Step 3: update mode C weights
      if(is.null(Cfixed) && is.na(cflag)){
        
        for(u in 1:nfac) BkrA[,u] <- kronecker(Bmat[,u], Amat[,u])
        if(mode == 3L){
          Qmat <- wts[1] * Xc %*% BkrA + wts[2] * projY %*% Dmat
          if(const[3] == "orthog"){
            Qsvd <- svd(Qmat)
            Cmat <- tcrossprod(Qsvd$u, Qsvd$v)
          } else {
            Tmat <- wts[1] * crossprod(Amat) * crossprod(Bmat) + wts[2] * crossprod(Dmat)
            Cmat <- Qmat %*% smpower(Tmat, power = -1)
          }
        } else {
          Cmat <- t(cmls(X = BkrA, Y = t(Xc), const = const[3], struc = Cstruc,
                         df = control$df[3], degree = control$degree[3],
                         intercept = control$intercept[3], 
                         XtX = crossprod(Bmat) * crossprod(Amat),
                         backfit = backfit, mode.range = Cmodes))
          if(any(colSums(abs(Cmat)) <= xtol[3])) cflag <- 2
        }
        
      } # end if(is.na(cflag))
      
      # Step 4: update regression coefficients
      if(is.null(Dfixed) && is.na(cflag)){
        if(mode == 1L){
          Dmat <- t(cmls(X = Amat, Y = dataY, const = const[4], struc = Dstruc,
                         df = control$df[4], degree = control$degree[4],
                         intercept = control$intercept[4],
                         backfit = backfit, mode.range = Dmodes))
        } else if(mode == 2L) {
          Dmat <- t(cmls(X = Bmat, Y = dataY, const = const[4], struc = Dstruc,
                         df = control$df[4], degree = control$degree[4],
                         intercept = control$intercept[4],
                         backfit = backfit, mode.range = Dmodes))
        } else if(mode == 3L) {
          Dmat <- t(cmls(X = Cmat, Y = dataY, const = const[4], struc = Dstruc,
                         df = control$df[4], degree = control$degree[4],
                         intercept = control$intercept[4],
                         backfit = backfit, mode.range = Dmodes))
          
        }
        if(any(colSums(abs(Dmat)) <= xtol[4])) cflag <- 2
      } # end if(is.na(cflag))
      
      # Step 5: check for convergence
      for(u in 1:nfac) CkrB[,u] <- kronecker(Cmat[,u], Bmat[,u])
      ssenew1 <- sum((Xa - tcrossprod(Amat, CkrB))^2)
      if(mode == 1L){
        ssenew2 <- sum((dataY - tcrossprod(Amat, Dmat))^2)
      } else if(mode == 2L){
        ssenew2 <- sum((dataY - tcrossprod(Bmat, Dmat))^2)
      } else {
        ssenew2 <- sum((dataY - tcrossprod(Cmat, Dmat))^2)
      }
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
    
    ### scale and order solution
    if(is.null(Afixed) && is.null(Bfixed) && is.null(Cfixed) && is.null(Dfixed)){
      
      # put the scale in Mode C
      adg <- colMeans(Amat^2)
      asg <- sign(colSums(Amat^3))
      Amat <- Amat %*% (diag(nfac) * (adg^-0.5) * asg)
      bdg <- colMeans(Bmat^2)
      bsg <- sign(colSums(Bmat^3))
      Bmat <- Bmat %*% (diag(nfac) * (bdg^-0.5) * bsg)
      Cmat <- Cmat %*% (diag(nfac) * ((adg*bdg)^0.5) * (asg*bsg))
      if(mode == 1L){
        Dmat <- Dmat %*% (diag(nfac) * (adg^0.5) * asg)
        Wmat <- Wmat %*% (diag(nfac) * (adg^-0.5) * asg)
      } else if(mode == 2L){
        Dmat <- Dmat %*% (diag(nfac) * (bdg^0.5) * bsg)
        Wmat <- Wmat %*% (diag(nfac) * (bdg^-0.5) * bsg)
      } else if(mode == 3L){
        Dmat <- Dmat %*% (diag(nfac) * ((adg*bdg)^-0.5) * (asg*bsg))
        Wmat <- Wmat %*% (diag(nfac) * ((adg*bdg)^0.5) * (asg*bsg))
      }
      
      # order according to sum-of-squares
      if(is.null(Astruc) && is.null(Bstruc) && is.null(Cstruc) && is.null(Dstruc)){
        fordr <- order(colSums(Cmat^2), decreasing = TRUE)
        Amat <- Amat[,fordr,drop=FALSE]
        Bmat <- Bmat[,fordr,drop=FALSE]
        Cmat <- Cmat[,fordr,drop=FALSE]
        Dmat <- Dmat[,fordr,drop=FALSE]
        Wmat <- Wmat[,fordr,drop=FALSE]
      }
      
    }
    
    ### collect results
    SSE <- c(X = ssenew1, Y = ssenew2)
    RsqX <- 1 - ssenew1 / ssx
    RsqY <- 1 - ssenew2 / ssy
    Rsq <- c(X = RsqX, Y = RsqY)
    if(is.na(cflag)) cflag <- ifelse(vtol <= ctol, 0, 1)
    fixed <- c(ifelse(is.null(Afixed), FALSE, TRUE), ifelse(is.null(Bfixed), FALSE, TRUE), 
               ifelse(is.null(Cfixed), FALSE, TRUE), ifelse(is.null(Dfixed), FALSE, TRUE))
    struc <- c(ifelse(is.null(Astruc), FALSE, TRUE), ifelse(is.null(Bstruc), FALSE, TRUE), 
               ifelse(is.null(Cstruc), FALSE, TRUE), ifelse(is.null(Dstruc), FALSE, TRUE))
    pfac <- list(A = Amat, B = Bmat, C = Cmat, D = Dmat, W = Wmat, 
                 LOSS = msenew, SSE = SSE, Rsq = Rsq, iter = iter, 
                 cflag = cflag, model = "parafac", const = const, 
                 control = control, weights = NULL, alpha = alpha, 
                 fixed = fixed, struc = struc)
    return(pfac)
    
  } # end mcr_parafac