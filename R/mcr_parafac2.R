mcr_parafac2 <- 
  function(dataX, dataY, nfac, alpha = 0.5, ssx = NULL, ssy = NULL, 
           const = rep("uncons", 4L), control = const.control(const), 
           Gfixed = NULL, Bfixed = NULL, Cfixed = NULL, Dfixed = NULL,
           Gstart = NULL, Bstart = NULL, Cstart = NULL, Dstart = NULL,
           Gstruc = NULL, Bstruc = NULL, Cstruc = NULL, Dstruc = NULL,
           Amodes = NULL, Bmodes = NULL, Cmodes = NULL, Dmodes = NULL,
           maxit = 500, ctol = 1e-4, projY = NULL, xsvd = NULL, backfit = FALSE){
    # Multi-way Covariates Regression (Parafac2 3-way balanced)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: April 19, 2018
    
    ### get dimensions
    xdims <- dim(dataX)
    ydims <- dim(dataY)
    
    ### check ssx and ssy
    if(is.null(ssx)) ssx <- sumsq(dataX)
    if(is.null(ssy)) ssy <- sumsq(dataY)
    
    ### initialize Khatri-Rao product matrices
    BkrA <- matrix(0, nrow = nfac * xdims[2], ncol = nfac)
    CkrA <- matrix(0, nrow = nfac * xdims[3], ncol = nfac)
    CkrB <- matrix(0, nrow = xdims[2] * xdims[3], ncol = nfac)
    
    ### initialize reshaped data matrices
    Xc <- matrix(aperm(dataX, perm = c(3,1,2)), nrow = xdims[3], ncol = xdims[1] * xdims[2])
    
    ### project response
    if(xdims[3] > prod(xdims[-3])){
      if(is.null(xsvd)) xsvd <- svd(Xc, nv = 0)
      if(is.null(projY)) projY <- xsvd$u %*% crossprod(xsvd$u, dataY)
    } else {
      projY <- dataY
    }
    
    ### initialize stuff for Mode A update
    Amat <- vector("list", xdims[3])
    Xtilde <- array(0, dim = c(nfac,xdims[2],xdims[3]))
    
    ### initialize matrices
    updateGmat <- FALSE
    oconst <- const
    Gmat <- Gfixed
    if(is.null(Gmat)){
      if(any(const[1] == c("orthog", "ortnon", "ortsmo",  "orsmpe"))){
        Gmat <- diag(nfac)
      } else {
        updateGmat <- TRUE
        Gmat <- Gstart
        if(is.null(Gmat)) Gmat <- matrix(rnorm(nfac^2),nfac,nfac)
        if(!is.null(Gstruc)) Gmat <- Gmat * Gstruc
        if(const[1] == "smooth") const[1] <- "ortsmo"
        if(const[1] == "smoper") const[1] <- "orsmpe"
      }
    }
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
    tGstruc <- Gstruc
    if(!is.null(Gstruc)) tGstruc <- t(Gstruc)
    if(!is.null(Bstruc)) Bstruc <- t(Bstruc)
    if(!is.null(Cstruc)) Cstruc <- t(Cstruc)
    if(!is.null(Dstruc)) Dstruc <- t(Dstruc)
    
    ### iterative update of matrices
    wts <- c(alpha / ssx, (1-alpha) / ssy)
    ssxssy <- 1 + ctol
    xtol <- c(nfac, xdims[2:3], ydims[2]) * .Machine$double.eps
    vtol <- mseold <- ssxssy
    iter <- 0
    cflag <- NA
    while((vtol > ctol) && (iter < maxit)) {
      
      # Step 1: update mode A weights
      if(is.na(cflag)){
        
        # 1a: update orthogonal projections
        for(kk in 1:xdims[3]){
          if(any(const[1] == c("ortnon", "ortsmo", "orsmpe"))){
            Bk <- Bmat %*% tcrossprod((diag(nfac)*Cmat[kk,]), Gmat)
            Amat[[kk]] <- t(cmls(X = Bk, Y = t(dataX[,,kk]), const = const[1], 
                                 df = control$df[1], degree = control$degree[1],
                                 intercept = control$intercept[1]))
            if(const[1] == "ortnon"){
              ssRkk <- colSums(Amat[[kk]]^2)
              newdf <- attr(Amat[[kk]], "df")
              for(ll in 1:nfac){
                if(ssRkk[ll] > 0){
                  Amat[[kk]][,ll] <- Amat[[kk]][,ll] / sqrt(ssRkk[ll])
                  Cmat[kk,ll] <- Cmat[kk,ll] * sqrt(ssRkk[ll])
                  newdf[ll] <- newdf[ll] - 1L
                }
              }
              attr(Amat[[kk]], "df") <- newdf
            } # end if(const[1] == "ortnon")
          } else {
            Rsvd <- svd(dataX[,,kk] %*% Bmat %*% tcrossprod((diag(nfac)*Cmat[kk,]), Gmat))
            Amat[[kk]] <- tcrossprod(Rsvd$u, Rsvd$v)
            attr(Amat[[kk]], "df") <- rep(xdims[kk] - (nfac + 1) / 2, nfac)
          } # end if(any(const[1] == c("ortnon", "ortsmo", "orsmpe")))
          Xtilde[,,kk] <- crossprod(Amat[[kk]], dataX[,,kk])
        } # end for(kk in 1:xdims[3])
        
        # 1b: update correlation matrix
        Xa <- matrix(Xtilde, nrow = nfac, ncol = xdims[2] * xdims[3])
        if(updateGmat){
          for(u in 1:nfac) CkrB[,u] <- kronecker(Cmat[,u], Bmat[,u])
          cpmat <- crossprod(Cmat) * crossprod(Bmat)
          Gmat <- t(cmls(X = CkrB, Y = t(Xa), const = "uncons", struc = tGstruc, XtX = cpmat))
          if(any(colSums(abs(Gmat)) <= xtol[1])) cflag <- 2
        }
        
      } # end if(is.na(cflag))
      
      # Step 2: update mode B weights
      if(is.null(Bfixed) && is.na(cflag)){
        Xb <- matrix(aperm(Xtilde, perm = c(2,1,3)), nrow = xdims[2], ncol = nfac * xdims[3])
        for(u in 1:nfac) CkrA[,u] <- kronecker(Cmat[,u], Gmat[,u])
        Bmat <- t(cmls(X = CkrA, Y = t(Xb), const = const[2], struc = Bstruc,
                       df = control$df[2], degree = control$degree[2],
                       intercept = control$intercept[2], 
                       XtX = crossprod(Cmat) * crossprod(Gmat),
                       backfit = backfit, mode.range = Bmodes))
        if(any(colSums(abs(Bmat)) <= xtol[2])) cflag <- 2
      } # end if(is.null(Bfixed) && is.na(cflag))
      
      # Step 3: update mode C weights
      if(is.null(Cfixed) && is.na(cflag)){
        Zc <- matrix(0,xdims[3],nfac)
        for(kk in 1:xdims[3]) Zc[kk,] <- Xc[kk,,drop=FALSE] %*% krprod(Bmat, Amat[[kk]] %*% Gmat)
        if(xdims[3] > prod(xdims[-3])) Zc <- xsvd$u %*% crossprod(xsvd$u, Zc)
        Qmat <- wts[1] * Zc + wts[2] * projY %*% Dmat
        if(const[3] == "orthog"){
          Qsvd <- svd(Qmat)
          Cmat <- tcrossprod(Qsvd$u, Qsvd$v)
        } else {
          Tmat <- wts[1] * crossprod(Gmat) * crossprod(Bmat) + wts[2] * crossprod(Dmat)
          Cmat <- Qmat %*% smpower(Tmat, power = -1)
        }
      } # end if(is.null(Cfixed) && is.na(cflag))
      
      # Step 4: update regression coefficients
      if(is.null(Dfixed) && is.na(cflag)){
        Dmat <- t(cmls(X = Cmat, Y = dataY, const = const[4], struc = Dstruc,
                       df = control$df[4], degree = control$degree[4],
                       intercept = control$intercept[4],
                       backfit = backfit, mode.range = Dmodes))
        if(any(colSums(abs(Dmat)) <= xtol[4])) cflag <- 2
      } # end if(is.null(Dfixed) && is.na(cflag))
      ssenew2 <- sum((dataY - tcrossprod(Cmat, Dmat))^2)
      
      # Step 5: check for convergence
      for(u in 1:nfac) CkrB[,u] <- kronecker(Cmat[,u], Bmat[,u])
      ssenew1 <- sum((Xa - tcrossprod(Gmat,CkrB))^2) - sum(Xtilde^2) + ssx
      msenew <- wts[1] * ssenew1 + wts[2] * ssenew2
      vtol <- (mseold - msenew) / (mseold + ctol)
      mseold <- msenew
      iter <- iter + 1
      
    } # end while(vtol>ctol && iter<maxit)
    
    ### get coefficients
    Wmat <- mpinv(Xc) %*% Cmat
    
    ### scale and order solution
    if(is.null(Gfixed) && is.null(Bfixed) && is.null(Cfixed) && is.null(Dfixed)){
      
      # put the scale in Mode C
      adg <- colSums(Gmat^2) / xdims[1]
      asg <- sign(colSums(Gmat^3))
      Gmat <- Gmat %*% (diag(nfac) * (adg^-0.5) * asg)
      bdg <- colMeans(Bmat^2)
      bsg <- sign(colSums(Bmat^3))
      Bmat <- Bmat %*% (diag(nfac) * (bdg^-0.5) * bsg)
      Cmat <- Cmat %*% (diag(nfac) * ((adg*bdg)^0.5) * (asg*bsg))
      Dmat <- Dmat %*% (diag(nfac) * ((adg*bdg)^-0.5) * (asg*bsg))
      Wmat <- Wmat %*% (diag(nfac) * ((adg*bdg)^0.5) * (asg*bsg))
      
      # order according to sum-of-squares
      if(is.null(Gstruc) && is.null(Bstruc) && is.null(Cstruc) && is.null(Dstruc)){
        fordr <- order(colSums(Cmat^2), decreasing = TRUE)
        Gmat <- Gmat[,fordr,drop=FALSE]
        Bmat <- Bmat[,fordr,drop=FALSE]
        Cmat <- Cmat[,fordr,drop=FALSE]
        Dmat <- Dmat[,fordr,drop=FALSE]
        Wmat <- Wmat[,fordr,drop=FALSE]
      }
      
    }
    
    ### make Amat list
    for(k in 1:xdims[3]) Amat[[k]] <- Amat[[k]] %*% Gmat
    Phi <- crossprod(Gmat)
    
    ### collect results
    SSE <- c(X = ssenew1, Y = ssenew2)
    RsqX <- 1 - ssenew1 / ssx
    RsqY <- 1 - ssenew2 / ssy
    Rsq <- c(X = RsqX, Y = RsqY)
    if(is.na(cflag)) cflag <- ifelse(vtol <= ctol, 0, 1)
    fixed <- c(ifelse(is.null(Gfixed), FALSE, TRUE), ifelse(is.null(Bfixed), FALSE, TRUE), 
               ifelse(is.null(Cfixed), FALSE, TRUE), ifelse(is.null(Dfixed), FALSE, TRUE))
    struc <- c(ifelse(is.null(Gstruc), FALSE, TRUE), ifelse(is.null(Bstruc), FALSE, TRUE), 
               ifelse(is.null(Cstruc), FALSE, TRUE), ifelse(is.null(Dstruc), FALSE, TRUE))
    pfac <- list(A = Amat, B = Bmat, C = Cmat, D = Dmat, W = Wmat,
                 LOSS = msenew, SSE = SSE, Rsq = Rsq, iter = iter,
                 cflag = cflag, model = "parafac2", const = oconst, 
                 control = control, weights = NULL, alpha = alpha, 
                 fixed = fixed, struc = struc, Phi = Phi)
    return(pfac)
    
  } # end mcr_parafac2