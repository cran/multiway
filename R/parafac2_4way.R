parafac2_4way <- 
  function(data, nfac, xcx = sumsq(data), const = rep("uncons", 4),
           control = const.control(const), maxit = 500, ctol = 1e-4, 
           Gfixed = NULL, Bfixed = NULL, Cfixed = NULL, Dfixed = NULL, 
           Gstart = NULL, Bstart = NULL, Cstart = NULL, Dstart = NULL, 
           Gstruc = NULL, Bstruc = NULL, Cstruc = NULL, Dstruc = NULL,
           Gmodes = NULL, Bmodes = NULL, Cmodes = NULL, Dmodes = NULL,
           backfit = FALSE){
    # 4-way Parallel Factor Analysis 2 (Parafac2)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 25, 2018
    
    ### initialize Khatri-Rao product matrices
    nx <- sapply(data, function(x) dim(x)[1])
    xdims <- rep(NA, 4)
    xdims[2] <- dim(data[[1]])[2]
    xdims[3] <- dim(data[[1]])[3]
    xdims[4] <- length(data)
    CkrB <- matrix(0, nrow = xdims[2] * xdims[3], ncol = nfac)
    DCkrB <- matrix(0, nrow = xdims[2] * xdims[3] * xdims[4], ncol = nfac)
    if(is.null(Bfixed)) DCkrA <- matrix(0, nrow = nfac * xdims[3] * xdims[4], ncol = nfac)
    if(is.null(Cfixed)) DBkrA <- matrix(0, nrow = nfac * xdims[2] * xdims[4], ncol = nfac)
    if(is.null(Dfixed)) CBkrA <- matrix(0, nrow = nfac * xdims[2] * xdims[3], ncol = nfac)
    
    ### initialize stuff for Mode A update
    Amat <- vector("list", xdims[4])
    Xtilde <- array(0, dim = c(nfac, xdims[2], xdims[3], xdims[4]))
    
    ### reshape raw data
    for(kk in 1:xdims[4]){
      mdim <- dim(data[[kk]])
      data[[kk]] <- matrix(data[[kk]], nrow = mdim[1], ncol = xdims[2] * xdims[3])
    }
    
    ### initialize parameter matrices
    updateGmat <- FALSE
    oconst <- const
    Gmat <- Gfixed
    if(is.null(Gmat)){
      if(any(const[1] == c("orthog", "ortnon", "ortsmo", "orsmpe"))){
        Gmat <- diag(nfac)
      } else {
        updateGmat <- TRUE
        Gmat <- Gstart
        if(is.null(Gmat)) Gmat <- matrix(rnorm(nfac^2), nrow = nfac, ncol = nfac)
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
      Dmat <- initcmls(nobs = xdims[4], nfac = nfac, const = const[4],
                       struc = Dstruc, df = control$df[4], 
                       degree = control$degree[4], 
                       intercept = control$intercept[4],
                       mode.range = Dmodes)
    }
    
    ### transpose struc arguments (if given)
    tGstruc <- NULL
    if(!is.null(Gstruc)) tGstruc <- t(Gstruc)
    if(!is.null(Bstruc)) Bstruc <- t(Bstruc)
    if(!is.null(Cstruc)) Cstruc <- t(Cstruc)
    if(!is.null(Dstruc)) Dstruc <- t(Dstruc)
    
    ### iterative update of matrices
    xtol <- c(nfac, xdims[2:4]) * .Machine$double.eps
    vtol <- sseold <- xcx + ctol
    iter <- 0
    cflag <- NA
    while((vtol > ctol) && (iter < maxit)) {
      
      ## Step 1: update mode A weights
      if(is.na(cflag)){
        
        # 1a: update orthogonal projections
        for(u in 1:nfac) CkrB[,u] <- kronecker(Cmat[,u], Bmat[,u])
        if(any(const[1] == c("uncons", "orthog"))){
          for(kk in 1:xdims[4]){
            xsvd <- svd(data[[kk]] %*% CkrB %*% tcrossprod((diag(nfac)*Dmat[kk,]), Gmat))
            Amat[[kk]] <- tcrossprod(xsvd$u, xsvd$v)
            attr(Amat[[kk]], "df") <- rep(nx[kk] - (nfac + 1) / 2, nfac)
            Xtilde[,,,kk] <- array(crossprod(Amat[[kk]], data[[kk]]), dim = c(nfac, xdims[2], xdims[3]))
          }
        } else {
          for(kk in 1:xdims[4]){
            Bk <- CkrB %*% tcrossprod((diag(nfac)*Dmat[kk,]), Gmat)
            Amat[[kk]] <- t(cmls(X = Bk, Y = t(data[[kk]]), const = const[1], 
                                 df = control$df[1], degree = control$degree[1],
                                 intercept = control$intercept[1]))
            if(const[1] == "ortnon"){
              ssRkk <- colSums(Amat[[kk]]^2)
              newdf <- attr(Amat[[kk]], "df")
              for(ll in 1:nfac){
                if(ssRkk[ll] > 0){
                  Amat[[kk]][,ll] <- Amat[[kk]][,ll] / sqrt(ssRkk[ll])
                  Dmat[kk,ll] <- Dmat[kk,ll] * sqrt(ssRkk[ll])
                  newdf[ll] <- newdf[ll] - 1L
                }
              }
              attr(Amat[[kk]], "df") <- newdf
            } # end if(const[1] == "ortnon")
            Xtilde[,,,kk] <- array(crossprod(Amat[[kk]], data[[kk]]), dim = c(nfac, xdims[2], xdims[3]))
          } # end for(kk in 1:xdims[4])
        } # end if(any(const[1] == c("uncons", "orthog")))
        
        # 1b: update correlation matrix
        if(updateGmat){
          Xa <- matrix(Xtilde, nrow = nfac, ncol = xdims[2] * xdims[3] * xdims[4])
          for(u in 1:nfac) DCkrB[,u] <- kronecker(Dmat[,u], kronecker(Cmat[,u], Bmat[,u]))
          cpmat <- crossprod(Dmat) * crossprod(Cmat) * crossprod(Bmat)
          Gmat <- t(cmls(X = DCkrB, Y = t(Xa), const = "uncons", struc = tGstruc, XtX = cpmat))
          if(any(colSums(abs(Gmat)) <= xtol[1])) cflag <- 2
        }
        
      } # end if(is.na(cflag))
      
      ## Step 2: update mode B weights
      if(is.null(Bfixed) && is.na(cflag)){
        Xb <- matrix(aperm(Xtilde, perm = c(2,1,3,4)), nrow = xdims[2], ncol = nfac * xdims[3] * xdims[4])
        for(u in 1:nfac) DCkrA[,u] <- kronecker(Dmat[,u], kronecker(Cmat[,u], Gmat[,u]))
        cpmat <- crossprod(Dmat) * crossprod(Cmat) * crossprod(Gmat)
        Bmat <- t(cmls(X = DCkrA, Y = t(Xb), const = const[2], struc = Bstruc,
                       df = control$df[2], degree = control$degree[2],
                       intercept = control$intercept[2], XtX = cpmat,
                       backfit = backfit, mode.range = Bmodes))
        if(any(colSums(abs(Bmat)) <= xtol[2])) cflag <- 2
      }
      
      ## Step 3: update mode C weights
      if(is.null(Cfixed) && is.na(cflag)){
        Xc <- matrix(aperm(Xtilde, perm = c(3,1,2,4)), nrow = xdims[3], ncol = nfac * xdims[2] * xdims[4])
        for(u in 1:nfac) DBkrA[,u] <- kronecker(Dmat[,u], kronecker(Bmat[,u], Gmat[,u]))
        cpmat <- crossprod(Dmat) * crossprod(Bmat) * crossprod(Gmat)
        Cmat <- t(cmls(X = DBkrA, Y = t(Xc), const = const[3], struc = Cstruc,
                       df = control$df[3], degree = control$degree[3],
                       intercept = control$intercept[3], XtX = cpmat,
                       backfit = backfit, mode.range = Cmodes))
        if(any(colSums(abs(Cmat)) <= xtol[3])) cflag <- 2
      }
      
      ## Step 4: update mode D weights
      if(is.null(Dfixed) && is.na(cflag)){
        Xd <- matrix(aperm(Xtilde, perm = c(4,1,2,3)), nrow = xdims[4], ncol = nfac * xdims[2] * xdims[3])
        for(u in 1:nfac) CBkrA[,u] <- kronecker(Cmat[,u], kronecker(Bmat[,u], Gmat[,u]))
        cpmat <- crossprod(Cmat) * crossprod(Bmat) * crossprod(Gmat)
        Dmat <- t(cmls(X = CBkrA, Y = t(Xd), const = const[4], struc = Dstruc,
                       df = control$df[4], degree = control$degree[4],
                       intercept = control$intercept[4], XtX = cpmat,
                       backfit = backfit, mode.range = Dmodes))
        if(any(colSums(abs(Dmat)) <= xtol[4])) cflag <- 2
      }
      
      ## Step 5: check for convergence
      for(u in 1:nfac) CkrB[,u] <- kronecker(Cmat[,u], Bmat[,u])
      ssenew <- 0
      for(kk in 1:xdims[4]){
        ssenew <- ssenew + sum((data[[kk]] - tcrossprod(Amat[[kk]] %*% Gmat %*% (diag(nfac)*Dmat[kk,]), CkrB))^2)
      }
      vtol <- (sseold - ssenew) / xcx
      sseold <- ssenew
      iter <- iter + 1
      
    }  # end while(vtol>ctol && iter<maxit)
    
    ### effective degrees of freedom (uncorrected)
    Adf <- sum(sapply(Amat, function(x) sum(attr(x, "df"))))
    Gdf <- ifelse(is.null(Gfixed), ifelse(updateGmat, nfac * (nfac + 1) / 2, nfac), 0)
    Bdf <- ifelse(is.null(Bfixed), sum(attr(Bmat, "df")), 0)
    Cdf <- ifelse(is.null(Cfixed), sum(attr(Cmat, "df")), 0)
    Ddf <- ifelse(is.null(Dfixed), sum(attr(Dmat, "df")), 0)
    
    ### correct Gmat (if structured)
    if(updateGmat && !is.null(Gstruc)){
      uptri <- upper.tri(Gstruc, diag = TRUE)
      Phistruc <- crossprod(Gstruc)
      Gdf <- Gdf - sum(Phistruc[uptri] == 0)
    }
    
    ### add back scale for orthogonal constraints
    const.orthog <- c("orthog", "ortsmo", "orsmpe")
    if(is.null(Bfixed) && any(const[2] == const.orthog)) Bdf <- Bdf + nfac
    if(is.null(Cfixed) && any(const[3] == const.orthog)) Cdf <- Cdf + nfac
    if(is.null(Dfixed) && any(const[4] == const.orthog)) Ddf <- Ddf + nfac
    
    ### correct effective degrees of freedom
    fixedID <- c(is.null(Gfixed), is.null(Bfixed), is.null(Cfixed), is.null(Dfixed))
    nfixed <- 4L - sum(fixedID)
    if(nfixed == 0L){
      Gdf <- Gdf - nfac
      Bdf <- Bdf - nfac
      Cdf <- Cdf - nfac
    } else if(nfixed == 1L){
      if(!is.null(Gfixed)){
        Bdf <- Bdf - nfac
        Cdf <- Cdf - nfac
      } else if(!is.null(Bfixed)){
        Gdf <- Gdf - nfac
        Cdf <- Cdf - nfac
      } else {
        Gdf <- Gdf - nfac
        Bdf <- Bdf - nfac
      }
    } else if(nfixed == 2L){
      if(is.null(Gfixed)){
        Gdf <- Gdf - nfac
      } else if(is.null(Bfixed)){
        Bdf <- Bdf - nfac
      } else {
        Cdf <- Cdf - nfac
      }
    }
    # no change needed if nfixed == 3L (or nfixed == 4L)
    
    ### scale and order solution
    fixedNULL <- is.null(Gfixed) && is.null(Bfixed) && is.null(Cfixed) && is.null(Dfixed)
    if(fixedNULL){
      
      # put the scale in Mode D
      adg <- colSums(Gmat^2)
      if(any(adg == 0)) adg[adg == 0] <- 1
      Gmat <- Gmat %*% (diag(nfac)*(adg^-0.5))
      bdg <- colMeans(Bmat^2)
      if(any(bdg == 0)) bdg[bdg == 0] <- 1
      Bmat <- Bmat %*% (diag(nfac)*(bdg^-0.5))
      cdg <- colMeans(Cmat^2)
      if(any(cdg == 0)) cdg[cdg == 0] <- 1
      Cmat <- Cmat %*% (diag(nfac)*(cdg^-0.5))
      Dmat <- Dmat %*% (diag(nfac)*((adg*bdg*cdg)^0.5))
      
      # order according to sum-of-squares
      strucNULL <- is.null(Gstruc) && is.null(Bstruc) && is.null(Cstruc) && is.null(Dstruc)
      modesNULL <- is.null(Gmodes) && is.null(Bmodes) && is.null(Cmodes) && is.null(Dmodes)
      if(strucNULL && modesNULL){
        fordr <- order(colSums(Dmat^2), decreasing=TRUE)
        Gmat <- Gmat[,fordr,drop=FALSE]
        Bmat <- Bmat[,fordr,drop=FALSE]
        Cmat <- Cmat[,fordr,drop=FALSE]
        Dmat <- Dmat[,fordr,drop=FALSE]
      }
      
    }
    
    ### GCV criterion
    edf <- c(Adf + Gdf, Bdf, Cdf, Ddf)
    names(edf) <- LETTERS[1:4]
    pxdim <- sum(nx) * prod(xdims[2:3])
    GCV <- (ssenew / pxdim) / (1 - sum(edf) / pxdim)^2
    Rsq <- 1 - ssenew / xcx
    
    ### collect results
    if(is.na(cflag)) cflag <- ifelse(vtol <= ctol, 0, 1)
    fixed <- c(ifelse(is.null(Gfixed), FALSE, TRUE), ifelse(is.null(Bfixed), FALSE, TRUE), 
               ifelse(is.null(Cfixed), FALSE, TRUE), ifelse(is.null(Dfixed), FALSE, TRUE))
    struc <- c(ifelse(is.null(Gstruc), FALSE, TRUE), ifelse(is.null(Bstruc), FALSE, TRUE), 
               ifelse(is.null(Cstruc), FALSE, TRUE), ifelse(is.null(Dstruc), FALSE, TRUE))
    for(k in 1:xdims[4]) Amat[[k]] <- Amat[[k]] %*% Gmat
    pfac <- list(A = Amat, B = Bmat, C = Cmat, D = Dmat, Phi = crossprod(Gmat),
                 SSE = ssenew, Rsq = Rsq, GCV = GCV, edf = edf, iter = iter,
                 cflag = cflag, const = oconst, control = control, 
                 fixed = fixed, struc = struc)
    return(pfac)
    
  }