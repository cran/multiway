parafac2_3wayna <- 
  function(data, nfac, naid = NULL, const = rep("uncons", 3),
           control = const.control(const), maxit = 500, ctol = 1e-4, 
           Gfixed = NULL, Bfixed = NULL, Cfixed = NULL, 
           Gstart = NULL, Bstart = NULL, Cstart = NULL,
           Gstruc = NULL, Bstruc = NULL, Cstruc = NULL,
           Gmodes = NULL, Bmodes = NULL, Cmodes = NULL, 
           backfit = FALSE){
    # 3-way Parallel Factor Analysis 2 (Parafac2)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 25, 2018
    
    ### initialize Khatri-Rao product matrices
    nx <- sapply(data, nrow)
    xdims <- rep(NA, 3)
    xdims[2] <- ncol(data[[1]])
    xdims[3] <- length(data)
    CkrB <- matrix(0, nrow = xdims[2] * xdims[3], ncol = nfac)
    if(is.null(Bfixed)) CkrA <- matrix(0, nrow = nfac * xdims[3], ncol = nfac)
    if(is.null(Cfixed)) BkrA <- matrix(0, nrow = nfac * xdims[2], ncol = nfac)
    
    ### initialize missing data
    xcx.nona <- sumsq(data, na.rm = TRUE)
    if(is.null(naid)) naid <- lapply(data, function(x) which(is.na(x)))
    nmiss <- sapply(naid, length)
    for(k in 1:xdims[3]){
      if(nmiss[k] > 0) data[[k]][naid[[k]]] <- rnorm(nmiss[k])
    }
    xcx <- sumsq(data)
    Xhat <- vector("list", xdims[3])
    
    ### initialize stuff for Mode A update
    Amat <- vector("list", xdims[3])
    Xtilde <- array(0, dim = c(nfac, xdims[2], xdims[3]))
    
    ### initialize parameter matrices
    updateGmat <- FALSE
    oconst <- const
    Gmat <- Gfixed
    if(is.null(Gmat)){
      if(any(const[1] == c("orthog", "ortnon", "ortsmo",  "orsmpe"))){
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
    
    ### transpose struc arguments (if given)
    tGstruc <- NULL
    if(!is.null(Gstruc)) tGstruc <- t(Gstruc)
    if(!is.null(Bstruc)) Bstruc <- t(Bstruc)
    if(!is.null(Cstruc)) Cstruc <- t(Cstruc)
    
    ### iterative update of matrices
    xtol <- c(nfac, xdims[2:3]) * .Machine$double.eps
    vtol <- sseold <- xcx + ctol
    iter <- 0
    cflag <- NA
    while((vtol > ctol) && (iter < maxit)) {
      
      ## Step 1: update mode A weights
      if(is.na(cflag)){
        
        # 1a: update orthogonal projections
        if(any(const[1] == c("uncons", "orthog"))){
          for(kk in 1:xdims[3]){
            xsvd <- svd(data[[kk]] %*% Bmat %*% tcrossprod((diag(nfac)*Cmat[kk,]), Gmat))
            Amat[[kk]] <- tcrossprod(xsvd$u, xsvd$v)
            attr(Amat[[kk]], "df") <- rep(nx[kk] - (nfac + 1) / 2, nfac)
            Xtilde[,,kk] <- crossprod(Amat[[kk]], data[[kk]])
          }
        } else {
          for(kk in 1:xdims[3]){
            Bk <- Bmat %*% tcrossprod((diag(nfac)*Cmat[kk,]), Gmat)
            Amat[[kk]] <- t(cmls(X = Bk, Y = t(data[[kk]]), const = const[1], 
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
            Xtilde[,,kk] <- crossprod(Amat[[kk]], data[[kk]])
          } # end for(kk in 1:xdims[3])
        } # end if(any(const[1] == c("uncons", "orthog")))
        
        # 1b: update correlation matrix
        if(updateGmat){
          Xa <- matrix(Xtilde, nrow = nfac, ncol = xdims[2] * xdims[3])
          for(u in 1:nfac) CkrB[,u] <- kronecker(Cmat[,u], Bmat[,u])
          cpmat <- crossprod(Cmat) * crossprod(Bmat)
          Gmat <- t(cmls(X = CkrB, Y = t(Xa), const = "uncons", struc = tGstruc, XtX = cpmat))
          if(any(colSums(abs(Gmat)) <= xtol[1])) cflag <- 2
        }
        
      } # end if(is.na(cflag))
      
      ## Step 2: update mode B weights
      if(is.null(Bfixed) && is.na(cflag)){
        Xb <- matrix(aperm(Xtilde, perm = c(2,1,3)), nrow = xdims[2], ncol = nfac * xdims[3])
        for(u in 1:nfac) CkrA[,u] <- kronecker(Cmat[,u], Gmat[,u])
        cpmat <- crossprod(Cmat) * crossprod(Gmat)
        Bmat <- t(cmls(X = CkrA, Y = t(Xb), const = const[2], struc = Bstruc,
                       df = control$df[2], degree = control$degree[2],
                       intercept = control$intercept[2], XtX = cpmat,
                       backfit = backfit, mode.range = Bmodes))
        if(any(colSums(abs(Bmat)) <= xtol[2])) cflag <- 2
      }
      
      ## Step 3: update mode C weights
      if(is.null(Cfixed) && is.na(cflag)){
        Xc <- matrix(aperm(Xtilde, perm = c(3,1,2)), nrow = xdims[3], ncol = nfac * xdims[2])
        for(u in 1:nfac) BkrA[,u] <- kronecker(Bmat[,u], Gmat[,u])
        cpmat <- crossprod(Bmat) * crossprod(Gmat)
        Cmat <- t(cmls(X = BkrA, Y = t(Xc), const = const[3], struc = Cstruc,
                       df = control$df[3], degree = control$degree[3],
                       intercept = control$intercept[3], XtX = cpmat,
                       backfit = backfit, mode.range = Cmodes))
        if(any(colSums(abs(Cmat)) <= xtol[3])) cflag <- 2
      }
      
      ## Step 4: check for convergence
      ssenew <- 0
      for(kk in 1:xdims[3]){
        Xhat[[kk]] <- tcrossprod(Amat[[kk]] %*% Gmat %*% (diag(nfac)*Cmat[kk,]), Bmat)
        ssenew <- ssenew + sum((data[[kk]] - Xhat[[kk]])^2)
      }
      vtol <- (sseold - ssenew) / xcx
      sseold <- ssenew
      iter <- iter + 1
      
      # impute missing data
      for(k in 1:xdims[3]){
        if(nmiss[k] > 0) data[[k]][naid[[k]]] <- Xhat[[k]][naid[[k]]]
      }
      xcx <- sumsq(data)
      
    }  # end while(vtol>ctol && iter<maxit)
    
    ### update SSE
    ssenew <- 0
    for(kk in 1:xdims[3]) ssenew <- ssenew + sum((data[[kk]] - Xhat[[kk]])^2)
    
    ### effective degrees of freedom (uncorrected)
    Adf <- sum(sapply(Amat, function(x) sum(attr(x, "df"))))
    Gdf <- ifelse(is.null(Gfixed), ifelse(updateGmat, nfac * (nfac + 1) / 2, nfac), 0)
    Bdf <- ifelse(is.null(Bfixed), sum(attr(Bmat, "df")), 0)
    Cdf <- ifelse(is.null(Cfixed), sum(attr(Cmat, "df")), 0)
    
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
    
    ### correct effective degrees of freedom
    fixedID <- c(is.null(Gfixed), is.null(Bfixed), is.null(Cfixed))
    nfixed <- 3L - sum(fixedID)
    if(nfixed == 0L){
      Gdf <- Gdf - nfac
      Bdf <- Bdf - nfac
    } else if(nfixed == 1L){
      if(!is.null(Gfixed)){
        Bdf <- Bdf - nfac
      } else {
        Gdf <- Gdf - nfac
      }
    }
    # no change needed if nfixed == 2L (or nfixed == 3L)
    
    ### scale and order solution
    fixedNULL <- is.null(Gfixed) & is.null(Bfixed) & is.null(Cfixed)
    if(fixedNULL){
      
      # put the scale in Mode C
      adg <- colSums(Gmat^2)
      if(any(adg == 0)) adg[adg == 0] <- 1
      Gmat <- Gmat %*% (diag(nfac)*(adg^-0.5))
      bdg <- colMeans(Bmat^2)
      if(any(bdg == 0)) bdg[bdg == 0] <- 1
      Bmat <- Bmat %*% (diag(nfac)*(bdg^-0.5))
      Cmat <- Cmat %*% (diag(nfac)*((adg*bdg)^0.5))
      
      # order according to sum-of-squares
      strucNULL <- is.null(Gstruc) && is.null(Bstruc) && is.null(Cstruc)
      modesNULL <- is.null(Gmodes) && is.null(Bmodes) && is.null(Cmodes)
      if(strucNULL && modesNULL){
        fordr <- order(colSums(Cmat^2), decreasing=TRUE)
        Gmat <- Gmat[,fordr,drop=FALSE]
        Bmat <- Bmat[,fordr,drop=FALSE]
        Cmat <- Cmat[,fordr,drop=FALSE]
      }
      
    }
    
    ### GCV criterion
    edf <- c(Adf + Gdf, Bdf, Cdf)
    names(edf) <- LETTERS[1:3]
    pxdim <- sum(nx) * xdims[2]
    GCV <- (ssenew / pxdim) / (1 - sum(edf) / pxdim)^2
    Rsq <- 1 - ssenew / xcx.nona
    
    ### collect results
    if(is.na(cflag)) cflag <- ifelse(vtol <= ctol, 0, 1)
    fixed <- c(ifelse(is.null(Gfixed), FALSE, TRUE), ifelse(is.null(Bfixed), FALSE, TRUE), ifelse(is.null(Cfixed), FALSE, TRUE))
    struc <- c(ifelse(is.null(Gstruc), FALSE, TRUE), ifelse(is.null(Bstruc), FALSE, TRUE), ifelse(is.null(Cstruc), FALSE, TRUE))
    for(k in 1:xdims[3]) Amat[[k]] <- Amat[[k]] %*% Gmat
    pfac <- list(A = Amat, B = Bmat, C = Cmat, Phi = crossprod(Gmat), 
                 SSE = ssenew, Rsq = Rsq, GCV = GCV, edf = edf,
                 iter = iter, cflag = cflag, const = oconst,
                 control = control, fixed = fixed, struc = struc)
    return(pfac)
    
  }