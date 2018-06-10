parafac_3way <-
  function(data, nfac, xcx = sumsq(data), const = rep("uncons", 3),
           maxit = 500, ctol = 1e-4, control = const.control(const),
           Afixed = NULL, Bfixed = NULL, Cfixed = NULL,
           Astart = NULL, Bstart = NULL, Cstart = NULL, 
           Astruc = NULL, Bstruc = NULL, Cstruc = NULL,
           Amodes = NULL, Bmodes = NULL, Cmodes = NULL,
           backfit = FALSE){
    # 3-way Parallel Factor Analysis (Parafac)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 25, 2018
    
    ### initialize Khatri-Rao product matrices
    xdims <- dim(data)
    CkrB <- matrix(0, nrow = xdims[2] * xdims[3], ncol = nfac)
    if(is.null(Bfixed)) CkrA <- matrix(0, nrow = xdims[1] * xdims[3], ncol = nfac)
    if(is.null(Cfixed)) BkrA <- matrix(0, nrow = xdims[1] * xdims[2], ncol = nfac)
    
    ### initialize reshaped data matrices
    Xa <- t(matrix(data, nrow = xdims[1], ncol = xdims[2] * xdims[3]))
    if(is.null(Bfixed)) Xb <- t(matrix(aperm(data, perm = c(2,1,3)), nrow = xdims[2], ncol = xdims[1] * xdims[3]))
    if(is.null(Cfixed)) Xc <- t(matrix(aperm(data, perm = c(3,1,2)), nrow = xdims[3], ncol = xdims[1] * xdims[2]))
    rm(data)
    
    ### initialize parameter matrices
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
    
    ### transpose struc arguments (if given)
    if(!is.null(Astruc)) Astruc <- t(Astruc)
    if(!is.null(Bstruc)) Bstruc <- t(Bstruc)
    if(!is.null(Cstruc)) Cstruc <- t(Cstruc)
    
    ### iterative update of matrices
    xtol <- xdims * .Machine$double.eps
    vtol <- sseold <- xcx + ctol
    iter <- 0
    cflag <- NA
    while((vtol > ctol) && (iter < maxit)) {
      
      # Step 1: update mode A weights
      if(is.null(Afixed) && is.na(cflag)){
        for(u in 1:nfac) CkrB[,u] <- kronecker(Cmat[,u], Bmat[,u])
        cpmat <- crossprod(Cmat) * crossprod(Bmat)
        Amat <- t(cmls(X = CkrB, Y = Xa, const = const[1], struc = Astruc,
                       df = control$df[1], degree = control$degree[1],
                       intercept = control$intercept[1], XtX = cpmat, 
                       backfit = backfit, mode.range = Amodes))
        if(any(colSums(abs(Amat)) <= xtol[1])) cflag <- 2
      }
      
      # Step 2: update mode B weights
      if(is.null(Bfixed) && is.na(cflag)){
        for(u in 1:nfac) CkrA[,u] <- kronecker(Cmat[,u], Amat[,u])
        cpmat <- crossprod(Cmat) * crossprod(Amat)
        Bmat <- t(cmls(X = CkrA, Y = Xb, const = const[2], struc = Bstruc,
                       df = control$df[2], degree = control$degree[2],
                       intercept = control$intercept[2], XtX = cpmat, 
                       backfit = backfit, mode.range = Bmodes))
        if(any(colSums(abs(Bmat)) <= xtol[2])) cflag <- 2
      }
      
      # Step 3: update mode C weights
      if(is.null(Cfixed) && is.na(cflag)){
        for(u in 1:nfac) BkrA[,u] <- kronecker(Bmat[,u], Amat[,u])
        cpmat <- crossprod(Bmat) * crossprod(Amat)
        Cmat <- t(cmls(X = BkrA, Y = Xc, const = const[3], struc = Cstruc,
                       df = control$df[3], degree = control$degree[3],
                       intercept = control$intercept[3], XtX = cpmat, 
                       backfit = backfit, mode.range = Cmodes))
        if(any(colSums(abs(Cmat)) <= xtol[3])) cflag <- 2
      }
      
      # Step 4: check for convergence
      for(u in 1:nfac) CkrB[,u] <- kronecker(Cmat[,u], Bmat[,u])
      #ssenew <- sum((Xa - tcrossprod(Amat, CkrB))^2)
      ssenew <- sum((Xa - tcrossprod(CkrB, Amat))^2)
      vtol <- (sseold - ssenew) / xcx
      sseold <- ssenew
      iter <- iter + 1
      
    } # end while(vtol>ctol && iter<maxit)
    
    ### effective degrees of freedom (uncorrected)
    Adf <- ifelse(is.null(Afixed), sum(attr(Amat, "df")), 0)
    Bdf <- ifelse(is.null(Bfixed), sum(attr(Bmat, "df")), 0)
    Cdf <- ifelse(is.null(Cfixed), sum(attr(Cmat, "df")), 0)
    
    ### add back scale for orthogonal constraints
    const.orthog <- c("orthog", "ortsmo", "orsmpe")
    if(is.null(Afixed) && any(const[1] == const.orthog)) Adf <- Adf + nfac
    if(is.null(Bfixed) && any(const[2] == const.orthog)) Bdf <- Bdf + nfac
    if(is.null(Cfixed) && any(const[3] == const.orthog)) Cdf <- Cdf + nfac
    
    ### correct effective degrees of freedom
    fixedID <- c(is.null(Afixed), is.null(Bfixed), is.null(Cfixed))
    nfixed <- 3L - sum(fixedID)
    if(nfixed == 0L){
      Adf <- Adf - nfac
      Bdf <- Bdf - nfac
    } else if(nfixed == 1L){
      if(!is.null(Afixed)){
        Bdf <- Bdf - nfac
      } else {
        Adf <- Adf - nfac
      }
    }
    # no change needed if nfixed == 2L (and assuming nfixed != 3L)
    
    ### scale and order solution
    fixedNULL <- is.null(Afixed) && is.null(Bfixed) && is.null(Cfixed)
    if(fixedNULL){
      
      # put the scale in Mode C
      adg <- colMeans(Amat^2)
      if(any(adg == 0)) adg[adg == 0] <- 1
      Amat <- Amat %*% (diag(nfac) * (adg^-0.5))
      bdg <- colMeans(Bmat^2)
      if(any(bdg == 0)) bdg[bdg == 0] <- 1
      Bmat <- Bmat %*% (diag(nfac) * (bdg^-0.5) )
      Cmat <- Cmat %*% (diag(nfac) * ((adg*bdg)^0.5))
      
      # order according to sum-of-squares
      strucNULL <- is.null(Astruc) && is.null(Bstruc) && is.null(Cstruc)
      modesNULL <- is.null(Amodes) && is.null(Bmodes) && is.null(Cmodes)
      if(strucNULL && modesNULL){
        fordr <- order(colSums(Cmat^2), decreasing = TRUE)
        Amat <- Amat[,fordr,drop=FALSE]
        Bmat <- Bmat[,fordr,drop=FALSE]
        Cmat <- Cmat[,fordr,drop=FALSE]
      }
      
    }
    
    ### GCV and R-squared
    edf <- c(Adf, Bdf, Cdf)
    names(edf) <- LETTERS[1:3]
    pxdim <- prod(xdims)
    GCV <- (ssenew / pxdim) / (1 - sum(edf) / pxdim)^2
    Rsq <- 1 - ssenew / xcx
    
    ### collect results
    if(is.na(cflag)) cflag <- ifelse(vtol <= ctol, 0, 1)
    fixed <- c(ifelse(is.null(Afixed), FALSE, TRUE), ifelse(is.null(Bfixed), FALSE, TRUE), ifelse(is.null(Cfixed), FALSE, TRUE))
    struc <- c(ifelse(is.null(Astruc), FALSE, TRUE), ifelse(is.null(Bstruc), FALSE, TRUE), ifelse(is.null(Cstruc), FALSE, TRUE))
    pfac <- list(A = Amat, B = Bmat, C = Cmat, SSE = ssenew, Rsq = Rsq, 
                 GCV = GCV, edf = edf, iter = iter, cflag = cflag, 
                 const = const, control = control, fixed = fixed, struc = struc)
    return(pfac)
    
  }