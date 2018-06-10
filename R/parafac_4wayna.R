parafac_4wayna <-
  function(data, nfac, naid = which(is.na(data)), const = rep("uncons", 4),
           maxit = 500, ctol = 1e-4, control = const.control(const),
           Afixed = NULL, Bfixed = NULL, Cfixed = NULL, Dfixed = NULL, 
           Astart = NULL, Bstart = NULL, Cstart = NULL, Dstart = NULL,
           Astruc = NULL, Bstruc = NULL, Cstruc = NULL, Dstruc = NULL,
           Amodes = NULL, Bmodes = NULL, Cmodes = NULL, Dmodes = NULL,
           backfit = FALSE){
    # 4-way Parallel Factor Analysis (Parafac)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 26, 2018
    
    ### initialize Khatri-Rao product matrices
    xdims <- dim(data)
    DCkrB <- matrix(0, nrow = xdims[2] * xdims[3] * xdims[4], ncol = nfac)
    if(is.null(Bfixed)) DCkrA <- matrix(0, nrow = xdims[1] * xdims[3] * xdims[4], ncol = nfac)
    if(is.null(Cfixed)) DBkrA <- matrix(0, nrow = xdims[1] * xdims[2] * xdims[4], ncol = nfac)
    if(is.null(Dfixed)) CBkrA <- matrix(0, nrow = xdims[1] * xdims[2] * xdims[3], ncol = nfac)
    
    ### initialize missing data
    xcx.nona <- sumsq(data, na.rm = TRUE)
    data[naid] <- rnorm(length(naid))
    xcx <- sum(data^2)
    
    ### initialize reshaped data matrices
    Xa <- matrix(data, nrow = xdims[1], ncol = xdims[2] * xdims[3] * xdims[4])
    if(is.null(Bfixed)) Xb <- matrix(aperm(data, perm = c(2,1,3,4)), nrow = xdims[2], ncol = xdims[1] * xdims[3] * xdims[4])
    if(is.null(Cfixed)) Xc <- matrix(aperm(data, perm = c(3,1,2,4)), nrow = xdims[3], ncol = xdims[1] * xdims[2] * xdims[4])
    if(is.null(Dfixed)) Xd <- matrix(aperm(data, perm = c(4,1,2,3)), nrow = xdims[4], ncol = xdims[1] * xdims[2] * xdims[3])
    
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
    if(!is.null(Astruc)) Astruc <- t(Astruc)
    if(!is.null(Bstruc)) Bstruc <- t(Bstruc)
    if(!is.null(Cstruc)) Cstruc <- t(Cstruc)
    if(!is.null(Dstruc)) Dstruc <- t(Dstruc)
    
    ### iterative update of matrices
    xtol <- xdims * .Machine$double.eps
    vtol <- sseold <- xcx + ctol
    iter <- 0
    cflag <- NA
    while((vtol > ctol) && (iter < maxit)) {
      
      # Step 1: update mode A weights
      if(is.null(Afixed) && is.na(cflag)){
        for(u in 1:nfac) DCkrB[,u] <- kronecker(Dmat[,u], kronecker(Cmat[,u], Bmat[,u]))
        cpmat <- crossprod(Dmat) * crossprod(Cmat) * crossprod(Bmat)
        Amat <- t(cmls(X = DCkrB, Y = t(Xa), const = const[1], struc = Astruc,
                       df = control$df[1], degree = control$degree[1],
                       intercept = control$intercept[1], XtX = cpmat, 
                       backfit = backfit, mode.range = Amodes))
        if(any(colSums(abs(Amat)) <= xtol[1])) cflag <- 2
      }
      
      # Step 2: update mode B weights
      if(is.null(Bfixed) && is.na(cflag)){
        for(u in 1:nfac) DCkrA[,u] <- kronecker(Dmat[,u], kronecker(Cmat[,u], Amat[,u]))
        cpmat <- crossprod(Dmat) * crossprod(Cmat) * crossprod(Amat)
        Bmat <- t(cmls(X = DCkrA, Y = t(Xb), const = const[2], struc = Bstruc,
                       df = control$df[2], degree = control$degree[2],
                       intercept = control$intercept[2], XtX = cpmat, 
                       backfit = backfit, mode.range = Bmodes))
        if(any(colSums(abs(Bmat)) <= xtol[2])) cflag <- 2
      }
      
      # Step 3: update mode C weights
      if(is.null(Cfixed) && is.na(cflag)){
        for(u in 1:nfac) DBkrA[,u] <- kronecker(Dmat[,u], kronecker(Bmat[,u], Amat[,u]))
        cpmat <- crossprod(Dmat) * crossprod(Bmat) * crossprod(Amat)
        Cmat <- t(cmls(X = DBkrA, Y = t(Xc), const = const[3], struc = Cstruc,
                       df = control$df[3], degree = control$degree[3],
                       intercept = control$intercept[3], XtX = cpmat, 
                       backfit = backfit, mode.range = Cmodes))
        if(any(colSums(abs(Cmat)) <= xtol[3])) cflag <- 2
      }
      
      # Step 4: update mode D weights
      if(is.null(Dfixed) && is.na(cflag)){
        for(u in 1:nfac) CBkrA[,u] <- kronecker(Cmat[,u], kronecker(Bmat[,u], Amat[,u]))
        cpmat <- crossprod(Cmat) * crossprod(Bmat) * crossprod(Amat)
        Dmat <- t(cmls(X = CBkrA, Y = t(Xd), const = const[4], struc = Dstruc,
                       df = control$df[4], degree = control$degree[4],
                       intercept = control$intercept[4], XtX = cpmat, 
                       backfit = backfit, mode.range = Dmodes))
        if(any(colSums(abs(Dmat)) <= xtol[4])) cflag <- 2
      }
      
      # Step 5: check for convergence
      for(u in 1:nfac) DCkrB[,u] <- kronecker(Dmat[,u], kronecker(Cmat[,u], Bmat[,u]))
      Xhat <- tcrossprod(Amat, DCkrB)
      ssenew <- sum((Xa - Xhat)^2)
      vtol <- (sseold - ssenew) / xcx
      sseold <- ssenew
      iter <- iter + 1
      
      # impute missing data
      data[naid] <- Xhat[naid]
      xcx <- sum(data^2)
      Xa <- matrix(data, nrow = xdims[1], ncol = xdims[2] * xdims[3] * xdims[4])
      if(is.null(Bfixed)) Xb <- matrix(aperm(data, perm = c(2,1,3,4)), nrow = xdims[2], ncol = xdims[1] * xdims[3] * xdims[4])
      if(is.null(Cfixed)) Xc <- matrix(aperm(data, perm = c(3,1,2,4)), nrow = xdims[3], ncol = xdims[1] * xdims[2] * xdims[4])
      if(is.null(Dfixed)) Xd <- matrix(aperm(data, perm = c(4,1,2,3)), nrow = xdims[4], ncol = xdims[1] * xdims[2] * xdims[3])
      
    }  # end while(vtol>ctol && iter<maxit)
    
    ### update SSE
    ssenew <- sum((Xa - Xhat)^2)
    
    ### effective degrees of freedom (uncorrected)
    Adf <- ifelse(is.null(Afixed), sum(attr(Amat, "df")), 0)
    Bdf <- ifelse(is.null(Bfixed), sum(attr(Bmat, "df")), 0)
    Cdf <- ifelse(is.null(Cfixed), sum(attr(Cmat, "df")), 0)
    Ddf <- ifelse(is.null(Dfixed), sum(attr(Dmat, "df")), 0)
    
    ### add back scale for orthogonal constraints
    const.orthog <- c("orthog", "ortsmo", "orsmpe")
    if(is.null(Afixed) && any(const[1] == const.orthog)) Adf <- Adf + nfac
    if(is.null(Bfixed) && any(const[2] == const.orthog)) Bdf <- Bdf + nfac
    if(is.null(Cfixed) && any(const[3] == const.orthog)) Cdf <- Cdf + nfac
    if(is.null(Dfixed) && any(const[4] == const.orthog)) Ddf <- Ddf + nfac
    
    ### correct effective degrees of freedom
    fixedID <- c(is.null(Afixed), is.null(Bfixed), is.null(Cfixed), is.null(Dfixed))
    nfixed <- 4L - sum(fixedID)
    if(nfixed == 0L){
      Adf <- Adf - nfac
      Bdf <- Bdf - nfac
      Cdf <- Cdf - nfac
    } else if(nfixed == 1L){
      if(!is.null(Afixed)){
        Bdf <- Bdf - nfac
        Cdf <- Cdf - nfac
      } else if(!is.null(Bfixed)){
        Adf <- Adf - nfac
        Cdf <- Cdf - nfac
      } else {
        Adf <- Adf - nfac
        Bdf <- Bdf - nfac
      }
    } else if(nfixed == 2L){
      if(is.null(Afixed)){
        Adf <- Adf - nfac
      } else if(is.null(Bfixed)){
        Bdf <- Bdf - nfac
      } else {
        Cdf <- Cdf - nfac
      }
    }
    # no change needed if nfixed == 3L (and assuming nfixed != 4L)
    
    ### scale and order solution
    fixedNULL <- is.null(Afixed) && is.null(Bfixed) && is.null(Cfixed) && is.null(Dfixed)
    if(fixedNULL){
      
      # put the scale in Mode D
      adg <- colMeans(Amat^2)
      if(any(adg == 0)) adg[adg == 0] <- 1
      Amat <- Amat %*% (diag(nfac) * (adg^-0.5))
      bdg <- colMeans(Bmat^2)
      if(any(bdg == 0)) bdg[bdg == 0] <- 1
      Bmat <- Bmat %*% (diag(nfac) * (bdg^-0.5))
      cdg <- colMeans(Cmat^2)
      if(any(cdg == 0)) cdg[cdg == 0] <- 1
      Cmat <- Cmat %*% (diag(nfac) * (cdg^-0.5))
      Dmat <- Dmat %*% (diag(nfac) * ((adg*bdg*cdg)^0.5))
      
      # order according to sum-of-squares
      strucNULL <- is.null(Astruc) && is.null(Bstruc) && is.null(Cstruc) && is.null(Dstruc)
      modesNULL <- is.null(Amodes) && is.null(Bmodes) && is.null(Cmodes) && is.null(Dmodes)
      if(strucNULL && modesNULL){
        fordr <- order(colSums(Dmat^2),decreasing=TRUE)
        Amat <- Amat[,fordr,drop=FALSE]
        Bmat <- Bmat[,fordr,drop=FALSE]
        Cmat <- Cmat[,fordr,drop=FALSE]
        Dmat <- Dmat[,fordr,drop=FALSE]
      }
      
    }
    
    ### GCV and R-squared
    edf <- c(Adf, Bdf, Cdf, Ddf)
    names(edf) <- LETTERS[1:4]
    pxdim <- prod(xdims)
    GCV <- (ssenew / pxdim) / (1 - sum(edf) / pxdim)^2
    Rsq <- 1 - ssenew / xcx.nona
    
    ### collect results
    if(is.na(cflag)) cflag <- ifelse(vtol <= ctol, 0, 1)
    fixed <- c(ifelse(is.null(Afixed), FALSE, TRUE), ifelse(is.null(Bfixed), FALSE, TRUE), 
               ifelse(is.null(Cfixed), FALSE, TRUE), ifelse(is.null(Dfixed), FALSE, TRUE))
    struc <- c(ifelse(is.null(Astruc), FALSE, TRUE), ifelse(is.null(Bstruc), FALSE, TRUE), 
               ifelse(is.null(Cstruc), FALSE, TRUE), ifelse(is.null(Dstruc), FALSE, TRUE))
    pfac <- list(A = Amat, B = Bmat, C = Cmat, D = Dmat, SSE = ssenew, Rsq = Rsq,
                 GCV = GCV, edf = edf, iter = iter, cflag = cflag, const = const,
                 control = control, fixed = fixed, struc = struc)
    return(pfac)
    
  }