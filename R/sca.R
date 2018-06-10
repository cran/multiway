sca <- 
  function(X, nfac, nstart = 10, maxit = 500,
           type = c("sca-p", "sca-pf2", "sca-ind", "sca-ecp"),
           rotation = c("none", "varimax", "promax"),
           ctol = 1e-4, parallel = FALSE, cl = NULL, verbose = TRUE){
    # Simultaneous Component Analysis
    # via alternating least squares (ALS) or closed-form solution
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 25, 2018
    
    # check 'X' input
    if(is.array(X)){
      xdim <- dim(X)
      lxdim <- length(xdim)
      if(lxdim!=3L){stop("Input 'X' must be 3-way array or list of 3-way arrays")}
      if(any(is.na(X)) | any(is.nan(X)) | any(is.infinite(X))){stop("Input 'X' cannot contain NA, NaN, or Inf values")}
      mylist <- vector("list",xdim[3])
      for(kk in 1:xdim[3]){mylist[[kk]] <- X[,,kk]}
      X <- mylist
      rm(mylist)
    } else if(is.list(X)){
      d1x <- dim(X[[1]])
      lxdim <- length(d1x) + 1L
      if( any(sapply(X,function(x) any(any(is.na(x)),any(is.nan(x)),any(is.infinite(x))))) ){stop("Input 'X' cannot contain NA, NaN, or Inf values")}
      if(lxdim==3L){
        xdim <- rep(NA,3)
        xdim[2] <- d1x[2]
        xdim[3] <- length(X)
        if(sum((sapply(X,ncol)-xdim[2])^2)>0L){stop("Input 'X' must be list of matrices with same number of columns.")}
      } else{stop("Input 'X' must be list of 2-way arrays.")}
    } else{stop("Input 'X' must be an array or list.")}
    xcx <- sumsq(X)
    
    # check 'nfac' and 'nstart' inputs
    nfac <- as.integer(nfac[1])
    if(nfac<1L){stop("Input 'nfac' must be positive integer")}
    nstart <- as.integer(nstart[1])
    if(nstart<1L){stop("Input 'nstart' must be positive integer")}
    
    # check 'maxit' and 'ctol' inputs
    maxit <- as.integer(maxit[1])
    if(maxit<1L){stop("Input 'maxit' must be positive integer")}
    ctol <- as.numeric(ctol[1])
    if(ctol<0L){stop("Input 'ctol' must be positive numeric")}
    
    # check 'type' and 'rotation'
    type <- type[1]
    if(is.na(match(type,c("sca-p","sca-pf2","sca-ind","sca-ecp",
                          "p","pf2","ind","ecp")))){stop("Input 'type' must be one of four specified options")}
    if(type=="p") type <- "sca-p"
    if(type=="pf2") type <- "sca-pf2"
    if(type=="ind") type <- "sca-ind"
    if(type=="ecp") type <- "sca-ecp"
    rotation <- rotation[1]
    if(is.na(match(rotation,c("none","varimax","promax")))){stop("Input 'rotation' must be one of three specified options")}
    
    # check 'parallel' and 'cl' inputs
    if(parallel && !any(class(cl)=="cluster")){
      stop("Input 'cl' must be cluster (created by makeCluster) when parallel=TRUE \n  See examples in documentation:  ?sca")
    }
    
    # determine type of model fit
    if(type[1] == "sca-p"){
      
      # fit model
      matdata <- X[[1]]
      for(kk in 2:xdim[3]) matdata <- rbind(matdata, X[[kk]])
      mysvd <- svd(matdata, nu = nfac, nv = nfac)
      Dmat <- mysvd$u %*% (diag(nfac)*mysvd$d[1:nfac] / sqrt(xdim[2]))
      Bmat <- mysvd$v * sqrt(xdim[2])
      ssr <- sum(rowSums(Dmat^2))
      
      # rotate solution
      if(rotation=="varimax"){
        Vrot <- varimax(Bmat)
        Bmat <- Bmat %*% Vrot$rotmat
        Dmat <- Dmat %*% Vrot$rotmat
      } else if(rotation=="promax"){
        Prot <- promax(Bmat)
        Bmat <- Bmat %*% Prot$rotmat
        Dmat <- Dmat %*% t(solve(Prot$rotmat))
      } 
      
      # get factor score std dev and sse
      Dmats <- vector("list", xdim[3])
      dfc <- sse <- 0
      Cmat <- matrix(0,xdim[3],nfac)
      for(kk in 1:xdim[3]){
        newdim <- dim(X[[kk]])[1]
        dinds <- 1:newdim + dfc
        Dmats[[kk]] <- as.matrix(Dmat[dinds,])
        Cmat[kk,] <- sqrt(colSums(Dmats[[kk]]^2) / newdim)
        sse <- sse + sum((X[[kk]] - tcrossprod(Dmats[[kk]], Bmat))^2)
        dfc <- dfc + newdim
      }
      rm(Dmat)
      
      # fit statistics
      Rsq <- 1 - sse / xcx
      iter <- 1
      cflag <- 0
      Phimat <- NULL
      ntotal <- nrow(matdata)
      Adf <- ntotal - xdim[3] * (nfac + 1) / 2
      Gdf <- xdim[3] * (nfac + 1) / 2
      Bdf <- xdim[2] - 1L
      Cdf <- 0
      edf <- nfac * c(Adf+Gdf,Bdf,Cdf)
      pxdim <- ntotal * xdim[2]
      GCV <- (sse / pxdim) / (1 - sum(edf) / pxdim)^2
      
    } else if(type[1] == "sca-pf2"){
      
      # fit model
      if(parallel){
        
        nstartlist <- vector("list", nstart)
        nstartlist[1:nstart] <- nfac
        scamod <- parLapply(cl = cl, X = nstartlist, fun = "parafac2_3way",
                            data = X, xcx = xcx, maxit = maxit, ctol = ctol)
        widx <- which.min(sapply(scamod, function(x) x$SSE))
        scamod <- scamod[[widx]]
        
      } else {
        
        if(verbose) pbar <- txtProgressBar(min = 0, max = nstart, style = 3)
        
        scamod <- parafac2_3way(data = X, nfac = nfac, xcx = xcx,
                                maxit = maxit, ctol = ctol)
        if(verbose) setTxtProgressBar(pbar, 1)
        if(nstart > 1L){
          for(j in 2:nstart){
            scanew <- parafac2_3way(data = X, nfac = nfac, xcx = xcx,
                                    maxit = maxit, ctol = ctol)
            if(scanew$SSE < scamod$SSE) scamod <- scanew
            if(verbose) setTxtProgressBar(pbar, j)
          } # end for(j in 2:nstart)
        } # end if(nstart > 1L)
        
        if(verbose) close(pbar)
        
      } # end if(parallel)
      
      # rescale parafac2 solution
      Bmat <- scamod$B
      Cmat <- matrix(0, xdim[3], nfac)
      Dmats <- vector("list", xdim[3])
      gsqrt <- sqrt(diag(scamod$Phi))
      for(kk in 1:xdim[3]){
        Dmats[[kk]] <- scamod$A[[kk]] %*% (diag(nfac)*scamod$C[kk,])
        Cmat[kk,] <- scamod$C[kk,] * gsqrt / sqrt(dim(scamod$A[[kk]])[1])
      }
      
      # fit statistics
      Rsq <- scamod$Rsq
      iter <- scamod$iter
      cflag <- scamod$cflag
      Smat <- (diag(nfac)/gsqrt)
      Phimat <- Smat %*% scamod$Phi %*% Smat
      GCV <- scamod$GCV
      edf <- scamod$edf
      
    } else if(type[1] == "sca-ind"){
      
      # fit model
      if(parallel){
        
        nstartlist <- vector("list", nstart)
        nstartlist[1:nstart] <- nfac
        scamod <- parLapply(cl = cl, X = nstartlist, fun = "parafac2_3way",
                            data = X, xcx = xcx, const = c("orthog", "uncons", "uncons"),
                            maxit = maxit, ctol = ctol)
        widx <- which.min(sapply(scamod, function(x) x$SSE))
        scamod <- scamod[[widx]]
        
      } else {
        
        if(verbose) pbar <- txtProgressBar(min = 0, max = nstart, style = 3)
        
        scamod <- parafac2_3way(data = X, nfac = nfac, xcx = xcx,
                                maxit = maxit, ctol = ctol,
                                const = c("orthog", "uncons", "uncons"))
        if(verbose) setTxtProgressBar(pbar, 1)
        if(nstart > 1L){
          for(j in 2:nstart){
            scanew <- parafac2_3way(data = X, nfac = nfac, xcx = xcx,
                                    maxit = maxit, ctol = ctol,
                                    const = c("orthog", "uncons", "uncons"))
            if(scanew$SSE < scamod$SSE) scamod <- scanew
            if(verbose) setTxtProgressBar(pbar, j)
          } # end for(j in 2:nstart)
        } # end if(nstart > 1L)
        
        if(verbose) close(pbar)
        
      } # end if(parallel)
      
      # rescale parafac2 solution
      dg <- sqrt(diag(scamod$Phi))
      Bmat <- scamod$B
      Cmat <- matrix(0,xdim[3],nfac)
      Dmats <- vector("list",xdim[3])
      for(kk in 1:xdim[3]){
        Dmats[[kk]] <- scamod$A[[kk]]%*%(diag(nfac)*scamod$C[kk,])
        Cmat[kk,] <- scamod$C[kk,]*(dg/sqrt(dim(scamod$A[[kk]])[1]))
      }
      
      # fit statistics
      Rsq <- scamod$Rsq
      iter <- scamod$iter
      cflag <- scamod$cflag
      Phimat <- diag(nfac)
      GCV <- scamod$GCV
      edf <- scamod$edf
      
    } else if(type[1] == "sca-ecp"){
      
      # make Cfixed
      nks <- rep(0, xdim[3])
      for(kk in 1:xdim[3]){nks[kk] <- sqrt(dim(X[[kk]])[1])}
      Cfixed <- matrix(nks,xdim[3],nfac)
      
      # fit model
      if(parallel){
        
        nstartlist <- vector("list", nstart)
        nstartlist[1:nstart] <- nfac
        scamod <- parLapply(cl = cl, X = nstartlist, fun = "parafac2_3way",
                            data = X, xcx = xcx, const = c("orthog","uncons","uncons"),
                            maxit = maxit, ctol = ctol, Cfixed = Cfixed)
        widx <- which.min(sapply(scamod,function(x) x$SSE))
        scamod <- scamod[[widx]]
        
      } else {
        
        if(verbose) pbar <- txtProgressBar(min = 0, max = nstart, style = 3)
        
        scamod <- parafac2_3way(data = X, nfac = nfac, xcx = xcx,
                                maxit = maxit, ctol = ctol,
                                const = c("orthog", "uncons", "uncons"))
        if(verbose) setTxtProgressBar(pbar, 1)
        if(nstart > 1L){
          for(j in 2:nstart){
            scanew <- parafac2_3way(data = X, nfac = nfac, xcx = xcx,
                                    maxit = maxit, ctol = ctol,
                                    const = c("orthog", "uncons", "uncons"),
                                    Cfixed = Cfixed)
            if(scanew$SSE < scamod$SSE) scamod <- scanew
            if(verbose) setTxtProgressBar(pbar, j)
          } # end for(j in 2:nstart)
        } # end if(nstart > 1L)
        
        if(verbose) close(pbar)
        
      } # end if(parallel)
      
      # order solution
      if(nfac > 1L){
        fordr <- order(colSums(scamod$B^2), decreasing = TRUE)
        scamod$A <- lapply(scamod$A, function(x) x[,fordr])
        scamod$B <- scamod$B[,fordr]
        scamod$C <- scamod$C[,fordr]
        scamod$Phi <- scamod$Phi[fordr,fordr]
      }
      
      # rotate solution
      Bmat <- scamod$B
      rotmat <- diag(nfac)
      if(rotation=="varimax"){
        rotmat <- varimax(Bmat)$rotmat
        Bmat <- Bmat %*% rotmat
      } else if(rotation=="promax"){
        rotmat <- promax(Bmat)$rotmat
        Bmat <- Bmat %*% rotmat
        rotmat <- t(solve(rotmat))
      }
      Dmats <- vector("list",xdim[3])
      for(kk in 1:xdim[3]){
        Dmats[[kk]] <- scamod$A[[kk]] %*% (diag(nfac)*scamod$C[kk,]) %*% rotmat
      }
      Cmat <- scamod$C / nks
      
      # fit statistics
      Rsq <- scamod$Rsq
      iter <- scamod$iter
      cflag <- scamod$cflag
      Phimat <- scamod$Phi
      ntotal <- sum(sapply(X,nrow))
      Adf <- ntotal - xdim[3] * (nfac + 1) / 2
      Gdf <- 0
      Bdf <- xdim[2]
      Cdf <- 0
      edf <- nfac * c(Adf+Gdf,Bdf,Cdf)
      pxdim <- ntotal * xdim[2]
      ssenew <- (1 - scamod$Rsq) * xcx
      GCV <- (ssenew / pxdim) / (1 - sum(edf) / pxdim)^2
    }
    names(edf) <- c("A", "B", "C")
    scafit <- list(D = Dmats, B = Bmat, C = Cmat, Phi = Phimat,
                   SSE = xcx * (1 - Rsq), Rsq = Rsq, GCV = GCV,
                   edf = edf, iter = iter, cflag = cflag,
                   type = type, rotation = rotation)
    class(scafit) <- "sca"
    return(scafit)
    
  }