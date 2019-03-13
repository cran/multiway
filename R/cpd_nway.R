cpd_nway <-
  function(data, nfac, xcx = sum(data^2), maxit = 500, ctol = 1e-4){
    # N-way Canonical Polyadic Decomposition (CPD)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: March 11, 2019
    
    ### initializations
    xdims <- dim(data)
    pxdim <- prod(xdims)
    nmode <- length(xdims)
    index <- 1:nmode
    Amats <- vector("list", nmode)
    CPmat <- matrix(1, nfac, nfac)
    Amats[[1]] <- matrix(1 / sqrt(xdims[1]), xdims[1], nfac)
    for(n in 2:nmode) {
      Amats[[n]] <- matrix(rnorm(xdims[n] * nfac), xdims[n], nfac)
      CPmat <- CPmat * crossprod(Amats[[n]])
    }
    
    ### iterative update of matrices
    vtol <- sseold <- xcx + ctol
    iter <- 0
    cflag <- NA
    while((vtol > ctol) && (iter < maxit)) {
      
      ## loop through modes
      for(n in 1:nmode){
        
        # index minus n
        nindex <- index[-n]
        
        # matricized version of data
        Xmat <- matrix(aperm(data, perm = c(n, nindex)), 
                       nrow = xdims[n], ncol = pxdim / xdims[n])
        
        # Khatri-Rao product matrix
        Zmat <- krprod(Amats[[nindex[2]]], Amats[[nindex[1]]])
        if(nmode > 3L) {
          for(m in 3:(nmode - 1L)) {
            Zmat <- krprod(Amats[[nindex[m]]], Zmat)
          }
        }
        
        # update weights
        CPmat <- CPmat / crossprod(Amats[[n]])
        Amats[[n]] <- Xmat %*% Zmat %*% smpower(CPmat, power = -1)
        CPmat <- CPmat * crossprod(Amats[[n]])
        
      } # end for(n in 1:nmode)
      
      ## convergence check
      ssenew <- sum((Xmat - Amats[[n]] %*% t(Zmat))^2)
      vtol <- (sseold - ssenew) / xcx
      sseold <- ssenew
      iter <- iter + 1
      
    } # end while((vtol > ctol) && (iter < maxit))
    
    ### put scale in N-th mode
    svec <- rep(1, nfac)
    for(n in 1:(nmode-1)){
      nvec <- sqrt(colMeans(Amats[[n]]^2))
      Amats[[n]] <- scale(Amats[[n]], center = FALSE, scale = nvec)
      svec <- svec * nvec
    }
    Amats[[nmode]] <- Amats[[nmode]] %*% diag(svec)
    
    ### order solution
    fordr <- order(colSums(Amats[[nmode]]^2), decreasing = TRUE)
    for(n in 1:nmode) Amats[[n]] <- Amats[[n]][,fordr]
    
    ### GCV and R-squared
    edf <- xdims * nfac - rep(c(nfac, 0), c(nmode - 1, 1))
    GCV <- (ssenew / pxdim) / (1 - sum(edf) / pxdim)^2
    Rsq <- 1 - ssenew / xcx
    
    ### collect results
    if(is.na(cflag)) cflag <- ifelse(vtol <= ctol, 0, 1)
    cpd <- list(A = Amats, SSE = ssenew, Rsq = Rsq, 
                GCV = GCV, edf = edf, iter = iter, cflag = cflag)
    return(cpd)
    
  } # end cdp_nway