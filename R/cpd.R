cpd <- 
  function(X, nfac, nstart = 10, maxit = 500, 
           ctol = 1e-4, parallel = FALSE, cl = NULL, 
           output = "best", verbose = TRUE){
    # N-way Canonical Polyadic Decomposition (CPD)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: March 11, 2019
    
    # check 'X' input
    xdim <- dim(X)
    lxdim <- length(xdim)
    if(any(is.nan(X)) | any(is.infinite(X))) stop("Input 'X' cannot contain NaN or Inf values") 
    if(any(is.na(X))){
      missingdata <- TRUE
      naid <- which(is.na(X))
    } else {
      missingdata <- FALSE
      xcx <- sumsq(X)
    }
    
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
    
    # check 'parallel' and 'cl' inputs
    if(parallel && !any(class(cl)=="cluster")){
      stop("Input 'cl' must be cluster (created by makeCluster) when parallel=TRUE \n  See examples in documentation:  ?parafac")
    }
    
    # cpd fitting
    if(parallel){
      nstartlist <- vector("list", nstart)
      nstartlist[1:nstart] <- nfac
      if(missingdata){
        pfaclist <- parLapply(cl = cl, X = nstartlist, fun = "cpd_nwayna",
                              data = X, naid = naid, maxit = maxit, ctol = ctol)
      } else {
        pfaclist <- parLapply(cl = cl, X = nstartlist, fun = "cpd_nway",
                              data = X, xcx = xcx, maxit = maxit, ctol = ctol)
      }
    } else {
      
      if(verbose) pbar <- txtProgressBar(min = 0, max = nstart, style = 3)
      
      if(output[1] == "best"){
        if(missingdata){
          pfac <- cpd_nwayna(data = X, nfac = nfac, naid = naid, 
                             maxit = maxit, ctol = ctol)
          if(verbose) setTxtProgressBar(pbar, 1)
          if(nstart > 1L){
            for(j in 2:nstart){
              pnew <- cpd_nwayna(data = X, nfac = nfac, naid = naid, 
                                 maxit = maxit, ctol = ctol)
              if(pnew$SSE < pfac$SSE) pfac <- pnew
              if(verbose) setTxtProgressBar(pbar, j)
            } # end for(j in 2:nstart)
          } # end if(nstart > 1L)
        } else {
          pfac <- cpd_nway(data = X, nfac = nfac, xcx = xcx,
                           maxit = maxit, ctol = ctol)
          if(verbose) setTxtProgressBar(pbar, 1)
          if(nstart > 1L){
            for(j in 2:nstart){
              pnew <- cpd_nway(data = X, nfac = nfac, xcx = xcx, 
                               maxit = maxit, ctol = ctol)
              if(pnew$SSE < pfac$SSE) pfac <- pnew
              if(verbose) setTxtProgressBar(pbar, j)
            } # end for(j in 2:nstart)
          } # end if(nstart > 1L)
        } # end if(missingdata)
      } else {
        
        pfaclist <- vector("list", nstart)
        if(missingdata){
          for(j in 1:nstart){
            pfaclist[[j]] <- cpd_nwayna(data = X, nfac = nfac, naid = naid, 
                                        maxit = maxit, ctol = ctol)
            if(verbose) setTxtProgressBar(pbar, j)
          }
        } else {
          for(j in 1:nstart){
            pfaclist[[j]] <- cpd_nway(data = X, nfac = nfac, xcx = xcx, 
                                      maxit = maxit, ctol = ctol)
            if(verbose) setTxtProgressBar(pbar, j)
          }
        } # end if(missingdata)
        
      } # end if(output[1] == "best")
      
      if(verbose) close(pbar)
      
    } # end if(parallel)
    
    # output results
    if(output[1] == "best"){
      if(parallel){
        SSE <- sapply(pfaclist,function(x) x$SSE)
        pfac <- pfaclist[[which.min(SSE)]]
      }
      class(pfac) <- "cpd"
      return(pfac)
    } else {
      pfaclist <- lapply(pfaclist, function(pfac) {
        class(pfac) <- "cpd"
        pfac
      })
      return(pfaclist)
    }
    
  } # end cpd