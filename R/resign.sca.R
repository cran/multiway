resign.sca <-
  function(x, mode="B", newsign=1, ...){
    # Resigns Weights of fit SCA nmodel
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 25, 2018
    
    # check mode
    mode <- mode[1]
    if(mode=="B"){
      absorb <- "C"
    } else if(mode=="C"){
      absorb <- "B"
    } else {
      stop("Incorrect input for 'mode'. Must set to 'B' or 'C' for SCA")
    }
    
    # check newsign
    nfac <- ncol(x$B)
    newsign <- sign(newsign)
    if(length(newsign)!=nfac) newsign <- rep(newsign[1],nfac)
    if(any(newsign == 0)) stop("Input 'newsign' must contain entries of c(-1, 1).")
    
    # resign factors
    if(mode=="B"){
      
      Bsign <- sign(colMeans(x$B^3))
      if(any(Bsign == 0)) Bsign[Bsign == 0] <- 1
      svec <- newsign*Bsign
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$B <- x$B %*% Smat
      x$C <- x$C %*% Smat
      x$D <- lapply(x$D, function(x) x %*% Smat)
      return(x)
      
    } else {
      
      Csign <- sign(colMeans(x$C^3))
      if(any(Csign == 0)) Csign[Csign == 0] <- 1
      svec <- newsign*Csign
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$C <- x$C %*% Smat
      x$D <- lapply(x$D, function(x) x %*% Smat)
      x$B <- x$B %*% Smat
      return(x)
      
    }
    
  }