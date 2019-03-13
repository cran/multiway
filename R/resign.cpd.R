resign.cpd <-
  function(x, mode=1, newsign=1, absorb=3, ...){
    # Resigns Weights of fit CPD model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: March 12, 2019
    
    # check mode and absorb
    mode <- as.integer(mode[1])
    absorb <- as.integer(absorb[1])
    if(mode == absorb) stop("Inputs 'mode' and 'absorb' must be different.")
    nmodes <- length(x$A)
    if(mode < 1 | mode > nmodes) stop("Incorrect input for 'mode'. Must be an integer between 1 and",nmodes)
    if(absorb < 1 | absorb > nmodes) stop("Incorrect input for 'absorb'. Must be an integer between 1 and",nmodes)
    
    # check newsign
    nfac <- ncol(x$A[[1]])
    newsign <- sign(newsign)
    if(length(newsign)!=nfac) newsign <- rep(newsign[1],nfac)
    if(any(newsign == 0)) stop("Input 'newsign' must contain entries of c(-1, 1).")
    
    # resign factors
    Asign <- sign(colMeans(x$A[[mode]]^3))
    if(any(Asign == 0)) Asign[Asign == 0] <- 1
    svec <- newsign*Asign
    if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
    x$A[[mode]] <- x$A[[mode]] %*% Smat
    x$A[[absorb]] <- x$A[[absorb]] %*% Smat
    return(x)
    
  }