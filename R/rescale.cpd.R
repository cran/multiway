rescale.cpd <-
  function(x, mode=1, newscale=1, absorb=3, ...){
    # Rescales Weights of fit CPD model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: March 12, 2019
    
    # check mode and absorb
    mode <- as.integer(mode[1])
    absorb <- as.integer(absorb[1])
    if(mode == absorb) stop("Inputs 'mode' and 'absorb' must be different.")
    nmodes <- length(x$A)
    if(mode < 1 | mode > nmodes) stop("Incorrect input for 'mode'. Must be an integer between 1 and",nmodes)
    if(absorb < 1 | absorb > nmodes) stop("Incorrect input for 'absorb'. Must be an integer between 1 and",nmodes)
    
    # check newscale
    nfac <- ncol(x$A[[1]])
    if(length(newscale)!=nfac) newscale <- rep(newscale[1],nfac)
    if(any(newscale <= 0)) stop("Input 'newscale' must contain positive values.")
    
    # rescale factors
    Ascale <- sqrt(colMeans(x$A[[mode]]^2))
    if(any(Ascale == 0)) Ascale[Ascale == 0] <- 1
    svec <- newscale / Ascale
    if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
    x$A[[mode]] <- x$A[[mode]] %*% Smat
    if(nfac==1L) { Smat <- matrix(1/svec) } else { Smat <- diag(1/svec) }
    x$A[[absorb]] <- x$A[[absorb]] %*% Smat
    return(x)
    
  }