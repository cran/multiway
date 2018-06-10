resign.indscal <-
  function(x, mode="B", newsign=1, ...){
    # Resigns Weights of fit INDSCAL model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 25, 2018
    
    # check newsign
    nfac <- ncol(x$B)
    newsign <- sign(newsign)
    if(length(newsign)!=nfac) newsign <- rep(newsign[1],nfac)
    if(any(newsign == 0)) stop("Input 'newsign' must contain entires of c(-1, 1).")
    
    # resign factors
    Bsign <- sign(colMeans(x$B^3))
    if(any(Bsign == 0)) Bsign[Bsign == 0] <- 1
    svec <- newsign*Bsign
    if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
    x$B <- x$B %*% Smat
    return(x)
    
  }