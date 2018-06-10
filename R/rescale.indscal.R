rescale.indscal <-
  function(x, mode = "B", newscale = 1, ...){
    # Rescales Weights of fit INDSCAL model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 25, 2018
    
    # check mode
    mode <- mode[1]
    if(!any(mode==c("B","C"))) stop("Incorrect input for 'mode'. Must set to 'B' or 'C' for INDSCAL.")
    
    # check newscale
    nfac <- ncol(x$B)
    if(length(newscale)!=nfac) newscale <- rep(newscale[1],nfac)
    if(any(newscale <= 0)) stop("Input 'newscale' must contain positive values.")
    
    # rescale factors
    if(mode=="B"){
      
      Bscale <- sqrt(colMeans(x$B^2))
      if(any(Bscale == 0)) Bscale[Bscale == 0] <- 1
      svec <- newscale/Bscale
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$B <- x$B %*% Smat
      if(nfac==1L) { Smat <- matrix(1/svec^2) } else { Smat <- diag(1/svec^2) }
      x$C <- x$C %*% Smat
      return(x)
      
    } else {
      
      Cscale <- sqrt(colMeans(x$C^2))
      if(any(Cscale == 0)) Cscale[Cscale == 0] <- 1
      svec <- newscale/Cscale
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$C <- x$C %*% Smat
      if(nfac==1L) { Smat <- matrix(1/svec^0.5) } else { Smat <- diag(1/svec^0.5) }
      x$B <- x$B %*% Smat
      return(x)
      
    } 
    
  }