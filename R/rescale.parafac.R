rescale.parafac <-
  function(x, mode="A", newscale=1, absorb="C", ...){
    # Rescales Weights of fit Parafac model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 25, 2018
    
    # check mode and absorb
    mode <- mode[1]
    absorb <- absorb[1]
    if(mode==absorb) stop("Inputs 'mode' and 'absorb' must be different.")
    if(is.null(x$D)){
      if(!any(mode==c("A","B","C"))) stop("Incorrect input for 'mode'. Must set to 'A', 'B', or 'C' for 3-way Parafac.")
      if(!any(absorb==c("A","B","C"))) stop("Incorrect input for 'absorb'. Must set to 'A', 'B', or 'C' for 3-way Parafac.")
    } else {
      if(!any(mode==c("A","B","C","D"))) stop("Incorrect input for 'mode'. Must set to 'A', 'B', 'C', or 'D' for 4-way Parafac.")
      if(!any(absorb==c("A","B","C","D"))) stop("Incorrect input for 'absorb'. Must set to 'A', 'B', 'C', or 'D' for 4-way Parafac.")
    }
    
    # check newscale
    nfac <- ncol(x$A)
    if(length(newscale)!=nfac) newscale <- rep(newscale[1],nfac)
    if(any(newscale <= 0)) stop("Input 'newscale' must contain positive values.")
    
    # rescale factors
    if(mode=="A"){
      
      Ascale <- sqrt(colMeans(x$A^2))
      if(any(Ascale == 0)) Ascale[Ascale == 0] <- 1
      svec <- newscale/Ascale
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$A <- x$A %*% Smat
      if(nfac==1L) { Smat <- matrix(1/svec) } else { Smat <- diag(1/svec) }
      if(absorb=="B") {
        x$B <- x$B %*% Smat
      } else if(absorb=="C"){
        x$C <- x$C %*% Smat
      } else {
        x$D <- x$D %*% Smat
      }
      return(x)
      
    } else if(mode=="B"){
      
      Bscale <- sqrt(colMeans(x$B^2))
      if(any(Bscale == 0)) Bscale[Bscale == 0] <- 1
      svec <- newscale/Bscale
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$B <- x$B %*% Smat
      if(nfac==1L) { Smat <- matrix(1/svec) } else { Smat <- diag(1/svec) }
      if(absorb=="A") {
        x$A <- x$A %*% Smat
      } else if(absorb=="C"){
        x$C <- x$C %*% Smat
      } else {
        x$D <- x$D %*% Smat
      }
      return(x)
      
    } else if(mode=="C"){
      
      Cscale <- sqrt(colMeans(x$C^2))
      if(any(Cscale == 0)) Cscale[Cscale == 0] <- 1
      svec <- newscale/Cscale
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$C <- x$C %*% Smat
      if(nfac==1L) { Smat <- matrix(1/svec) } else { Smat <- diag(1/svec) }
      if(absorb=="A") {
        x$A <- x$A %*% Smat
      } else if(absorb=="B"){
        x$B <- x$B %*% Smat
      } else {
        x$D <- x$D %*% Smat
      }
      return(x)
      
    } else if(mode=="D"){
      
      Dscale <- sqrt(colMeans(x$D^2))
      if(any(Dscale == 0)) Dscale[Dscale == 0] <- 1
      svec <- newscale/Dscale
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$D <- x$D %*% Smat
      if(nfac==1L) { Smat <- matrix(1/svec) } else { Smat <- diag(1/svec) }
      if(absorb=="A") {
        x$A <- x$A %*% Smat
      } else if(absorb=="B"){
        x$B <- x$B %*% Smat
      } else {
        x$C <- x$C %*% Smat
      }
      return(x)
      
    }
    
  }