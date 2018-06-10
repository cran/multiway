rescale.mcr <-
  function(x, mode = "A", newscale = 1, absorb = "C", ...){
    # Rescales Weights of fit MCR model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 25, 2018
    
    ### check mode and absorb
    mode <- mode[1]
    absorb <- absorb[1]
    if(mode == absorb) stop("Inputs 'mode' and 'absorb' must be different.")
    if(!any(mode == c("A","B","C"))) stop("Incorrect input for 'mode'. Must set to 'A', 'B', or 'C'.")
    if(!any(absorb == c("A","B","C"))) stop("Incorrect input for 'absorb'. Must set to 'A', 'B', or 'C'.")
    
    ### check newscale
    if(x$model == "tucker"){
      nfac <- ncol(x[[match(mode, LETTERS[1:3])]])
    } else {
      nfac <- ncol(x$B)
    }
    if(length(newscale) != nfac) newscale <- rep(newscale[1], nfac)
    if(any(newscale <= 0)) stop("Input 'newscale' must contain positive values.")
    
    ### MCR-Parafac
    if(x$model == "parafac"){
      
      # rescale mode
      if(mode == "A"){
        Ascale <- sqrt(colMeans(x$A^2))
        if(any(Ascale == 0)) Ascale[Ascale == 0] <- 1
        svec <- newscale / Ascale
        if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        x$A <- x$A %*% Smat
      } else if(mode == "B"){
        Bscale <- sqrt(colMeans(x$B^2))
        if(any(Bscale == 0)) Bscale[Bscale == 0] <- 1
        svec <- newscale / Bscale
        if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        x$B <- x$B %*% Smat
      } else {
        Cscale <- sqrt(colMeans(x$C^2))
        if(any(Cscale == 0)) Cscale[Cscale == 0] <- 1
        svec <- newscale / Cscale
        if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        x$C <- x$C %*% Smat
        x$W <- x$W %*% Smat
        if(nfac==1L) { Sinv <- matrix(1/svec) } else { Sinv <- diag(1/svec) }
        x$D <- x$D %*% Sinv
      } # end if(mode == "A")
      
      # rescale absorb
      if(nfac==1L) { Sinv <- matrix(1/svec) } else { Sinv <- diag(1/svec) }
      if(absorb == "A"){
        x$A <- x$A %*% Sinv
      } else if(absorb == "B"){
        x$B <- x$B %*% Sinv
      } else {
        x$C <- x$C %*% Sinv
        x$W <- x$W %*% Sinv
        x$D <- x$D %*% Smat
      } # end if(absorb == "A")
      
      return(x)
      
    } # end if(x$model == "parafac")
    
    ### MCR-Parafac2
    if(x$model == "parafac2"){
      
      # rescale mode
      if(mode == "A"){
        Ascale <- sqrt(diag(x$Phi) / mean(sapply(x$A,nrow)))
        if(any(Ascale == 0)) Ascale[Ascale == 0] <- 1
        svec <- newscale / Ascale
        if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        for(k in 1:length(x$A)) x$A[[k]] <- x$A[[k]] %*% Smat
        x$Phi <- Smat %*% x$Phi %*% Smat
      } else if(mode == "B"){
        Bscale <- sqrt(colMeans(x$B^2))
        if(any(Bscale == 0)) Bscale[Bscale == 0] <- 1
        svec <- newscale / Bscale
        if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        x$B <- x$B %*% Smat
      } else {
        Cscale <- sqrt(colMeans(x$C^2))
        if(any(Cscale == 0)) Cscale[Cscale == 0] <- 1
        svec <- newscale / Cscale
        if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        x$C <- x$C %*% Smat
        x$W <- x$W %*% Smat
        if(nfac==1L) { Sinv <- matrix(1/svec) } else { Sinv <- diag(1/svec) }
        x$D <- x$D %*% Sinv
      } # end if(mode == "A")
      
      # rescale absorb
      if(nfac==1L) { Sinv <- matrix(1/svec) } else { Sinv <- diag(1/svec) }
      if(absorb == "A"){
        for(k in 1:length(x$A)) x$A[[k]] <- x$A[[k]] %*% Sinv
        x$Phi <- Sinv %*% x$Phi %*% Sinv
      } else if(absorb == "B"){
        x$B <- x$B %*% Sinv
      } else {
        x$C <- x$C %*% Sinv
        x$W <- x$W %*% Sinv
        x$D <- x$D %*% Smat
      } # end if(absorb == "A")
      
      return(x)
      
    } # end if(x$model == "parafac2")
    
    
    ### MCR-Tucker
    if(x$model == "tucker"){
      
      mydim <- c(ncol(x$A), ncol(x$B), ncol(x$C))
      apdim <- 1:3
      if(mode == "A"){
        Ascale <- sqrt(colMeans(x$A^2))
        if(any(Ascale == 0)) Ascale[Ascale == 0] <- 1
        svec <- newscale / Ascale
        if(mydim[1]==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        x$A <- x$A %*% Smat
        if(mydim[1]==1L) { Sinv <- matrix(1/svec) } else { Sinv <- diag(1/svec) }
        Gmat <- matrix(x$G, mydim[1], prod(mydim[-1]))
        x$G <- array(Sinv %*% Gmat, dim=mydim)
      } else if(mode == "B"){
        permvec <- c(apdim[2],apdim[-2])
        Bscale <- sqrt(colMeans(x$B^2))
        if(any(Bscale == 0)) Bscale[Bscale == 0] <- 1
        svec <- newscale / Bscale
        if(mydim[2]==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        x$B <- x$B %*% Smat
        if(mydim[2]==1L) { Sinv <- matrix(1/svec) } else { Sinv <- diag(1/svec) }
        Gmat <- matrix(aperm(x$G, permvec), mydim[2], prod(mydim[-2]))
        x$G <- aperm(array(Sinv %*% Gmat, dim=c(mydim[2],mydim[-2])), sort(permvec,index=T)$ix)
      } else {
        permvec <- c(apdim[3],apdim[-3])
        Cscale <- sqrt(colMeans(x$C^2))
        if(any(Cscale == 0)) Cscale[Cscale == 0] <- 1
        svec <- newscale / Cscale
        if(mydim[3]==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        x$C <- x$C %*% Smat
        x$W <- x$W %*% Smat
        if(mydim[3]==1L) { Sinv <- matrix(1/svec) } else { Sinv <- diag(1/svec) }
        x$D <- x$D %*% Sinv
        Gmat <- matrix(aperm(x$G, permvec), mydim[3], prod(mydim[-3]))
        x$G <- aperm(array(Sinv %*% Gmat, dim=c(mydim[3],mydim[-3])), sort(permvec,index=T)$ix)
      } # end if(mode == "A")
      
      return(x)
      
    } # end if(x$model == "tucker")
    
  }