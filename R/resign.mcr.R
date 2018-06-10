resign.mcr <-
  function(x, mode = "A", newsign = 1, absorb = "C", ...){
    # Resigns Weights of fit MCR model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 25, 2018
    
    ### check mode and absorb
    mode <- mode[1]
    absorb <- absorb[1]
    if(mode == absorb) stop("Inputs 'mode' and 'absorb' must be different.")
    if(!any(mode == c("A","B","C"))) stop("Incorrect input for 'mode'. Must set to 'A', 'B', or 'C'.")
    if(!any(absorb == c("A","B","C"))) stop("Incorrect input for 'absorb'. Must set to 'A', 'B', or 'C'.")
    
    ### check newsign
    if(x$model == "tucker"){
      nfac <- ncol(x[[match(mode, LETTERS[1:3])]])
    } else {
      nfac <- ncol(x$B)
    }
    newsign <- sign(newsign)
    if(any(newsign == 0)) stop("Input 'newsign' must contain entries of c(-1, 1).")
    
    ### MCR-Parafac
    if(x$model == "parafac"){
      
      if(length(newsign) != nfac) newsign <- rep(newsign[1], nfac)
      
      # resign factors
      if(mode == "A"){
        Asign <- sign(colMeans(x$A^3))
        if(any(Asign == 0)) Asign[Asign == 0] <- 1
        svec <- newsign*Asign
        if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        x$A <- x$A %*% Smat
        if(absorb == "B"){
          x$B <- x$B %*% Smat
        } else {
          x$C <- x$C %*% Smat
          x$D <- x$D %*% Smat
          x$W <- x$W %*% Smat
        }
      } else if(mode == "B"){
        Bsign <- sign(colMeans(x$B^3))
        if(any(Bsign == 0)) Bsign[Bsign == 0] <- 1
        svec <- newsign*Bsign
        if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        x$B <- x$B %*% Smat
        if(absorb == "A"){
          x$A <- x$A %*% Smat
        } else {
          x$C <- x$C %*% Smat
          x$D <- x$D %*% Smat
          x$W <- x$W %*% Smat
        }
      } else {
        Csign <- sign(colMeans(x$C^3))
        if(any(Csign == 0)) Csign[Csign == 0] <- 1
        svec <- newsign*Csign
        if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        x$C <- x$C %*% Smat
        x$D <- x$D %*% Smat
        x$W <- x$W %*% Smat
        if(absorb == "A"){
          x$A <- x$A %*% Smat
        } else {
          x$B <- x$B %*% Smat
        }
      } # end if(mode == "A")
      
      return(x)
      
    } # end if(x$model == "parafac")
    
    ### MCR-Parafac2
    if(x$model == "parafac2"){
      
      if(length(newsign) != nfac) newsign <- rep(newsign[1], nfac)
      
      # resign factors
      if(mode == "A"){
        Asign <- rep(0, nfac)
        for(k in 1:length(x$A)) Asign <- Asign + colMeans(x$A[[k]]^3)
        Asign <- sign(Asign)
        if(any(Asign == 0)) Asign[Asign == 0] <- 1
        svec <- newsign*Asign
        if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        for(k in 1:length(x$A)) x$A[[k]] <- x$A[[k]] %*% Smat
        x$Phi <- Smat %*% x$Phi %*% Smat
        if(absorb=="B") {
          x$B <- x$B %*% Smat
        } else {
          x$C <- x$C %*% Smat
          x$D <- x$D %*% Smat
          x$W <- x$W %*% Smat
        }
      } else if(mode == "B"){
        Bsign <- sign(colMeans(x$B^3))
        if(any(Bsign == 0)) Bsign[Bsign == 0] <- 1
        svec <- newsign*Bsign
        if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        x$B <- x$B %*% Smat
        if(absorb == "A"){
          for(k in 1:length(x$A)) x$A[[k]] <- x$A[[k]] %*% Smat
          x$Phi <- Smat %*% x$Phi %*% Smat
        } else {
          x$C <- x$C %*% Smat
          x$D <- x$D %*% Smat
          x$W <- x$W %*% Smat
        }
      } else {
        Csign <- sign(colMeans(x$C^3))
        if(any(Csign == 0)) Csign[Csign == 0] <- 1
        svec <- newsign*Csign
        if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        x$C <- x$C %*% Smat
        x$D <- x$D %*% Smat
        x$W <- x$W %*% Smat
        if(absorb=="A") {
          for(k in 1:length(x$A)) x$A[[k]] <- x$A[[k]] %*% Smat
          x$Phi <- Smat %*% x$Phi %*% Smat
        } else {
          x$B <- x$B %*% Smat
        }
      } # end if(mode == "A")
      
      return(x)
      
    } # end if(x$model == "parafac2")
    
    
    ### MCR-Tucker
    if(x$model == "tucker"){
      
      mydim <- c(ncol(x$A), ncol(x$B), ncol(x$C))
      apdim <- 1:3
      if(mode == "A"){
        if(length(newsign) != mydim[1]) newsign <- rep(newsign[1], mydim[1])
        Asign <- sign(colMeans(x$A^3))
        if(any(Asign == 0)) Asign[Asign == 0] <- 1
        svec <- newsign*Asign
        if(mydim[1]==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        x$A <- x$A %*% Smat
        Gmat <- matrix(x$G, mydim[1], prod(mydim[-1]))
        x$G <- array(Smat %*% Gmat, dim = mydim)
      } else if(mode == "B"){
        if(length(newsign) != mydim[2]) newsign <- rep(newsign[1], mydim[2])
        permvec <- c(apdim[2], apdim[-2])
        Bsign <- sign(colMeans(x$B^3))
        if(any(Bsign == 0)) Bsign[Bsign == 0] <- 1
        svec <- newsign*Bsign
        if(mydim[2]==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        x$B <- x$B %*% Smat
        Gmat <- matrix(aperm(x$G, permvec), mydim[2], prod(mydim[-2]))
        x$G <- aperm(array(Smat %*% Gmat, dim = c(mydim[2], mydim[-2])), sort(permvec, index = TRUE)$ix)
      } else {
        if(length(newsign) != mydim[3]) newsign <- rep(newsign[1], mydim[3])
        permvec <- c(apdim[3], apdim[-3])
        Csign <- sign(colMeans(x$C^3))
        if(any(Csign == 0)) Csign[Csign == 0] <- 1
        svec <- newsign*Csign
        if(mydim[3]==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
        x$C <- x$C %*% Smat
        x$D <- x$D %*% Smat
        x$W <- x$W %*% Smat
        Gmat <- matrix(aperm(x$G, permvec), mydim[3], prod(mydim[-3]))
        x$G <- aperm(array(Smat %*% Gmat, dim = c(mydim[3], mydim[-3])), sort(permvec, index = TRUE)$ix)
      } # end if(mode == "A")
      
      return(x)
      
    } # end if(x$model == "tucker")
    
  }