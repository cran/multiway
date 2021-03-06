resign.tucker <-
  function(x, mode="A", newsign=1, ...){
    # Resigns Weights of fit Tucker model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 25, 2018
    
    # check mode and dimensions
    mode <- mode[1]
    mydim <- c(ncol(x$A),ncol(x$B),ncol(x$C))
    apdim <- 1:3
    if(is.null(x$D)){
      if(!any(mode==c("A","B","C"))) stop("Incorrect input for 'mode'. Must set to 'A', 'B', or 'C' for 3-way Tucker")
    } else {
      mydim <- c(mydim,ncol(x$D))
      apdim <- c(apdim,4)
      if(!any(mode==c("A","B","C","D"))) stop("Incorrect input for 'mode'. Must set to 'A', 'B', 'C', or 'D' for 4-way Parafac.")
    }
    
    # resign factors
    newsign <- sign(newsign)
    if(any(newsign == 0)) stop("Input 'newsign' must contain entries of c(-1, 1).")
    if(mode=="A"){
      
      if(length(newsign)!=mydim[1]) newsign <- rep(newsign[1],mydim[1])
      Asign <- sign(colMeans(x$A^3))
      if(any(Asign == 0)) Asign[Asign == 0] <- 1
      svec <- newsign*Asign
      if(mydim[1]==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$A <- x$A %*% Smat
      Gmat <- matrix(x$G, mydim[1], prod(mydim[-1]))
      x$G <- array(Smat %*% Gmat, dim=mydim)
      return(x)
      
    } else if(mode=="B"){
      
      if(length(newsign)!=mydim[2]) newsign <- rep(newsign[1],mydim[2])
      permvec <- c(apdim[2],apdim[-2])
      Bsign <- sign(colMeans(x$B^3))
      if(any(Bsign == 0)) Bsign[Bsign == 0] <- 1
      svec <- newsign*Bsign
      if(mydim[2]==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$B <- x$B %*% Smat
      Gmat <- matrix(aperm(x$G, permvec), mydim[2], prod(mydim[-2]))
      x$G <- aperm(array(Smat %*% Gmat, dim=c(mydim[2],mydim[-2])), sort(permvec,index=T)$ix)
      return(x)
      
    } else if(mode=="C"){
      
      if(length(newsign)!=mydim[3]) newsign <- rep(newsign[1],mydim[3])
      permvec <- c(apdim[3],apdim[-3])
      Csign <- sign(colMeans(x$C^3))
      if(any(Csign == 0)) Csign[Csign == 0] <- 1
      svec <- newsign*Csign
      if(mydim[3]==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$C <- x$C %*% Smat
      Gmat <- matrix(aperm(x$G, permvec), mydim[3], prod(mydim[-3]))
      x$G <- aperm(array(Smat %*% Gmat, dim=c(mydim[3],mydim[-3])), sort(permvec,index=T)$ix)
      return(x)
      
    } else if(mode=="D"){
      
      if(length(newsign)!=mydim[4]) newsign <- rep(newsign[1],mydim[4])
      permvec <- c(apdim[4],apdim[-4])
      Dsign <- sign(colMeans(x$D^3))
      if(any(Dsign == 0)) Dsign[Dsign == 0] <- 1
      svec <- newsign*Dsign
      if(mydim[4]==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$D <- x$D %*% Smat
      Gmat <- matrix(aperm(x$G, permvec), mydim[4], prod(mydim[-4]))
      x$G <- aperm(array(Smat %*% Gmat, dim=c(mydim[4],mydim[-4])), sort(permvec,index=T)$ix)
      return(x)
      
    }
    
  }