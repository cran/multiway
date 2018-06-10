reorder.mcr <- 
  function(x, neworder,  mode = "A", ...){
    # Reorder Factors of fit MCR model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: February 2, 2018
    
    # check mode
    if(x$model == "tucker"){
      mode <- mode[1]
      if(!any(mode == c("A","B","C"))) stop("Incorrect input for 'mode'. Must set to 'A', 'B', or 'C' for MCR-Tucker")
      nfac <- ncol(x[[match(mode, LETTERS[1:3])]])
    } else {
      nfac <- ncol(x$B)
    }
    
    # check neworder
    neworder <- as.integer(neworder)
    if(length(neworder) != nfac) stop("Incorrect input for 'neworder'. Must have length equal to number of factors.")
    if(!identical(seq(1L,nfac), sort(neworder))) stop(paste("Incorrect input for 'neworder'. Must be unique integers in range of 1 to",nfac))
    
    # check model
    if(x$model == "parafac"){
      x$A <- x$A[,neworder]
      x$B <- x$B[,neworder]
      x$C <- x$C[,neworder]
      x$D <- x$D[,neworder]
      x$W <- x$W[,neworder]
      return(x)
    } else if(x$model == "parafac2"){
      for(k in 1:length(x$A)) x$A[[k]] <- x$A[[k]][,neworder]
      x$Phi <- x$Phi[neworder,neworder]
      x$B <- x$B[,neworder]
      x$C <- x$C[,neworder]
      x$D <- x$D[,neworder]
      x$W <- x$W[,neworder]
      return(x)
    } else {
      if(mode == "A"){
        x$A <- x$A[,neworder]
        x$G <- x$G[neworder,,]
      } else if(mode == "B"){
        x$B <- x$B[,neworder]
        x$G <- x$G[,neworder,]
      } else{
        x$C <- x$C[,neworder]
        x$D <- x$D[,neworder]
        x$W <- x$W[,neworder]
        x$G <- x$G[,,neworder]
      }
      return(x)
    }
    
  }