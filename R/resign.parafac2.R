resign.parafac2 <-
  function(x, mode="A", newsign=1, absorb="C", method="pearson", ...){
    # Resigns Weights of fit Parafac2 model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 25, 2018
    
    # check mode and absorb and method
    mode <- mode[1]
    absorb <- absorb[1]
    method <- method[1]
    if(mode==absorb) stop("Inputs 'mode' and 'absorb' must be different.")
    if(is.null(x$D)){
      if(!any(mode==c("A","B","C"))) stop("Incorrect input for 'mode'. Must set to 'A', 'B', or 'C' for 3-way Parafac2.")
      if(!any(absorb==c("A","B","C"))) stop("Incorrect input for 'absorb'. Must set to 'A', 'B', or 'C' for 3-way Parafac2.")
    } else {
      if(!any(mode==c("A","B","C","D"))) stop("Incorrect input for 'mode'. Must set to 'A', 'B', 'C', or 'D' for 4-way Parafac2.")
      if(!any(absorb==c("A","B","C","D"))) stop("Incorrect input for 'absorb'. Must set to 'A', 'B', 'C', or 'D' for 4-way Parafac2.")
    }
    if(!any(method==c("pearson", "kendall", "spearman"))) stop("Incorrect input for 'method'. Must be set to 'pearson', 'kendall', or 'spearman'.")
    
    # check newsign
    nfac <- ncol(x$B)
    if(mode=="A" & is.list(newsign)){
      nmat <- ifelse(is.null(x$D),nrow(x$C),nrow(x$D))
      if(length(newsign)!=nmat) stop(paste("Input 'newsign' must be a list of length",nmat,"where each element contains a vector of covariates for resigning Mode A."))
      signlen <- sapply(newsign, function(x) length(c(x)))
      if(!identical(signlen, sapply(x$A, nrow))) stop("Incorrect length for an element of 'newsign'. Need length(newsign[[k]]) = nrow(x$A[[k]]) for all k.")
    } else {
      newsign <- sign(newsign)
      if(length(newsign)!=nfac) newsign <- rep(newsign[1],nfac)
      if(any(newsign == 0)) stop("Input 'newsign' must contain entries of c(-1, 1).")
    }
    
    # resign factors
    if(mode=="A"){
      
      if(is.list(newsign)){
        
        ksign <- rep(0,nmat)
        for(k in 1:nmat){
          ksign[k] <- sign( cor( (x$A[[k]][,1L]), c(newsign[[k]]), method=method ) )
          if(ksign[k] == 0) ksign[k] <- 1
          x$A[[k]] <- ksign[k] * x$A[[k]]
        }
        if(is.null(x$D)) { x$C <- diag(ksign) %*% x$C } else { x$D <- diag(ksign) %*% x$D }
        
      } else {
        
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
        } else if(absorb=="C"){
          x$C <- x$C %*% Smat
        } else {
          x$D <- x$D %*% Smat
        }
        
      }
      return(x)
      
    } else if(mode=="B"){
      
      Bsign <- sign(colMeans(x$B^3))
      if(any(Bsign == 0)) Bsign[Bsign == 0] <- 1
      svec <- newsign*Bsign
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$B <- x$B %*% Smat
      if(absorb=="A") {
        for(k in 1:length(x$A)) x$A[[k]] <- x$A[[k]] %*% Smat
        x$Phi <- Smat %*% x$Phi %*% Smat
      } else if(absorb=="C"){
        x$C <- x$C %*% Smat
      } else {
        x$D <- x$D %*% Smat
      }
      return(x)
      
    } else if(mode=="C"){
      
      Csign <- sign(colMeans(x$C^3))
      if(any(Csign == 0)) Csign[Csign == 0] <- 1
      svec <- newsign*Csign
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$C <- x$C %*% Smat
      if(absorb=="A") {
        for(k in 1:length(x$A)) x$A[[k]] <- x$A[[k]] %*% Smat
        x$Phi <- Smat %*% x$Phi %*% Smat
      } else if(absorb=="B"){
        x$B <- x$B %*% Smat
      } else {
        x$D <- x$D %*% Smat
      }
      return(x)
      
    } else if(mode=="D"){
      
      Dsign <- sign(colMeans(x$D^3))
      if(any(Dsign == 0)) Dsign[Dsign == 0] <- 1
      svec <- newsign*Dsign
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$D <- x$D %*% Smat
      if(absorb=="A") {
        for(k in 1:length(x$A)) x$A[[k]] <- x$A[[k]] %*% Smat
        x$Phi <- Smat %*% x$Phi %*% Smat
      } else if(absorb=="B"){
        x$B <- x$B %*% Smat
      } else {
        x$C <- x$C %*% Smat
      }
      return(x)
      
    }
    
  }