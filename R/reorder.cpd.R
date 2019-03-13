reorder.cpd <-
  function(x, neworder, ...){
    # Reorder Factors of fit CPD model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: March 12, 2019
    
    # check neworder
    neworder <- as.integer(neworder)
    nfac <- ncol(x$A[[1]])
    if(length(neworder) != nfac) stop("Incorrect input for 'neworder'. Must have length equal to number of factors.")
    if(!identical(seq(1L,nfac), sort(neworder))) stop(paste("Incorrect input for 'neworder'. Must be unique integers in range of 1 to",nfac))
    
    # reorder factors
    for(n in 1:length(x$A)){
      x$A[[n]] <- x$A[[n]][,neworder]
    }
    return(x)
    
  }