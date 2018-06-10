sumsq <- 
  function(X, na.rm = FALSE){
    # sum of squares of input object
    
    if(is.list(X)){
      xss <- sum(sapply(X, sumsq, na.rm = na.rm))
    } else {
      xss <- sum(X^2, na.rm = na.rm)
    }
    
    xss
    
  }