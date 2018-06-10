meansq <- 
  function(X, na.rm = FALSE){
    
    if(is.list(X)){
      xms <- mean(unlist(X)^2, na.rm = na.rm)
    } else {
      xms <- mean(X^2, na.rm = na.rm)
    }
    
    xms
    
  }