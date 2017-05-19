rfunc <- 
  function(x, df=7L, degree=3L, type=c("smo","mon","uni","per"), nonneg=FALSE){
    # Generate Random Function
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 12, 2017
    
    ### INITIAL CHECKS ###
    
    # check 'x' and define 'n'
    x <- as.numeric(x)
    n <- length(x)
    
    # check 'df'
    df <- as.integer(df[1])
    if(df < 4L) stop("Input 'df' must be an integer >= 4.")
    
    # check 'degree'
    degree <- as.integer(degree[1])
    if(degree < 1L) stop("Input 'degree' must be an integer >= 1.")
    
    # check 'type'
    type <- as.character(type[1])
    if(!any(type == c("smo","mon","uni","per"))) stop("Input 'type' needs to be one of the four options: \n'smo', 'mon', 'uni', or 'per'.")
    
    # check 'nonneg'
    nonneg <- as.logical(nonneg[1])
    
    
    ### GENERATE FUNCTION ###
    
    if(type == "smo"){
      X <- MsplineBasis(x, df=df, degree=degree, intercept=TRUE)$X
      if(nonneg){
        return(as.numeric(X %*% runif(df)))
      } else {
        return(as.numeric(X %*% rnorm(df)))
      }
    } else if (type == "per"){
      X <- cbind(1, MsplineBasis(x, df=df, degree=degree, intercept=FALSE)$X[,-df])
      if(nonneg){
        return(as.numeric(X %*% runif(df)))
      } else {
        return(as.numeric(X %*% rnorm(df)))
      }
    } else if (type =="mon"){
      X <- cbind(1, IsplineBasis(x, df=df-1L, degree=degree)$X)
      if(nonneg){
        return(as.numeric(X %*% runif(df)))
      } else {
        return(as.numeric(X %*% c(rnorm(1), runif(df-1L))))
      }
    } else if (type =="uni"){
      X <- cbind(1, IsplineBasis(x, df=df-1L, degree=degree)$X)
      modeIndex <- sample.int(df-1L, size=1)
      signvec <- c(1,rep(c(1,-1), c(modeIndex, df-modeIndex-1)))
      fhat <- X %*% (c(rnorm(1), runif(df-1L)) * signvec)
      if(nonneg){
        minfhat <- min(fhat)
        if(minfhat < 0) fhat <- fhat - minfhat + runif(1)
      }
      return(as.numeric(fhat))
    }
    
  }