const.control <-
  function(const, df = NULL, degree = NULL, nonneg = NULL){
    
    # check nway and const
    nway <- length(const)
    if(!any(nway == c(3L, 4L))) stop("Input 'const' must be a vector of length 3 or 4.")
    const <- as.integer(const)
    if(any(const < 0L) | any(const > 6L)) stop("Input 'const' must contain integers between 0 and 6.")
    names(const) <- LETTERS[1:nway]
    ix <- which(const > 2L)
    
    # check if any constraints with options
    if(length(ix) > 0L){
      
      # check df
      if(is.null(df)) df <- rep(NA, nway)
      df <- as.integer(df)
      if(length(df) != nway) df <- rep(df[1], nway)
      for(k in ix){
        if(is.na(df[k]) | is.nan(df[k]) | is.null(df[k]) | is.infinite(df[k])){
          df[k] <- 7
        } else {
          if(df[k] < 4L) stop(paste0("Input 'df' must be of a vector of length ",nway," with elements >= 4."))
        }
      }
      df[-ix] <- NA
      names(df) <- LETTERS[1:nway]
      
      # check degree
      if(is.null(degree)) degree <- rep(NA, nway)
      degree <- as.integer(degree)
      if(length(degree) != nway) degree <- rep(degree[1], nway)
      for(k in ix){
        if(is.na(degree[k]) | is.nan(degree[k]) | is.null(degree[k]) | is.infinite(degree[k])){
          degree[k] <- 3
        } else {
          if(degree[k] < 1L) stop(paste0("Input 'degree' must be of a vector of length ",nway," with elements >= 1."))
        }
      }
      degree[-ix] <- NA
      names(degree) <- LETTERS[1:nway]
      
      # check nonneg
      if(is.null(nonneg)) nonneg <- rep(FALSE, nway)
      nonneg <- as.logical(nonneg)
      if(length(nonneg) != nway) nonneg <- rep(nonneg[1], nway)
      nonneg[-ix] <- NA
      names(nonneg) <- LETTERS[1:nway]
      
    } else {
      
      df <- degree <- nonneg <- rep(NA, nway)
      
    } # end if(length(ix) > 2L)
    
    return(list(const=const, df=df, degree=degree, nonneg=nonneg))
    
  }