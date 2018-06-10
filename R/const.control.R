const.control <-
  function(const, df = NULL, degree = NULL, intercept = NULL){
    
    # check nway
    nway <- length(const)
    if(!any(nway == c(3L, 4L))) stop("Input 'const' must be a vector of length 3 or 4.")
    names(const) <- LETTERS[1:nway]
    
    # const types (old and new)
    const.oldtypes <- c("uncons", "orthog", "nonneg", 
                        "unismo", "monsmo", "smoper", "smooth")
    const.newtypes <- c("uncons", "nonneg", "period", "pernon",
                        "smooth", "smonon", "smoper", "smpeno",
                        "orthog", "ortnon", "ortsmo", "orsmpe",
                        "moninc", "monnon", "monsmo", "mosmno",
                        "unimod", "uninon", "uniper", "unpeno", 
                        "unismo", "unsmno", "unsmpe", "unsmpn")
    const.smooth <- c("smooth", "smonon", "smoper", "smpeno",
                      "ortsmo", "orsmpe", "monsmo", "mosmno", 
                      "unismo", "unsmno", "unsmpe", "unsmpn")
    
    # check const
    const <- c(const)
    if(is.integer(const) | is.numeric(const)){
      const <- as.integer(const)
      if(any(is.na(const))) const[is.na(const)] <- 0L
      if(any(const < 0L) | any(const > 6L)) stop("Input 'const' must contain characters or integers between 0 and 6.")
      const <- const.oldtypes[const + 1L]
    } else {
      const <- as.character(const)
      if(any(is.na(const))) const[is.na(const)] <- "uncons"
      cid <- pmatch(const, const.newtypes, duplicates.ok = TRUE)
      if(any(is.na(cid))) stop("Input 'const' must be a character vector of length 3\n with each element matching one of the 18 available options.")
      const <- const.newtypes[cid]
    }
    
    # check if any constraints with options
    ix <- which(!is.na(match(const, const.smooth)))
    if(length(ix) > 0L){
      
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
      
      # check intercept
      if(is.null(intercept)) intercept <- rep(TRUE, nway)
      intercept <- as.logical(intercept)
      if(length(intercept) != nway) intercept <- rep(intercept[1], nway)
      intercept[-ix] <- NA
      names(intercept) <- LETTERS[1:nway]
      
      # check df
      if(is.null(df)) df <- rep(NA, nway)
      df <- as.integer(df)
      if(length(df) != nway) df <- rep(df[1], nway)
      for(k in ix){
        mindf <- degree[k] + 1 - !intercept[k]
        if(is.na(df[k]) | is.nan(df[k]) | is.null(df[k]) | is.infinite(df[k])){
          df[k] <- max(7, mindf)
        } else {
          if(df[k] < mindf) {
            warning("Input 'df' is too small: need df[k] >= degree[k] + 1 - !intercept[k].\n Resetting to minimum possible df.")
            df[k] <- mindf
          }
        }
      }
      df[-ix] <- NA
      names(df) <- LETTERS[1:nway]
      
    } else {
      
      df <- degree <- intercept <- rep(NA, nway)
      
    } # end if(length(ix) > 2L)
    
    # return results
    constlist <- list(const = const, df = df, 
                      degree = degree, intercept = intercept)
    return(constlist)
    
  }