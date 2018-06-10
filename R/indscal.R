indscal <- 
  function(X, nfac, nstart = 10, const = NULL, control = NULL,
           type = c("dissimilarity", "similarity"),
           Bfixed = NULL, Bstart = NULL, Bstruc = NULL, Bmodes = NULL,
           Cfixed = NULL, Cstart = NULL, Cstruc = NULL, Cmodes = NULL,
           maxit = 500, ctol = 1e-4, parallel = FALSE, cl = NULL,
           output = c("best", "all"), verbose = TRUE, backfit = FALSE){
    # Individual Differences Scaling (INDSCAL)
    # via alternating least squares (ALS) with optional constraints
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 25, 2018
    
    # check 'X' input
    if(is.list(X)){
      X <- lapply(X,as.matrix)
      px <- ncol(X[[1]])
      nx <- length(X)
      if(sum((sapply(X,dim)-matrix(px,2,nx))^2)>0L){
        stop("Input 'X' must be list of square matrices with same dimension.")
      }
      xdim <- c(px,px,nx)
      X <- array(unlist(X),dim=xdim)
    } else{
      xdim <- dim(X)
      if(length(xdim)!=3L){stop("Input 'X' must be 3-way array")}
      if(xdim[1]!=xdim[2]){stop("Input 'X' must be 3-way array with dim(X)[1]==dim(X)[2].")}
    }
    if(any(is.na(X)) | any(is.nan(X)) | any(is.infinite(X))){stop("Input 'X' cannot contain NA, NaN, or Inf values")}
    if(type[1]=="dissimilarity") { X <- array(apply(X,3,ed2sp),dim=xdim) }
    xcx <- sumsq(X)
    
    # check 'nfac' and 'nstart' inputs
    nfac <- as.integer(nfac[1])
    if(nfac<1L){stop("Input 'nfac' must be positive integer")}
    nstart <- as.integer(nstart[1])
    if(nstart<1L){stop("Input 'nstart' must be positive integer")}
    
    # check 'maxit' and 'ctol' inputs
    maxit <- as.integer(maxit[1])
    if(maxit<1L){stop("Input 'maxit' must be positive integer")}
    ctol <- as.numeric(ctol[1])
    if(ctol<0L){stop("Input 'ctol' must be positive numeric")}
    
    # const types: old, new, spline, and non-negative
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
    const.nonneg <- c("nonneg", "pernon", "smonon", "smpeno", 
                      "ortnon", "monnon", "mosmno", "uninon", 
                      "unpeno", "unsmno", "unsmpn")
    const.unimod <- c("unimod", "uninon", "uniper", "unpeno", 
                      "unismo", "unsmno", "unsmpe", "unsmpn")
    
    # check 'const' input
    if(is.null(const)){
      const <- c("uncons", "nonneg")
    } else {
      if(length(const) != 2L) stop("Input 'const' must be a two element vector specifying constraint for each mode")
      if(is.integer(const) | is.numeric(const)){
        const <- as.integer(const)
        if(any(is.na(const))) const[is.na(const)] <- 0L
        if(any(const < 0L) | any(const > 6L)) stop("Input 'const' must contain characters (preferred) or integers between 0 and 6.")
        const <- const.oldtypes[const + 1L]
      } else {
        const <- as.character(const)
        if(any(is.na(const))) const[is.na(const)] <- "uncons"
        cid <- pmatch(const, const.newtypes, duplicates.ok = TRUE)
        if(any(is.na(cid))) stop("Input 'const' must be a character vector of length 2\n with each element matching one of the 24 available options.")
        const <- const.newtypes[cid]
      }
      if(!any(const[2] == const.nonneg)){
        warning("Mode C weights must be one of the 11 possible non-negative options.\n  Resetting const[2] to 'nonneg' for current fitting.")
        const[2] <- "nonneg"
      } 
    } # end if(is.null(const))
    
    # check 'control' input
    const3 <- c(rep(const[1],2), const[2])
    if(!is.null(control)){
      if(!is.null(control$df)) control$df <- c(rep(control$df[1],2), control$df[2])
      if(!is.null(control$degree)) control$degree <- c(rep(control$degree[1],2), control$degree[2])
      if(!is.null(control$intercept)) control$intercept <- c(rep(control$intercept[1],2), control$intercept[2])
    }
    convars <- list(const=const3, df=control$df, degree=control$degree, intercept=control$intercept)
    control <- do.call("const.control", convars)
    for(j in 1:length(const)){
      if(any(const[j] == const.smooth)){
        if(control$df[j] >= xdim[j]) {
          warning(paste0("Input control$df[",j,"] = ",control$df[j]," >= dim(X)[",j,"] = ",xdim[j],"\nResetting control$df[",j,"] = ",xdim[j]-1,"."))
          control$df[j] <- xdim[j] - 1
        }
      }
    }
    
    # check 'Bfixed' and 'Cfixed' inputs
    if(!is.null(Bfixed)){
      Bfixed <- as.matrix(Bfixed)
      if(nrow(Bfixed)!=xdim[2]){stop("Input 'Bfixed' must have the same number of rows as dim(X)[2]")}
      if(ncol(Bfixed)!=nfac){stop("Input 'Bfixed' must have 'nfac' columns")}
    }
    if(!is.null(Cfixed)){
      Cfixed <- as.matrix(Cfixed)
      if(nrow(Cfixed)!=xdim[3]){stop("Input 'Cfixed' must have the same number of rows as dim(X)[3]")}
      if(ncol(Cfixed)!=nfac){stop("Input 'Cfixed' must have 'nfac' columns")}
      if(any(Cfixed < 0)) stop("Input 'Cfixed' contains negative values, which are not allowed.")
    }
    
    # check 'Bstart' and 'Cstart' inputs
    if(!is.null(Bstart)){
      Bstart <- as.matrix(Bstart)
      if(nrow(Bstart)!=xdim[2]){stop("Input 'Bstart' must have the same number of rows as dim(X)[2]")}
      if(ncol(Bstart)!=nfac){stop("Input 'Bstart' must have 'nfac' columns")}
    }
    if(!is.null(Cstart)){
      Cstart <- as.matrix(Cstart)
      if(nrow(Cstart)!=xdim[3]){stop("Input 'Cstart' must have the same number of rows as dim(X)[3]")}
      if(ncol(Cstart)!=nfac){stop("Input 'Cstart' must have 'nfac' columns")}
      if(any(Cstart < 0)) stop("Input 'Cstart' contains negative values, which are not allowed.")
    }
    
    # check 'Bstruc' and 'Cstruc' inputs
    if(!is.null(Bstruc)){
      Bstruc <- as.matrix(Bstruc)
      if(nrow(Bstruc)!=xdim[2]){stop("Input 'Bstruc' must have the same number of rows as dim(X)[2]")}
      if(ncol(Bstruc)!=nfac){stop("Input 'Bstruc' must have 'nfac' columns")}
      if(any("logical"!=c(apply(Bstruc,1:2,class)))){stop("Input 'Bstruc' must be a matrix with logical (TRUE/FALSE) entries.")}
    }
    if(!is.null(Cstruc)){
      Cstruc <- as.matrix(Cstruc)
      if(nrow(Cstruc)!=xdim[3]){stop("Input 'Cstruc' must have the same number of rows as dim(X)[3]")}
      if(ncol(Cstruc)!=nfac){stop("Input 'Cstruc' must have 'nfac' columns")}
      if(any("logical"!=c(apply(Cstruc,1:2,class)))){stop("Input 'Cstruc' must be a matrix with logical (TRUE/FALSE) entries.")}
    }
    
    # check 'Bmodes' and 'Cmodes' inputs
    if(any(const[2] == const.unimod) && !is.null(Bmodes)){
      Bmodes <- as.matrix(Bmodes)
      if(nrow(Bmodes) != 2L | ncol(Bmodes) != nfac) stop("Input 'Bmodes' must be a 2 x nfac matrix giving the\n minimum (row 1) and maximum (row 2) allowable mode for each factor.")
      Bmodes <- matrix(as.integer(Bmodes), nrow = 2L, ncol = nfac)
      if(any(Bmodes < 1L)) stop("First row of 'Bmodes' must contain integers greater than or equal to one.")
      if(any(Bmodes > xdim[2])) stop("Second row of 'Bmodes' must contain integers less than or equal to dim(X)[2].")
      if(any((Bmodes[2,] - Bmodes[1,]) < 0)) stop("Input 'Bmodes' must satisfy:  Bmodes[1,r] <= Bmodes[2,r]")
    }
    if(any(const[3] == const.unimod) && !is.null(Cmodes)){
      Cmodes <- as.matrix(Cmodes)
      if(nrow(Cmodes) != 2L | ncol(Cmodes) != nfac) stop("Input 'Cmodes' must be a 2 x nfac matrix giving the\n minimum (row 1) and maximum (row 2) allowable mode for each factor.")
      Cmodes <- matrix(as.integer(Cmodes), nrow = 2L, ncol = nfac)
      if(any(Cmodes < 1L)) stop("First row of 'Cmodes' must contain integers greater than or equal to one.")
      if(any(Cmodes > xdim[3])) stop("Second row of 'Cmodes' must contain integers less than or equal to dim(X)[3].")
      if(any((Cmodes[2,] - Cmodes[1,]) < 0)) stop("Input 'Cmodes' must satisfy:  Cmodes[1,r] <= Cmodes[2,r]")
    }
    
    # check 'parallel' and 'cl' inputs
    if(parallel && !any(class(cl)=="cluster")){
      stop("Input 'cl' must be cluster (created by makeCluster) when parallel=TRUE \n  See examples in documentation:  ?parafac")
    }
    
    # indscal fitting
    if(parallel){
      nstartlist <- vector("list",nstart)
      nstartlist[1:nstart] <- nfac
      pfaclist <- parLapply(cl = cl, X = nstartlist, fun = "parafac_3way", data = X, xcx = xcx,
                            const = const3, control = control, maxit = maxit, ctol = ctol, 
                            Afixed = Bfixed, Bfixed = Bfixed, Cfixed = Cfixed, 
                            Astart = Bstart, Bstart = Bstart, Cstart = Cstart,
                            Astruc = Bstruc, Bstruc = Bstruc, Cstruc = Cstruc,
                            Amodes = Bmodes, Bmodes = Bmodes, Cmodes = Cmodes,
                            backfit = backfit)
    } else {
      
      if(verbose) pbar <- txtProgressBar(min = 0, max = nstart, style = 3)
      
      if(output[1] == "best"){
        
        pfac <- parafac_3way(data = X, nfac = nfac, xcx = xcx, const = const3, 
                             control = control, maxit = maxit, ctol = ctol, 
                             Afixed = Bfixed, Bfixed = Bfixed, Cfixed = Cfixed, 
                             Astart = Bstart, Bstart = Bstart, Cstart = Cstart,
                             Astruc = Bstruc, Bstruc = Bstruc, Cstruc = Cstruc,
                             Amodes = Bmodes, Bmodes = Bmodes, Cmodes = Cmodes,
                             backfit = backfit)
        if(verbose) setTxtProgressBar(pbar, 1)
        if(nstart > 1L){
          for(j in 2:nstart){
            pnew <- parafac_3way(data = X, nfac = nfac, xcx = xcx, const = const3,
                                 control = control, maxit = maxit, ctol = ctol, 
                                 Afixed = Bfixed, Bfixed = Bfixed, Cfixed = Cfixed, 
                                 Astart = Bstart, Bstart = Bstart, Cstart = Cstart,
                                 Astruc = Bstruc, Bstruc = Bstruc, Cstruc = Cstruc,
                                 Amodes = Bmodes, Bmodes = Bmodes, Cmodes = Cmodes,
                                 backfit = backfit)
            if(pnew$SSE < pfac$SSE) pfac <- pnew
            if(verbose) setTxtProgressBar(pbar, j)
          } # end for(j in 2:nstart)
        } # end if(nstart > 1L)
        
      } else {
        
        pfaclist <- vector("list", nstart)
        for(j in 1:nstart){
          pfaclist[[j]] <- parafac_3way(data = X, nfac = nfac, xcx = xcx, const = const3,
                                        control = control, maxit = maxit, ctol = ctol, 
                                        Afixed = Bfixed, Bfixed = Bfixed, Cfixed = Cfixed, 
                                        Astart = Bstart, Bstart = Bstart, Cstart = Cstart,
                                        Astruc = Bstruc, Bstruc = Bstruc, Cstruc = Cstruc,
                                        Amodes = Bmodes, Bmodes = Bmodes, Cmodes = Cmodes,
                                        backfit = backfit)
          if(verbose) setTxtProgressBar(pbar, j)
        }
        
      } # end if(output[1] == "best")
      
      if(verbose) close(pbar)
      
    } # end if(parallel)
    
    # output results
    if(output[1] == "best"){
      if(parallel){
        SSE <- sapply(pfaclist, function(x) x$SSE)
        pfac <- pfaclist[[which.min(SSE)]]
      }
      pfac <- pfac[-1]
      pfac$SSE <- sum((X - array(tcrossprod(pfac$B, krprod(pfac$C, pfac$B)), dim = xdim))^2)
      pfac$GCV <- (pfac$SSE / prod(xdim)) / (1 - sum(pfac$edf)/prod(xdim))^2
      pfac$const <- const
      class(pfac) <- "indscal"
      return(pfac)
    } else {
      ifun <- function(x){
        x <- x[-1]
        x$SSE <- sum((X - array(tcrossprod(x$B, krprod(x$C, x$B)), dim = xdim))^2)
        x$GCV <- (x$SSE / prod(xdim)) / (1 - sum(x$edf)/prod(xdim))^2
        x$const <- const
        class(x) <- "indscal"
        x
      }
      pfaclist <- lapply(pfaclist, ifun)
      return(pfaclist)
    } # end if(output[1] == "best")
    
  } # end indscal.R