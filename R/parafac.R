parafac <- 
  function(X, nfac, nstart = 10, const = NULL, control = NULL,
           Afixed = NULL, Bfixed = NULL, Cfixed = NULL, Dfixed = NULL,
           Astart = NULL, Bstart = NULL, Cstart = NULL, Dstart = NULL,
           Astruc = NULL, Bstruc = NULL, Cstruc = NULL, Dstruc = NULL,
           Amodes = NULL, Bmodes = NULL, Cmodes = NULL, Dmodes = NULL,
           maxit = 500, ctol = 1e-4, parallel = FALSE, cl = NULL, 
           output = c("best", "all"), verbose = TRUE, backfit = FALSE){
    # 3-way or 4-way Parallel Factor Analysis (Parafac)
    # via alternating least squares (ALS) with optional constraints
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 25, 2018
    
    # check 'X' input
    xdim <- dim(X)
    lxdim <- length(xdim)
    if((lxdim < 3L) | (lxdim > 4L)) stop("Input 'X' must be 3-way or 4-way array") 
    if(any(is.nan(X)) | any(is.infinite(X))) stop("Input 'X' cannot contain NaN or Inf values") 
    if(any(is.na(X))){
      missingdata <- TRUE
      naid <- which(is.na(X))
    } else {
      missingdata <- FALSE
      xcx <- sumsq(X)
    }
    
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
    
    # const types: old, new, and spline
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
    const.unimod <- c("unimod", "uninon", "uniper", "unpeno", 
                      "unismo", "unsmno", "unsmpe", "unsmpn")
    
    
    # check 'const' input
    if(is.null(const)){
      const <- rep("uncons", lxdim)
    } else {
      if(length(const) != lxdim) stop("Input 'const' must be ",lxdim," element vector specifying constraint for each mode")
      if(is.integer(const) | is.numeric(const)){
        const <- as.integer(const)
        if(any(is.na(const))) const[is.na(const)] <- 0L
        if(any(const < 0L) | any(const > 6L)) stop("Input 'const' must contain characters (preferred) or integers between 0 and 6.")
        const <- const.oldtypes[const + 1L]
      } else {
        const <- as.character(const)
        if(any(is.na(const))) const[is.na(const)] <- "uncons"
        cid <- pmatch(const, const.newtypes, duplicates.ok = TRUE)
        if(any(is.na(cid))) stop("Input 'const' must be a character vector of length ",lxdim,"\n with each element matching one of the 24 available options.")
        const <- const.newtypes[cid]
      }
    } # end if(is.null(const))
    
    # check 'control' input
    convars <- list(const=const, df=control$df, degree=control$degree, intercept=control$intercept)
    control <- do.call("const.control", convars)
    for(j in 1:length(const)){
      if(any(const[j] == const.smooth)){
        if(control$df[j] >= xdim[j]) {
          warning(paste0("Input control$df[",j,"] = ",control$df[j]," >= dim(X)[",j,"] = ",xdim[j],"\nResetting control$df[",j,"] = ",xdim[j]-1,"."))
          control$df[j] <- xdim[j] - 1
        }
      }
    }
    
    # check 'Afixed' and 'Bfixed' and 'Cfixed' inputs
    if(!is.null(Afixed)){
      Afixed <- as.matrix(Afixed)
      if(nrow(Afixed)!=xdim[1]){stop("Input 'Afixed' must have the same number of rows as dim(X)[1]")}
      if(ncol(Afixed)!=nfac){stop("Input 'Afixed' must have 'nfac' columns")}
    }
    if(!is.null(Bfixed)){
      Bfixed <- as.matrix(Bfixed)
      if(nrow(Bfixed)!=xdim[2]){stop("Input 'Bfixed' must have the same number of rows as dim(X)[2]")}
      if(ncol(Bfixed)!=nfac){stop("Input 'Bfixed' must have 'nfac' columns")}
    }
    if(!is.null(Cfixed)){
      Cfixed <- as.matrix(Cfixed)
      if(nrow(Cfixed)!=xdim[3]){stop("Input 'Cfixed' must have the same number of rows as dim(X)[3]")}
      if(ncol(Cfixed)!=nfac){stop("Input 'Cfixed' must have 'nfac' columns")}
    }
    
    # check 'Astart' and 'Bstart' and 'Cstart' inputs
    if(!is.null(Astart)){
      Astart <- as.matrix(Astart)
      if(nrow(Astart)!=xdim[1]){stop("Input 'Astart' must have the same number of rows as dim(X)[1]")}
      if(ncol(Astart)!=nfac){stop("Input 'Astart' must have 'nfac' columns")}
    }
    if(!is.null(Bstart)){
      Bstart <- as.matrix(Bstart)
      if(nrow(Bstart)!=xdim[2]){stop("Input 'Bstart' must have the same number of rows as dim(X)[2]")}
      if(ncol(Bstart)!=nfac){stop("Input 'Bstart' must have 'nfac' columns")}
    }
    if(!is.null(Cstart)){
      Cstart <- as.matrix(Cstart)
      if(nrow(Cstart)!=xdim[3]){stop("Input 'Cstart' must have the same number of rows as dim(X)[3]")}
      if(ncol(Cstart)!=nfac){stop("Input 'Cstart' must have 'nfac' columns")}
    }
    
    # check 'Astruc' and 'Bstruc' and 'Cstruc' inputs
    if(!is.null(Astruc)){
      Astruc <- as.matrix(Astruc)
      if(nrow(Astruc)!=xdim[1]){stop("Input 'Astruc' must have the same number of rows as dim(X)[1]")}
      if(ncol(Astruc)!=nfac){stop("Input 'Astruc' must have 'nfac' columns")}
      if(any("logical"!=c(apply(Astruc,1:2,class)))){stop("Input 'Astruc' must be a matrix with logical (TRUE/FALSE) entries.")}
    }
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
    
    # check 'Amodes' and 'Bmodes' and 'Cmodes' inputs
    if(any(const[1] == const.unimod) && !is.null(Amodes)){
      Amodes <- as.matrix(Amodes)
      if(nrow(Amodes) != 2L | ncol(Amodes) != nfac) stop("Input 'Amodes' must be a 2 x nfac matrix giving the\n minimum (row 1) and maximum (row 2) allowable mode for each factor.")
      Amodes <- matrix(as.integer(Amodes), nrow = 2L, ncol = nfac)
      if(any(Amodes < 1L)) stop("First row of 'Amodes' must contain integers greater than or equal to one.")
      if(any(Amodes > xdim[1])) stop("Second row of 'Amodes' must contain integers less than or equal to dim(X)[1].")
      if(any((Amodes[2,] - Amodes[1,]) < 0)) stop("Input 'Amodes' must satisfy:  Amodes[1,r] <= Amodes[2,r]")
    }
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
    
    # parafac fitting
    if(lxdim==3L){
      if(parallel){
        nstartlist <- vector("list",nstart)
        nstartlist[1:nstart] <- nfac
        if(missingdata){
          pfaclist <- parLapply(cl = cl, X = nstartlist, fun = "parafac_3wayna",
                                data = X, naid = naid, const = const, 
                                control = control, maxit = maxit, ctol = ctol,
                                Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed,
                                Astart = Astart, Bstart = Bstart, Cstart = Cstart,
                                Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc,
                                Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes,
                                backfit = backfit)
        } else {
          pfaclist <- parLapply(cl = cl, X = nstartlist, fun = "parafac_3way",
                                data = X, xcx = xcx, const = const,
                                control = control, maxit = maxit, ctol = ctol,
                                Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed,
                                Astart = Astart, Bstart = Bstart, Cstart = Cstart,
                                Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc,
                                Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes,
                                backfit = backfit)
        } # end if(missingdata)
      } else {
        
        if(verbose) pbar <- txtProgressBar(min = 0, max = nstart, style = 3)
        
        if(output[1] == "best"){
          if(missingdata){
            pfac <- parafac_3wayna(data = X, nfac = nfac, naid = naid, const = const,
                                   control = control, maxit = maxit, ctol = ctol,
                                   Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed,
                                   Astart = Astart, Bstart = Bstart, Cstart = Cstart,
                                   Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc,
                                   Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes,
                                   backfit = backfit)
            if(verbose) setTxtProgressBar(pbar, 1)
            if(nstart > 1L){
              for(j in 2:nstart){
                pnew <- parafac_3wayna(data = X, nfac = nfac, naid = naid, const = const,
                                       control = control, maxit = maxit, ctol = ctol,
                                       Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed,
                                       Astart = Astart, Bstart = Bstart, Cstart = Cstart,
                                       Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc,
                                       Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes,
                                       backfit = backfit)
                if(pnew$SSE < pfac$SSE) pfac <- pnew
                if(verbose) setTxtProgressBar(pbar, j)
              } # end for(j in 2:nstart)
            } # end if(nstart > 1L)
          } else {
            pfac <- parafac_3way(data = X, nfac = nfac, xcx = xcx, const = const,
                                 control = control, maxit = maxit, ctol = ctol,
                                 Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed,
                                 Astart = Astart, Bstart = Bstart, Cstart = Cstart,
                                 Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc,
                                 Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes,
                                 backfit = backfit)
            if(verbose) setTxtProgressBar(pbar, 1)
            if(nstart > 1L){
              for(j in 2:nstart){
                pnew <- parafac_3way(data = X, nfac = nfac, xcx = xcx, const = const,
                                     control = control, maxit = maxit, ctol = ctol,
                                     Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed,
                                     Astart = Astart, Bstart = Bstart, Cstart = Cstart,
                                     Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc,
                                     Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes,
                                     backfit = backfit)
                if(pnew$SSE < pfac$SSE) pfac <- pnew
                if(verbose) setTxtProgressBar(pbar, j)
              } # end for(j in 2:nstart)
            } # end if(nstart > 1L)
          } # end if(missingdata)
          
        } else {
          
          pfaclist <- vector("list",nstart)
          if(missingdata){
            for(j in 1:nstart){
              pfaclist[[j]] <- parafac_3wayna(data = X, nfac = nfac, naid = naid, const = const,
                                              control = control, maxit = maxit, ctol = ctol,
                                              Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed,
                                              Astart = Astart, Bstart = Bstart, Cstart = Cstart,
                                              Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc,
                                              Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes,
                                              backfit = backfit)
              if(verbose) setTxtProgressBar(pbar, j)
            }
          } else {
            for(j in 1:nstart){
              pfaclist[[j]] <- parafac_3way(data = X, nfac = nfac, xcx = xcx, const = const,
                                            control = control, maxit = maxit, ctol = ctol,
                                            Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed,
                                            Astart = Astart, Bstart = Bstart, Cstart = Cstart,
                                            Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc,
                                            Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes,
                                            backfit = backfit)
              if(verbose) setTxtProgressBar(pbar, j)
            }
          } # end if(missingdata)
        } # end if(output[1] == "best")
        
        if(verbose) close(pbar)
        
      } # end if(parallel)
      
    } else if(lxdim==4L){
      # check 'Dfixed', 'Dstart', 'Dstruc', and 'Dmodes' inputs
      if(!is.null(Dfixed)){
        Dfixed <- as.matrix(Dfixed)
        if(nrow(Dfixed)!=xdim[4]){stop("Input 'Dfixed' must have the same number of rows as dim(X)[4]")}
        if(ncol(Dfixed)!=nfac){stop("Input 'Dfixed' must have 'nfac' columns")}
      }
      if(!is.null(Dstart)){
        Dstart <- as.matrix(Dstart)
        if(nrow(Dstart)!=xdim[4]){stop("Input 'Dstart' must have the same number of rows as dim(X)[4]")}
        if(ncol(Dstart)!=nfac){stop("Input 'Dstart' must have 'nfac' columns")}
      }
      if(!is.null(Dstruc)){
        Dstruc <- as.matrix(Dstruc)
        if(nrow(Dstruc)!=xdim[4]){stop("Input 'Dstruc' must have the same number of rows as dim(X)[4]")}
        if(ncol(Dstruc)!=nfac){stop("Input 'Dstruc' must have 'nfac' columns")}
        if(any("logical"!=c(apply(Dstruc,1:2,class)))){stop("Input 'Dstruc' must be a matrix with logical (TRUE/FALSE) entries.")}
      }
      if(any(const[4] == const.unimod) && !is.null(Dmodes)){
        Dmodes <- as.matrix(Dmodes)
        if(nrow(Dmodes) != 2L | ncol(Dmodes) != nfac) stop("Input 'Dmodes' must be a 2 x nfac matrix giving the\n minimum (row 1) and maximum (row 2) allowable mode for each factor.")
        Dmodes <- matrix(as.integer(Dmodes), nrow = 2L, ncol = nfac)
        if(any(Dmodes < 1L)) stop("First row of 'Dmodes' must contain integers greater than or equal to one.")
        if(any(Dmodes > xdim[4])) stop("Second row of 'Dmodes' must contain integers less than or equal to dim(X)[4].")
        if(any((Dmodes[2,] - Dmodes[1,]) < 0)) stop("Input 'Dmodes' must satisfy:  Dmodes[1,r] <= Dmodes[2,r]")
      }
      
      # fit models
      if(parallel){
        nstartlist <- vector("list",nstart)
        nstartlist[1:nstart] <- nfac
        if(missingdata){
          pfaclist <- parLapply(cl = cl, X = nstartlist, fun = "parafac_4wayna", data = X, naid = naid,
                                const = const, control = control, maxit = maxit, ctol = ctol,
                                Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                                Astart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                                Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc, Dstruc = Dstruc,
                                Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes, Dmodes = Dmodes,
                                backfit = backfit)
        } else {
          pfaclist <- parLapply(cl = cl, X = nstartlist, fun = "parafac_4way", data = X, xcx = xcx,
                                const = const, control = control, maxit = maxit, ctol = ctol,
                                Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                                Astart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                                Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc, Dstruc = Dstruc,
                                Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes, Dmodes = Dmodes,
                                backfit = backfit)
        } # end if(missingdata)
      } else {
        
        if(verbose) pbar <- txtProgressBar(min = 0, max = nstart, style = 3)
        
        if(output[1] == "best"){
          
          if(missingdata){
            pfac <- parafac_4wayna(data = X, nfac = nfac, naid = naid, const = const, 
                                   control = control, maxit = maxit, ctol = ctol,
                                   Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                                   Astart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                                   Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc, Dstruc = Dstruc,
                                   Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes, Dmodes = Dmodes,
                                   backfit = backfit)
            if(verbose) setTxtProgressBar(pbar, 1)
            if(nstart > 1L){
              for(j in 2:nstart){
                pnew <- parafac_4wayna(data = X, nfac = nfac, naid = naid, const = const, 
                                       control = control, maxit = maxit, ctol = ctol,
                                       Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                                       Astart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                                       Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc, Dstruc = Dstruc,
                                       Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes, Dmodes = Dmodes,
                                       backfit = backfit)
                if(pnew$SSE < pfac$SSE) pfac <- pnew
                if(verbose) setTxtProgressBar(pbar, j)
              } # end for(j in 2:nstart)
            } # end if(nstart > 1L)
          } else {
            pfac <- parafac_4way(data = X, nfac = nfac, xcx = xcx, const = const, 
                                 control = control, maxit = maxit,ctol = ctol,
                                 Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                                 Astart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                                 Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc, Dstruc = Dstruc,
                                 Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes, Dmodes = Dmodes,
                                 backfit = backfit)
            if(verbose) setTxtProgressBar(pbar, 1)
            if(nstart > 1L){
              for(j in 2:nstart){
                pnew <- parafac_4way(data = X, nfac = nfac, xcx = xcx, const = const, 
                                     control = control, maxit = maxit,ctol = ctol,
                                     Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                                     Astart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                                     Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc, Dstruc = Dstruc,
                                     Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes, Dmodes = Dmodes,
                                     backfit = backfit)
                if(pnew$SSE < pfac$SSE) pfac <- pnew
                if(verbose) setTxtProgressBar(pbar, j)
              } # end for(j in 2:nstart)
            } # end if(nstart > 1L)
          } # end if(missingdata)
          
        } else {
          
          pfaclist <- vector("list",nstart)
          if(missingdata){
            for(j in 1:nstart){
              pfaclist[[j]] <- parafac_4wayna(data = X, nfac = nfac, naid = naid, const = const, 
                                              control = control, maxit = maxit, ctol = ctol,
                                              Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                                              Astart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                                              Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc, Dstruc = Dstruc,
                                              Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes, Dmodes = Dmodes,
                                              backfit = backfit)
              if(verbose) setTxtProgressBar(pbar, j)
            }
          } else {
            for(j in 1:nstart){
              pfaclist[[j]] <- parafac_4way(data = X, nfac = nfac, xcx = xcx, const = const, 
                                            control = control, maxit = maxit,ctol = ctol,
                                            Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                                            Astart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                                            Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc, Dstruc = Dstruc,
                                            Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes, Dmodes = Dmodes,
                                            backfit = backfit)
              if(verbose) setTxtProgressBar(pbar, j)
            }
          } # end if(missingdata)
          
        } # end if(output[1] == "best")
        
        if(verbose) close(pbar)
        
      } # end if(parallel)
    } # end if(lxdim==3L)
    
    # output results
    if(output[1] == "best"){
      if(parallel){
        SSE <- sapply(pfaclist,function(x) x$SSE)
        pfac <- pfaclist[[which.min(SSE)]]
      }
      class(pfac) <- "parafac"
      return(pfac)
    } else {
      pfaclist <- lapply(pfaclist, function(pfac) {
        class(pfac) <- "parafac"
        pfac
        })
      return(pfaclist)
    }
    
  } # end parafac.R