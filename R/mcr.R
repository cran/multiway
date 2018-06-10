mcr <- 
  function(X, Y, nfac = 1, alpha = 0.5, nstart = 10, 
           model = c("parafac", "parafac2", "tucker"),
           const = NULL, control = NULL, weights = NULL,
           Afixed = NULL, Bfixed = NULL, Cfixed = NULL, Dfixed = NULL,
           Astart = NULL, Bstart = NULL, Cstart = NULL, Dstart = NULL,
           Astruc = NULL, Bstruc = NULL, Cstruc = NULL, Dstruc = NULL,
           Amodes = NULL, Bmodes = NULL, Cmodes = NULL, Dmodes = NULL,
           maxit = 500, ctol = 1e-4, parallel = FALSE, cl = NULL, 
           output = c("best", "all"), verbose = TRUE, backfit = FALSE){
    # Multi-way Covariates Regression
    #   Smilde AK, Kiers HAL, Multiway covariates regression models,
    #   Journal of Chemometrics, 1999, 13, 31-48.
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: May 25, 2018
    
    
    ### check 'X'
    X <- as.array(X)
    xdim <- dim(X)
    lxdim <- length(xdim)
    if(lxdim != 3L) stop("Input 'X' must be 3-way array")
    
    ### check 'Y'
    Y <- as.matrix(Y)
    ydim <- dim(Y)
    lydim <- length(ydim)
    if(lydim != 2) stop("Input 'Y' must be vector or matrix")
    if(ydim[1] != xdim[3]) stop("Inputs 'X' and 'Y' incompatible: need dim(X)[3] == nrow(Y).")
    
    ### check for missing data
    if(any(any(is.na(X)) | any(is.nan(X)) | any(is.infinite(X)))){stop("Input 'X' cannot contain NA, NaN, or Inf values")}
    if(any(any(is.na(Y)) | any(is.nan(Y)) | any(is.infinite(Y)))){stop("Input 'Y' cannot contain NA, NaN, or Inf values")}
    ssx <- sum(X^2)
    ssy <- sum(Y^2)
    
    ### check 'alpha'
    alpha <- as.numeric(alpha[1])
    if(alpha < 0 | alpha > 1) stop("Input alpha must satisfy: 0 <= alpha <= 1")
    
    ### check 'model'
    model <- as.character(model[1])
    if(!any(model == c("parafac", "parafac2", "tucker"))) stop("Input 'model' must be 'tucker', 'parafac', or 'parafac2'.")
    
    ### check 'nfac' input
    nfac <- as.integer(nfac)
    if(any(nfac < 1L)) stop("Input 'nfac' must be positive integer")
    if(model == "parafac" | model == "parafac2") {
      nfac <- nfac[1]
    } else if(model == "tucker"){
      if(length(nfac) != 3L) nfac <- rep(nfac[1], 3)
    }
    
    ### check 'nstart' input
    nstart <- as.integer(nstart[1])
    if(nstart < 1L) stop("Input 'nstart' must be positive integer")
    
    ### check 'const' input
    const.unimod <- c("unimod", "uninon", "uniper", "unpeno", 
                      "unismo", "unsmno", "unsmpe", "unsmpn")
    if(is.null(const)){
      const <- rep("uncons", 4L)
    } else {
      const <- as.character(const)
      if(length(const) != 4L) const <- rep(const[1], 4L)
      if(any(is.na(const))) const[is.na(const)] <- "uncons"
      const.types <- c("uncons", "nonneg", "period", "pernon",
                       "smooth", "smonon", "smoper", "smpeno",
                       "orthog", "ortnon", "ortsmo", "orsmpe",
                       "moninc", "monnon", "monsmo", "mosmno",
                       "unimod", "uninon", "uniper", "unpeno", 
                       "unismo", "unsmno", "unsmpe", "unsmpn")
      cid <- pmatch(const, const.types, duplicates.ok = TRUE)
      if(any(is.na(cid))) stop("Input 'const' must be a character vector of length 3\n with each element matching one of the 18 available options.")
      const <- const.types[cid]
      pf2con <- c("uncons","orthog","ortnon","ortsmo","orsmpe","smooth","smoper")
      if(model == "parafac2" && !any(const[1] == pf2con)){
        stop("If model = 'parafac2', mode A constraint 'const[1]' must be one of:\n 'uncons', 'orthog', 'ortnon', 'ortsmo', 'orsmpe', 'smooth', 'smoper'")
      }
    }
    
    ### check 'control' input
    convars <- list(const=const, df=control$df, degree=control$degree, intercept=control$intercept)
    control <- do.call("const.control", convars)
    for(j in 1:length(const)){
      if(!is.na(control$df[j]) && control$df[j] >= xdim[j]){
        warning(paste0("Input control$df[",j,"] = ",control$df[j]," >= dim(X)[",j,"] = ",xdim[j],"\nResetting control$df[",j,"] = ",xdim[j]-1,"."))
        control$df[j] <- xdim[j] - 1
      }
    }
    
    ### check 'Afixed', 'Bfixed', 'Cfixed', and 'Dfixed'
    if(!is.null(Afixed)){
      Afixed <- as.matrix(Afixed)
      if(model == "parafac2"){
        if(nrow(Afixed) != nfac) stop("Need nrow(Afixed) == nfac when model == 'parafac2'.")
      } else {
        if(nrow(Afixed) != xdim[1]) stop("Input 'Afixed' is incompatible with input 'X'.\nNeed nrow(Afixed) == dim(X)[1]")
      }
      if(ncol(Afixed) != nfac) stop("Input 'Afixed' must have 'nfac' columns.")
    }
    if(!is.null(Bfixed)){
      Bfixed <- as.matrix(Bfixed)
      if(nrow(Bfixed) != xdim[2]) stop("Input 'Bfixed' is incompatible with input 'X'.\nNeed nrow(Bfixed) == dim(X)[2]")
      if(ncol(Bfixed) != nfac) stop("Input 'Bfixed' must have 'nfac' columns.")
    }
    if(!is.null(Cfixed)){
      Cfixed <- as.matrix(Cfixed)
      if(nrow(Cfixed) != xdim[3]) stop("Input 'Cfixed' is incompatible with input 'X'.\nNeed nrow(Cfixed) == dim(X)[3]")
      if(ncol(Cfixed) != nfac) stop("Input 'Cfixed' must have 'nfac' columns.")
    }
    if(!is.null(Dfixed)){
      Dfixed <- as.matrix(Dfixed)
      if(nrow(Dfixed) != ydim[2]) stop("Input 'Dfixed' is incompatible with input 'Y'.\nNeed nrow(Dfixed) == ncol(Y)")
      if(ncol(Dfixed) != nfac) stop("Input 'Dfixed' must have 'nfac' columns.")
    }
    
    ### check 'Astart', 'Bstart', 'Cstart', and 'Dstart'
    if(!is.null(Astart)){
      Astart <- as.matrix(Astart)
      if(model == "parafac2"){
        if(nrow(Astart) != nfac) stop("Need nrow(Astart) == nfac when model === 'parafac2'.")
      } else {
        if(nrow(Astart) != xdim[1]) stop("Input 'Astart' is incompatible with input 'X'.\nNeed nrow(Astart) == dim(X)[1]")
      }
      if(ncol(Astart) != nfac) stop("Input 'Astart' must have 'nfac' columns.")
    }
    if(!is.null(Bstart)){
      Bstart <- as.matrix(Bstart)
      if(nrow(Bstart) != xdim[2]) stop("Input 'Bstart' is incompatible with input 'X'.\nNeed nrow(Bstart) == dim(X)[2]")
      if(ncol(Bstart) != nfac) stop("Input 'Bstart' must have 'nfac' columns.")
    }
    if(!is.null(Cstart)){
      Cstart <- as.matrix(Cstart)
      if(nrow(Cstart) != xdim[3]) stop("Input 'Cstart' is incompatible with input 'X'.\nNeed nrow(Cstart) == dim(X)[3]")
      if(ncol(Cstart) != nfac) stop("Input 'Cstart' must have 'nfac' columns.")
    }
    if(!is.null(Dstart)){
      Dstart <- as.matrix(Dstart)
      if(nrow(Dstart) != ydim[2]) stop("Input 'Dstart' is incompatible with input 'Y'.\nNeed nrow(Dstart) == ncol(Y)")
      if(ncol(Dstart) != nfac) stop("Input 'Dstart' must have 'nfac' columns.")
    }
    
    ### check 'Astruc', 'Bstruc', and 'Dstruc'
    if(!is.null(Astruc)){
      Astruc <- as.matrix(Astruc)
      if(model == "parafac2"){
        if(nrow(Astruc) != nfac) stop("Need nrow(Astruc) == nfac when model === 'parafac2'.")
      } else {
        if(nrow(Astruc) != xdim[1]) stop("Input 'Astruc' is incompatible with input 'X'.\nNeed nrow(Astruc) == dim(X)[1]")
      }
      if(ncol(Astruc) != nfac) stop("Input 'Astruc' must have 'nfac' columns.")
      if(any("logical" != c(apply(Astruc, 1:2, class)))) stop("Input 'Astruc' must be a matrix with logical (TRUE/FALSE) entries.")
    }
    if(!is.null(Bstruc)){
      Bstruc <- as.matrix(Bstruc)
      if(nrow(Bstruc) != xdim[2]) stop("Input 'Bstruc' is incompatible with input 'X'.\nNeed nrow(Bstruc) == dim(X)[2]")
      if(ncol(Bstruc) != nfac) stop("Input 'Bstruc' must have 'nfac' columns.")
      if(any("logical" != c(apply(Bstruc, 1:2, class)))) stop("Input 'Bstruc' must be a matrix with logical (TRUE/FALSE) entries.")
    }
    Cstruc <- NULL
    # if(!is.null(Cstruc)){
    #   Cstruc <- as.matrix(Cstruc)
    #   if(nrow(Cstruc) != xdim[3]) stop("Input 'Cstruc' is incompatible with input 'X'.\nNeed nrow(Cstruc) == dim(X)[3]")
    #   if(ncol(Cstruc) != nfac) stop("Input 'Cstruc' must have 'nfac' columns.")
    #   if(any("logical" != c(apply(Cstruc, 1:2, class)))) stop("Input 'Cstruc' must be a matrix with logical (TRUE/FALSE) entries.")
    # }
    if(!is.null(Dstruc)){
      Dstruc <- as.matrix(Dstruc)
      if(nrow(Dstruc) != ydim[2]) stop("Input 'Dstruc' is incompatible with input 'Y'.\nNeed nrow(Dstruc) == ncol(Y)")
      if(ncol(Dstruc) != nfac) stop("Input 'Dstruc' must have 'nfac' columns.")
      if(any("logical" != c(apply(Dstruc, 1:2, class)))) stop("Input 'Dstruc' must be a matrix with logical (TRUE/FALSE) entries.")
    }
    
    # check 'Amodes' and 'Bmodes' and 'Dmodes' inputs
    if(any(const[1] == const.unimod) && !is.null(Amodes)){
      if(model == "parafac"){
        Amodes <- as.matrix(Amodes)
        if(nrow(Amodes) != 2L | ncol(Amodes) != nfac) stop("Input 'Amodes' must be a 2 x nfac matrix giving the\n minimum (row 1) and maximum (row 2) allowable mode for each factor.")
        Amodes <- matrix(as.integer(Amodes), nrow = 2L, ncol = nfac)
        if(any(Amodes < 1L)) stop("First row of 'Amodes' must contain integers greater than or equal to one.")
        if(any(Amodes > xdim[1])) stop("Second row of 'Amodes' must contain integers less than or equal to dim(X)[1].")
        if(any((Amodes[2,] - Amodes[1,]) < 0)) stop("Input 'Amodes' must satisfy:  Amodes[1,r] <= Amodes[2,r]")
      }
    }
    if(any(const[2] == const.unimod) && !is.null(Bmodes)){
      if(any(model == c("parafac", "parafac2"))){
        Bmodes <- as.matrix(Bmodes)
        if(nrow(Bmodes) != 2L | ncol(Bmodes) != nfac) stop("Input 'Bmodes' must be a 2 x nfac matrix giving the\n minimum (row 1) and maximum (row 2) allowable mode for each factor.")
        Bmodes <- matrix(as.integer(Bmodes), nrow = 2L, ncol = nfac)
        if(any(Bmodes < 1L)) stop("First row of 'Bmodes' must contain integers greater than or equal to one.")
        if(any(Bmodes > xdim[2])) stop("Second row of 'Bmodes' must contain integers less than or equal to dim(X)[2].")
        if(any((Bmodes[2,] - Bmodes[1,]) < 0)) stop("Input 'Bmodes' must satisfy:  Bmodes[1,r] <= Bmodes[2,r]")
      }
    }
    if(any(const[4] == const.unimod) && !is.null(Dmodes)){
      if(any(model == c("parafac", "parafac2"))){
        Dmodes <- as.matrix(Dmodes)
        if(nrow(Dmodes) != 2L | ncol(Dmodes) != nfac) stop("Input 'Dmodes' must be a 2 x nfac matrix giving the\n minimum (row 1) and maximum (row 2) allowable mode for each factor.")
        Dmodes <- matrix(as.integer(Dmodes), nrow = 2L, ncol = nfac)
        if(any(Dmodes < 1L)) stop("First row of 'Dmodes' must contain integers greater than or equal to one.")
        if(any(Dmodes > ydim[2])) stop("Second row of 'Dmodes' must contain integers less than or equal to dim(Y)[2].")
        if(any((Dmodes[2,] - Dmodes[1,]) < 0)) stop("Input 'Dmodes' must satisfy:  Dmodes[1,r] <= Dmodes[2,r]")
      }
    }
    
    ### check 'parallel' and 'cl' inputs
    if(parallel && !any(class(cl)=="cluster")){
      stop("Input 'cl' must be cluster (created by makeCluster) when parallel=TRUE \n  See examples in documentation:  ?mcr")
    }
    
    ### check 'maxit' and 'ctol' inputs
    maxit <- as.integer(maxit[1])
    if(maxit < 1L) stop("Input 'maxit' must be positive integer")
    ctol <- as.numeric(ctol[1])
    if(ctol < 0L) stop("Input 'ctol' must be positive numeric")
    
    ### check weights
    if(!is.null(weights)){
      weights <- as.numeric(weights)
      if(length(weights) != xdim[3]) stop("Input weights must have length equal to nrow(Y).")
      if(any(weights <= 0)) stop("Input weights must be positive scalars.")
      wsqrt <- sqrt(weights)
      for(k in 1:xdim[3]){
        X[,,k] <- wsqrt[k] * X[,,k]
        Y[k,] <- wsqrt[k] * Y[k,]
        if(!is.null(Cfixed)) Cfixed[k,] <-  wsqrt[k] * Cfixed[k,]
      }
      ssx.now <- ssx
      ssy.now <- ssy
      ssx <- sum(X^2)
      ssy <- sum(Y^2)
    }
    
    ### projection Y
    projY <- NULL
    if(xdim[3] > prod(xdim[-3])){
      xsvd <- svd(matrix(aperm(X, perm = c(3,1,2)), nrow = xdim[3], ncol = xdim[1] * xdim[2]), nv = 0)
      projY <- xsvd$u %*% crossprod(xsvd$u, Y)
    } 
    
    if(!parallel && verbose) pbar <- txtProgressBar(min = 0, max = nstart, style = 3)
    
    ### fit model
    if(model == "parafac"){
      
      # check if parallel
      if(parallel){
        nstartlist <- vector("list", nstart)
        nstartlist[1:nstart] <- nfac
        modlist <- parLapply(cl = cl, X = nstartlist, fun = "mcr_parafac",
                             dataX = X, dataY = Y, alpha = alpha,
                             const = const, control = control, ssx = ssx, ssy = ssy,
                             Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                             Astart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                             Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc, Dstruc = Dstruc,
                             Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes, Dmodes = Dmodes,
                             maxit = maxit, ctol = ctol, projY = projY, backfit = backfit)
      } else {
        
        if(output[1] == "best"){
          
          # fit model
          mod <- mcr_parafac(dataX = X, dataY = Y, nfac = nfac, alpha = alpha,
                             const = const, control = control, ssx = ssx, ssy = ssy,
                             Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                             Astart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                             Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc, Dstruc = Dstruc,
                             Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes, Dmodes = Dmodes,
                             maxit = maxit, ctol = ctol, projY = projY, backfit = backfit)
          if(verbose) setTxtProgressBar(pbar, 1)
          
          # fit model (additional random starts)
          if(nstart > 1L){
            for(i in 2:nstart){
              newmod <- mcr_parafac(dataX = X, dataY = Y, nfac = nfac, alpha = alpha,
                                    const = const, control = control, ssx = ssx, ssy = ssy, 
                                    Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                                    Astart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                                    Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc, Dstruc = Dstruc,
                                    Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes, Dmodes = Dmodes,
                                    maxit = maxit, ctol = ctol, projY = projY, backfit = backfit)
              if(verbose) setTxtProgressBar(pbar, i)
              if(newmod$LOSS < mod$LOSS) mod <- newmod
            } # end for(i in 2:nstart)
            rm(newmod)
          } # end if(nstart > 1L)
          
        } else {
          
          modlist <- vector("list",nstart)
          for(i in 1:nstart){
            modlist[[i]] <- mcr_parafac(dataX = X, dataY = Y, nfac = nfac, alpha = alpha,
                                        const = const, control = control, ssx = ssx, ssy = ssy, 
                                        Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                                        Astart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                                        Astruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc, Dstruc = Dstruc,
                                        Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes, Dmodes = Dmodes,
                                        maxit = maxit, ctol = ctol, projY = projY, backfit = backfit)
            if(verbose) setTxtProgressBar(pbar, i)
          } # end for(i in 1:nstart)
          
        } # end if(output[1] == "best")
        
      } # end if(parallel)
      
    } else if(model == "parafac2"){
      
      # check if parallel
      if(parallel){
        nstartlist <- vector("list", nstart)
        nstartlist[1:nstart] <- nfac
        modlist <- parLapply(cl = cl, X = nstartlist, fun = "mcr_parafac2",
                             dataX = X, dataY = Y, alpha = alpha,
                             const = const, control = control, ssx = ssx, ssy = ssy, 
                             Gfixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                             Gstart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                             Gstruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc, Dstruc = Dstruc,
                             Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes, Dmodes = Dmodes,
                             maxit = maxit, ctol = ctol, projY = projY, xsvd = xsvd, backfit = backfit)
      } else {
        
        if(output[1] == "best"){
          
          # fit model
          mod <- mcr_parafac2(dataX = X, dataY = Y, nfac = nfac, alpha = alpha,
                              const = const, control = control, ssx = ssx, ssy = ssy, 
                              Gfixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                              Gstart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                              Gstruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc, Dstruc = Dstruc,
                              Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes, Dmodes = Dmodes,
                              maxit = maxit, ctol = ctol, projY = projY, xsvd = xsvd, backfit = backfit)
          if(verbose) setTxtProgressBar(pbar, 1)
          
          # fit model (additional random starts)
          if(nstart > 1L){
            for(i in 2:nstart){
              newmod <- mcr_parafac2(dataX = X, dataY = Y, nfac = nfac, alpha = alpha,
                                     const = const, control = control, ssx = ssx, ssy = ssy, 
                                     Gfixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                                     Gstart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                                     Gstruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc, Dstruc = Dstruc,
                                     Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes, Dmodes = Dmodes,
                                     maxit = maxit, ctol = ctol, projY = projY, xsvd = xsvd, backfit = backfit)
              if(verbose) setTxtProgressBar(pbar, i)
              if(newmod$LOSS < mod$LOSS) mod <- newmod
            } # end for(i in 2:nstart)
            rm(newmod)
          } # end if(nstart > 1L)
          
        } else {
          
          modlist <- vector("list",nstart)
          for(i in 1:nstart){
            modlist[[i]] <- mcr_parafac2(dataX = X, dataY = Y, nfac = nfac, alpha = alpha,
                                         const = const, control = control, ssx = ssx, ssy = ssy, 
                                         Gfixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                                         Gstart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                                         Gstruc = Astruc, Bstruc = Bstruc, Cstruc = Cstruc, Dstruc = Dstruc,
                                         Amodes = Amodes, Bmodes = Bmodes, Cmodes = Cmodes, Dmodes = Dmodes,
                                         maxit = maxit, ctol = ctol, projY = projY, xsvd = xsvd, backfit = backfit)
            if(verbose) setTxtProgressBar(pbar, i)
          } # end for(i in 1:nstart)
          
        } # end if(output[1] == "best")
        
      } # end if(parallel)
      
    } else if(model == "tucker"){
      
      # check if parallel
      if(parallel){
        nstartlist <- vector("list", nstart)
        nstartlist[1:nstart] <- nfac
        modlist <- parLapply(cl = cl, X = nstartlist, fun = "mcr_tucker",
                             dataX = X, dataY = Y, alpha = alpha,
                             ssx = ssx, ssy = ssy, maxit = maxit, ctol = ctol, 
                             Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                             Astart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                             projY = projY)
      } else {
        
        if(output[1] == "best"){
          
          # fit model
          mod <- mcr_tucker(dataX = X, dataY = Y, nfac = nfac, alpha = alpha,
                            ssx = ssx, ssy = ssy, maxit = maxit, ctol = ctol, 
                            Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                            Astart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                            projY = projY)
          if(verbose) setTxtProgressBar(pbar, 1)
          
          # fit model (additional random starts)
          if(nstart > 1L){
            for(i in 2:nstart){
              newmod <- mcr_tucker(dataX = X, dataY = Y, nfac = nfac, alpha = alpha,
                                   ssx = ssx, ssy = ssy, maxit = maxit, ctol = ctol, 
                                   Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                                   Astart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                                   projY = projY)
              if(verbose) setTxtProgressBar(pbar, i)
              if(newmod$LOSS < mod$LOSS) mod <- newmod
            } # end for(i in 2:nstart)
            rm(newmod)
          } # end if(nstart > 1L)
          
        } else {
          
          modlist <- vector("list",nstart)
          for(i in 1:nstart){
            modlist[[i]] <- mcr_tucker(dataX = X, dataY = Y, nfac = nfac, alpha = alpha,
                                       ssx = ssx, ssy = ssy, maxit = maxit, ctol = ctol, 
                                       Afixed = Afixed, Bfixed = Bfixed, Cfixed = Cfixed, Dfixed = Dfixed,
                                       Astart = Astart, Bstart = Bstart, Cstart = Cstart, Dstart = Dstart,
                                       projY = projY)
            if(verbose) setTxtProgressBar(pbar, i)
          } # end for(i in 1:nstart)
          
        } # end if(output[1] == "best")
        
      } # end if(parallel)
      
    } # end if(model == "parafac")
    
    if(!parallel && verbose) close(pbar)
    
    ### output results
    if(output[1] == "best"){
      
      ## get best model (for parallel)
      if(parallel){
        LOSS <- sapply(modlist,function(x) x$LOSS)
        mod <- modlist[[which.min(LOSS)]]
      }
      
      ## check for weights
      if(!is.null(weights)){
        
        # insert weights
        mod$weights <- weights
        
        # retransform values
        for(k in 1:xdim[3]) {
          X[,,k] <- X[,,k] / wsqrt[k]
          Y[k,] <- Y[k,] / wsqrt[k]
          mod$C[k,] <- mod$C[k,] / wsqrt[k]
        }
        
        # fitted values for X
        if(model == "parafac"){
          Xhat <- array(mod$A %*% t(krprod(mod$C, mod$B)), dim = xdim)
        } else if(model == "parafac2"){
          Xhat <- array(0, dim = xdim)
          for(k in 1:xdim[3]) Xhat[,,k] <- mod$A[[k]] %*% (diag(nfac) * mod$C[k,]) %*% t(mod$B)
        } else if(model == "tucker"){
          Ga <- matrix(mod$G, nrow = nfac[1], ncol = prod(nfac[-1]))
          Xhat <- array(mod$A %*% Ga %*% t(kronecker(mod$C, mod$B)), dim = xdim)
        }
        
        # fitted values for Y
        Yhat <- tcrossprod(mod$C, mod$D)
        
        # get R-squared
        RsqX <- 1 - sum((X - Xhat)^2) / ssx.now
        RsqY <- 1 - sum((Y - Yhat)^2) / ssy.now
        mod$Rsq <- c(X = RsqX, Y = RsqY)
        
      } # end if(!is.null(weights))
      
      # define class and return
      class(mod) <- "mcr"
      return(mod)
      
    } # end if(output[1] == "best")
    
    # or output all?
    if(is.null(weights)){
      
      mcrfun <- function(mod){
        class(mod) <- "mcr"
        mod
      } # end mcrfun
      
    } else {
      
      mcrfun <- function(mod){
        
        # insert weights
        mod$weights <- weights
        
        # retransform values
        for(k in 1:xdim[3]) {
          X[,,k] <- X[,,k] / wsqrt[k]
          Y[k,] <- Y[k,] / wsqrt[k]
          mod$C[k,] <- mod$C[k,] / wsqrt[k]
        }
        
        # fitted values for X
        if(model == "parafac"){
          Xhat <- array(mod$A %*% t(krprod(mod$C, mod$B)), dim = xdim)
        } else if(model == "parafac2"){
          Xhat <- array(0, dim = xdim)
          for(k in 1:xdim[3]) Xhat[,,k] <- mod$A[[k]] %*% (diag(nfac) * mod$C[k,]) %*% t(mod$B)
        } else if(model == "tucker"){
          Ga <- matrix(mod$G, nrow = nfac[1], ncol = prod(nfac[-1]))
          Xhat <- array(mod$A %*% Ga %*% t(kronecker(mod$C, mod$B)), dim = xdim)
        }
        
        # fitted values for Y
        Yhat <- tcrossprod(mod$C, mod$D)
        
        # get R-squared
        RsqX <- 1 - sum((X - Xhat)^2) / ssx.now
        RsqY <- 1 - sum((Y - Yhat)^2) / ssy.now
        mod$Rsq <- c(X = RsqX, Y = RsqY)
        
        # define class and return
        class(mod) <- "mcr"
        mod
        
      } # end mcrfun
      
    } # end if(is.null(weights))
    
    modlist <- lapply(modlist, mcrfun)
    return(modlist)
    
  } # end mcr