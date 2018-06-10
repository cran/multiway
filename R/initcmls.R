initcmls <-
  function(nobs, nfac = 1, const = "uncons", struc = NULL, 
           df = 10, degree = 3, intercept = TRUE,
           mode.range = NULL, standardize = TRUE){
    # Randomly Initialize a Constrained Matrix (for cmls)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: Apr 19, 2018
    
    
    ######******######   CONSTRAINT OPTIONS   ######******######
    
    # __________________________________________________________
    # uncons = unconstrained
    # nonneg = non-negative
    # period = periodicity
    # pernon = periodicity and non-negativity
    # __________________________________________________________
    # smooth = smoothness
    # smonon = smoothness and non-negative
    # smoper = smoothness and periodic
    # smpeno = smoothness and periodic and non-negative
    # __________________________________________________________
    # orthog = orthogonal
    # ortnon = orthogonal and non-negative
    # ortsmo = orthogonal and smooth
    # orsmpe = orthogonal and smooth and periodic
    # __________________________________________________________
    # moninc = monotonic increasing
    # monnon = monotonic increasing and non-negative
    # monsmo = monotonic increasing and smooth
    # mosmno = monotonic increasing and smooth and non-negative
    # __________________________________________________________
    # unimod = unimodal
    # uninon = unimodal and non-negative
    # uniper = unimodal and periodic
    # unpeno = unimodal and periodic and non-negative
    # __________________________________________________________
    # unismo = unimodal and smooth
    # unsmno = unimodal and smooth and and non-negative
    # unsmpe = unimodal and smooth and periodic
    # unsmpn = unimodal and smooth and periodic and non-negative
    
    
    ######******######   INITITAL CHECKS   ######******######
    
    # check 'nobs' input
    nobs <- as.integer(nobs[1])
    if(nobs < 1L) stop("Input 'nobs' must be a positive integer.")
    
    # check 'nfac' input
    nfac <- as.integer(nfac[1])
    if(nfac < 1L) stop("Input 'nfac' must be a positive integer.")
    
    # check 'const' input
    if(is.na(const) | is.null(const) | missing(const)) const <- "uncons"
    const <- as.character(const[1])
    const.types <- c("uncons", "nonneg", "period", "pernon",
                     "smooth", "smonon", "smoper", "smpeno",
                     "orthog", "ortnon", "ortsmo", "orsmpe",
                     "moninc", "monnon", "monsmo", "mosmno",
                     "unimod", "uninon", "uniper", "unpeno", 
                     "unismo", "unsmno", "unsmpe", "unsmpn")
    if(!any(const == const.types)) stop("Invalid 'const' input.")
    
    # check struc
    if(!is.null(struc)){
      struc <- as.matrix(struc)
      if(nrow(struc) != nobs | ncol(struc) != nfac) stop("Input 'struc' must satisfy:  nrow(struc) == nobs  &&  ncol(struc) == nfac")
      schck <- max(abs( ((struc == TRUE) + (struc == FALSE)) - 1 ))
      if(schck > 0) stop("Input 'struc' must be a matrix of logicals (TRUE/FALSE).")
      if(any(const == c("orthog", "ortnon", "ortsmo", "orsmpe"))){
        StS <- crossprod(struc)
        maxLT <- max(abs(StS[lower.tri(StS)]))
        if(maxLT > 0) stop("When using orthogonality constraints, the input 'struc' must\n  contain orthogonal columns (i.e., one or less TRUE in each column).")
      }
    }
    
    # check inputs relevant to smoothness updates
    const.smooth <- c("smooth", "smonon", "smoper", "smpeno",
                      "ortsmo", "orsmpe", "monsmo", "mosmno", 
                      "unismo", "unsmno", "unsmpe", "unsmpn")
    if(any(const == const.smooth)){
      
      # check df
      df <- as.integer(df[1])
      
      # check degree
      degree <- as.integer(degree[1])
      if(degree < 1L) stop("Input 'degree' should be greater than 0.") 
      
      # check intercept
      intercept <- intercept[1]
      if(!is.logical(intercept)) stop("Input 'intercept' must be a logical (TRUE/FALSE).")
      
      # check df (relative to degree and intercept)
      mindf <- degree + 1 - !intercept
      if(any(const == c("monsmo", "mosmno")) && intercept) mindf <- mindf + 1L
      if(df < mindf){
        df <- mindf
        warning("Input 'df' is too small!\n  Resetting 'df' to minimum possible value = ", mindf)
      } 
      
    } # end if(any(const == const.smooth))
    
    # check inputs relevant to unimodality updates
    const.unimod <- c("unimod", "uninon", "uniper", "unpeno", 
                      "unismo", "unsmno", "unsmpe", "unsmpn")
    if(any(const == const.unimod)){
      periodic <- ifelse(any(const == c("uniper", "unpeno", "unsmpe", "unsmpn")), TRUE, FALSE)
      if(is.null(mode.range)){
        if(is.null(struc)){
          if(periodic){
            mode.range <- matrix(c(2, nobs - 1L), nrow = 2, ncol = nfac)
          } else {
            mode.range <- matrix(c(1, nobs), nrow = 2, ncol = nfac)
          }
        } else {
          mode.range <- apply(struc, 2, function(x) range(which(x)))
          if(periodic){
            mode.range[1,] <- mode.range[1,] + 1L
            mode.range[2,] <- mode.range[2,] - 1L
          }
        } # end if(is.null(struc))
      } else {
        mode.range <- as.matrix(mode.range)
        if(nrow(mode.range) != 2L | ncol(mode.range) != nfac) stop("Input 'mode.range' must be a 2 x p matrix giving the\n minimum (row 1) and maximum (row 2) allowable mode for column row of B.")
        mode.range <- matrix(as.integer(mode.range), nrow = 2L, ncol = nfac)
        if(any(mode.range[1,] < 1L)) stop("First row of 'mode.range' must contain integers greater than or equal to one.")
        if(any(mode.range[2,] > nobs)) stop("Second row of 'mode.range' must contain integers less than or equal to m.")
        if(any((mode.range[2,] - mode.range[1,]) < 0)) stop("Input 'mode.range' must satisfy:  mode.range[1,j] <= mode.range[2,j]")
        if(!is.null(struc)){
          for(jj in 1:nfac){
            mseq <- seq(mode.range[1,jj], mode.range[2,jj])
            if(any(!struc[mseq,jj])) stop("Inputs 'mode.range' and 'struc' are incompatible.\n Need struc[,j] to have TRUE all values between mode.range[1,j] and mode.range[2,j].")
            if(periodic){
              if(mode.range[1,jj] == 1L) stop("Unimodal and periodic functions must have a mode range between 2 and m-1.")
              if(mode.range[2,jj] == nobs) stop("Unimodal and periodic functions must have a mode range between 2 and m-1.")
            }
          }
        } # end if(!is.null(struc))
      } # end if(is.null(mode.range))
    } # end if(any(const == const.unimod))
    
    
    ######******######   UNCONSTRAINED   ######******######
    
    # unconstrained: unstructured and structured
    if(const == "uncons"){
      B <- matrix(rnorm(nobs * nfac), nrow = nobs, ncol = nfac)
      if(!is.null(struc)) B <- B * struc
      if(standardize) B <- scale(B, center = FALSE, scale = sqrt(colMeans(B^2)))
      return(matrix(B, nrow = nobs, ncol = nfac))
    } # end if(const == "uncons")
    
    
    ######******######   NON-NEGATIVE   ######******######
    
    # non-negative: unstructured and structured
    if(const == "nonneg"){
      B <- matrix(runif(nobs * nfac) + 0.5, nrow = nobs, ncol = nfac)
      if(!is.null(struc)) B <- B * struc
      if(standardize) B <- scale(B, center = FALSE, scale = sqrt(colMeans(B^2)))
      return(matrix(B, nrow = nobs, ncol = nfac))
    } # end if(const == "nonneg")
    
    
    ######******######   PERIODIC   ######******######
    
    # period and pernon: unstructured and structured
    if(any(const == c("period", "pernon"))){
      nonneg <- ifelse(const == "pernon", TRUE, FALSE)
      z <- seq(0, 1, length.out = nobs)
      if(is.null(struc)){
        Z <- MsplineBasis(z, df = nobs, intercept = TRUE, periodic = TRUE)
        if(nonneg){
          coefs <- matrix(abs(rnorm(nobs * nfac)), nrow = nobs, ncol = nfac)
        } else {
          coefs <- matrix(rnorm(nobs * nfac), nrow = nobs, ncol = nfac)
        }
        B <- Z %*% coefs
      } else {
        mindf <- 4L
        indx <- c(1, nobs)
        B <- matrix(0, nrow = nobs, ncol = nfac)
        for(jj in 1:nfac){
          zjj <- z[struc[,jj]]
          ndf <- max(round(length(zjj) / nobs), mindf)
          Z <- matrix(0, nrow = nobs, ncol = ndf)
          Z[struc[,jj],] <- MsplineBasis(x = zjj, df = ndf, intercept = TRUE, periodic = TRUE)
          if(struc[1,jj] != struc[nobs,jj]){
            Fid <- which(struc[indx,jj])
            Z[indx[Fid],] <- Z[indx[-Fid],]
          }
          if(nonneg){
            B[,jj] <- Z %*% abs(rnorm(ndf))
          } else {
            B[,jj] <- Z %*% rnorm(ndf)
          }
        } # end for(jj in 1:nfac)
      } # end if(is.null(struc))
      if(standardize) B <- scale(B, center = FALSE, scale = sqrt(colMeans(B^2)))
      return(matrix(B, nrow = nobs, ncol = nfac))
    } # end if(any(const == c("period", "pernon"))
    
    
    ######******######   SMOOTHNESS   ######******######
    
    # smooth and smonon: unstructured and structured
    if(any(const == c("smooth", "smonon"))){
      nonneg <- ifelse(const == "smonon", TRUE, FALSE)
      z <- seq(0, 1, length.out = nobs)
      if(is.null(struc)){
        Z <- MsplineBasis(z, df = df, degree = degree, intercept = intercept)
        if(nonneg){
          coefs <- matrix(abs(rnorm(df * nfac)), nrow = df, ncol = nfac)
        } else {
          coefs <- matrix(rnorm(df * nfac), nrow = df, ncol = nfac)
        }
        B <- Z %*% coefs
      } else {
        mindf <- degree + 1 - !intercept
        B <- matrix(0, nrow = nobs, ncol = nfac)
        for(jj in 1:nfac){
          zjj <- z[struc[,jj]]
          ndf <- max(round(df * length(zjj) / nobs), mindf)
          Z <- matrix(0, nrow = nobs, ncol = ndf)
          Z[struc[,jj],] <- MsplineBasis(x = zjj, df = ndf, degree = degree, intercept = intercept)
          if(nonneg){
            B[,jj] <- Z %*% abs(rnorm(ndf))
          } else {
            B[,jj] <- Z %*% rnorm(ndf)
          }
        } # end for(jj in 1:nfac)
      } # end if(is.null(struc))
      if(standardize) B <- scale(B, center = FALSE, scale = sqrt(colMeans(B^2)))
      return(matrix(B, nrow = nobs, ncol = nfac))
    } # end if(any(const == c("smooth", "smonon")))
    
    # smoper and smpeno: unstructured and structured
    if(any(const == c("smoper", "smpeno"))){
      nonneg <- ifelse(const == "smpeno", TRUE, FALSE)
      z <- seq(0, 1, length.out = nobs)
      if(is.null(struc)){
        Z <- MsplineBasis(z, df = df, degree = degree, 
                          intercept = intercept, periodic = TRUE)
        if(nonneg){
          coefs <- matrix(abs(rnorm(df * nfac)), nrow = df, ncol = nfac)
        } else {
          coefs <- matrix(rnorm(df * nfac), nrow = df, ncol = nfac)
        }
        B <- Z %*% coefs
      } else {
        mindf <- degree + 1 - !intercept
        indx <- c(1, nobs)
        B <- matrix(0, nrow = nobs, ncol = nfac)
        for(jj in 1:nfac){
          zjj <- z[struc[,jj]]
          ndf <- max(round(df * length(zjj) / nobs), mindf)
          Z <- matrix(0, nrow = nobs, ncol = ndf)
          Z[struc[,jj],] <- MsplineBasis(x = zjj, df = ndf, degree = degree, 
                                         intercept = intercept, periodic = TRUE)
          if(struc[1,jj] != struc[nobs,jj]){
            Fid <- which(struc[indx,jj])
            Z[indx[Fid],] <- Z[indx[-Fid],]
          }
          if(nonneg){
            B[,jj] <- Z %*% abs(rnorm(ndf))
          } else {
            B[,jj] <- Z %*% rnorm(ndf)
          }
        } # end for(jj in 1:nfac)
      } # end if(is.null(struc))
      if(standardize) B <- scale(B, center = FALSE, scale = sqrt(colMeans(B^2)))
      return(matrix(B, nrow = nobs, ncol = nfac))
    } # end if(any(const == c("smoper", "smpeno")))
    
    
    ######******######   ORTHOGONALITY   ######******######
    
    # orthog: unconstrained and structured
    if(const == "orthog"){
      if(is.null(struc)){
        B <- matrix(rnorm(nobs * nfac), nrow = nobs, ncol = nfac)
        B <- svd(B, nv = 0)$u
      } else {
        B <- matrix(rnorm(nobs * nfac), nrow = nobs, ncol = nfac) * struc
        B <- scale(B, center = FALSE, scale = sqrt(colSums(B^2)))
      }
      return(matrix(B, nrow = nobs, ncol = nfac))
    } # end if(const == "orthog")
    
    # ortnon: unconstrained and structured
    if(const == "ortnon"){
      B <- matrix(0, nrow = nobs, ncol = nfac)
      while(any(colSums(B) <= 0)){
        B <- matrix(0, nrow = nobs, ncol = nfac)
        if(is.null(struc)){
          for(ii in 1:nobs) B[ii,sample.int(nfac, size = 1)] <- abs(rnorm(1))
        } else {
          for(ii in 1:nobs) B[ii,which(struc[ii,])] <- abs(rnorm(1))
        } 
      }
      B <- scale(B, center = FALSE, scale = sqrt(colSums(B^2)))
      return(matrix(B, nrow = nobs, ncol = nfac))
    } # end if(const == "ortnon")
    
    # ortsmo and orsmpe: unconstrained and structured
    if(any(const == c("ortsmo", "orsmpe"))){
      periodic <- ifelse(const == "orsmpe", TRUE, FALSE)
      z <- seq(0, 1, length.out = nobs)
      if(is.null(struc)){
        Z <- MsplineBasis(z, df = df, degree = degree, 
                          intercept = intercept, periodic = periodic)
        B <- Z %*% matrix(rnorm(df * nfac), nrow = df, ncol = nfac)
        R <- svd(matrix(rnorm(nfac^2), nrow = nfac, ncol = nfac), nv = 0)$u
        B <- svd(B, nv = 0)$u %*% R
      } else {
        mindf <- degree + 1 - !intercept
        indx <- c(1, nobs)
        B <- matrix(0, nrow = nobs, ncol = nfac)
        for(jj in 1:nfac){
          zjj <- z[struc[,jj]]
          ndf <- max(round(df * length(zjj) / nobs), mindf)
          Z <- matrix(0, nrow = nobs, ncol = ndf)
          Z[struc[,jj],] <- MsplineBasis(x = zjj, df = ndf, degree = degree, 
                                         intercept = intercept, periodic = periodic)
          if(periodic && struc[1,jj] != struc[nobs,jj]){
            Fid <- which(struc[indx,jj])
            Z[indx[Fid],] <- Z[indx[-Fid],]
          }
          B[,jj] <- Z %*% rnorm(ndf)
        } # end for(jj in 1:nfac)
        B <- scale(B, center = FALSE, scale = sqrt(colSums(B^2)))
      } # end if(is.null(struc))
      return(matrix(B, nrow = nobs, ncol = nfac))
    } # end if(any(const == c("ortsmo", "orsmpe")))
    
    
    ######******######   MONOTONICITY   ######******######
    
    # moninc or monnon: unstructured or structured
    if(any(const == c("moninc", "monnon"))){
      nonneg <- ifelse(const == "monnon", TRUE, FALSE)
      z <- seq(0, 1, length.out = nobs)
      Z <- IsplineBasis(z, df = nobs)
      B <- matrix(0, nrow = nobs, ncol = nfac)
      for(jj in 1:nfac) {
        coefs <- rnorm(nobs)^4
        coefs <- coefs / sum(coefs)
        B[,jj] <- Z %*% coefs + ifelse(nonneg, runif(1), runif(1) - 0.5)
      }
      # structural constraints
      if(!is.null(struc)){
        B <- B * struc
        for(jj in 1:nfac){
          difjj <- struc[2:nobs,jj] - struc[1:(nobs-1),jj]
          ndiff <- sum(abs(difjj))
          if(ndiff == 1L){
            if(sum(difjj) < 0){
              indx <- which(difjj < 0)
              if(nonneg){
                B[1:indx,jj] <- 0
              } else {
                if(B[indx,jj] > B[indx+1,jj]) B[1:indx,jj] <- B[1:indx,jj] - B[indx,jj] - runif(1)
              } # end if(nonneg)
            } else {
              indx <- which(difjj > 0)
              B[(indx+1):nobs,jj] <- B[(indx+1):nobs,jj] - B[indx+1,jj] + runif(1)
              if(nonneg) B[(indx+1):nobs,jj] <- pmax(B[(indx+1):nobs,jj], 0)
            } # end if(sum(difjj) < 0)
          } else if(ndiff == 2L) {
            # T to F
            indx <- which(difjj < 0)
            if(nonneg){
              B[1:indx,jj] <- 0
            } else {
              if(B[indx,jj] > B[indx+1,jj]) B[1:indx,jj] <- B[1:indx,jj] - B[indx,jj] - runif(1)
            } # end if(nonneg)
            # F to T
            indx <- which(difjj > 0)
            B[(indx+1):nobs,jj] <- B[(indx+1):nobs,jj] - B[indx+1,jj] + runif(1)
            if(nonneg) B[(indx+1):nobs,jj] <- pmax(B[(indx+1):nobs,jj], 0)
          } # end if(ndiff == 1L)
        } # end for(jj in 1:nfac)
      } # end if(!is.null(struc))
      if(standardize) {
        scales <- sqrt(colMeans(B^2))
        if(any(scales == 0)) scales[scales == 0] <- 1
        B <- scale(B, center = FALSE, scale = scales)
      }
      return(matrix(B, nrow = nobs, ncol = nfac))
    } # end if(any(const == c("moninc", "monnon")))
    
    # monsmo or mosmno: unstructured or structured
    if(any(const == c("monsmo", "mosmno"))){
      nonneg <- ifelse(const == "mosmno", TRUE, FALSE)
      z <- seq(0, 1, length.out = nobs)
      B <- matrix(0, nrow = nobs, ncol = nfac)
      if(is.null(struc)){
        Z <- IsplineBasis(z, df = df, degree = degree, intercept = intercept)
        for(jj in 1:nfac){
          coefs <- rnorm(df)^4
          B[,jj] <- Z %*% coefs + ifelse(nonneg, runif(1), runif(1) - 0.5)
        }
      } else {
        # structural constraints
        mindf <- 4L
        B <- matrix(0, nrow = nobs, ncol = nfac)
        for(jj in 1:nfac){
          zjj <- z[struc[,jj]]
          ndf <- max(round(df * length(zjj) / nobs), mindf)
          Z <- matrix(0, nrow = nobs, ncol = ndf)
          Z[struc[,jj],] <- IsplineBasis(x = zjj, df = ndf, degree = degree, intercept = intercept)
          coefs <- rnorm(ndf)^4
          coefs <- coefs / sum(coefs)
          B[,jj] <- Z %*% coefs + ifelse(nonneg, runif(1), runif(1) - 0.5) * struc[,jj]
          difjj <- struc[2:nobs,jj] - struc[1:(nobs-1),jj]
          ndiff <- sum(abs(difjj))
          if(ndiff == 1L){
            if(sum(difjj) < 0){
              indx <- which(difjj < 0)
              if(nonneg){
                B[1:indx,jj] <- 0
              } else {
                if(B[indx,jj] > B[indx+1,jj]) B[1:indx,jj] <- B[1:indx,jj] - B[indx,jj] - runif(1, max = 0.1)
              } # end if(nonneg)
            } else {
              indx <- which(difjj > 0)
              B[(indx+1):nobs,jj] <- B[(indx+1):nobs,jj] - B[indx+1,jj] + runif(1, max = 0.1)
              if(nonneg) B[(indx+1):nobs,jj] <- pmax(B[(indx+1):nobs,jj], 0)
            } # end if(sum(difjj) < 0)
          } else if(ndiff == 2L) {
            # T to F
            indx <- which(difjj < 0)
            if(nonneg){
              B[1:indx,jj] <- 0
            } else {
              if(B[indx,jj] > B[indx+1,jj]) B[1:indx,jj] <- B[1:indx,jj] - B[indx,jj] - runif(1, max = 0.1)
            } # end if(nonneg)
            # F to T
            indx <- which(difjj > 0)
            B[(indx+1):nobs,jj] <- B[(indx+1):nobs,jj] - B[indx+1,jj] + runif(1, max = 0.1)
            if(nonneg) B[(indx+1):nobs,jj] <- pmax(B[(indx+1):nobs,jj], 0)
          } # end if(ndiff == 1L)
        } # end for(jj in 1:nfac)
      } # end if(is.null(struc))
      if(standardize) {
        scales <- sqrt(colMeans(B^2))
        if(any(scales == 0)) scales[scales == 0] <- 1
        B <- scale(B, center = FALSE, scale = scales)
      }
      return(matrix(B, nrow = nobs, ncol = nfac))
    } # end if(any(const == c("monsmo", "mosmno")))
    
    
    ######******######   UNIMODALITY   ######******######
    
    # unimod or uninon: unstructured or structured
    if(any(const == const.unimod)){
      B <- matrix(0, nrow = nobs, ncol = nfac)
      while(any(colSums(abs(B)) <= 0)){
        if(is.null(struc)) z <- seq(0, 1, length.out = nobs)
        B <- matrix(0, nrow = nobs, ncol = nfac)
        for(jj in 1:nfac){
          mseq <- seq(mode.range[1,jj], mode.range[2,jj])
          if(!is.null(struc)){
            wstruc <- which(struc[,jj])
            ntrue <- length(wstruc)
            z <- seq(0, 1, length.out = ntrue)
          }
          if(length(mseq) == 1){
            modejj <- mseq
          } else {
            modejj <- sample(mseq, size = 1)
          }
          if(is.null(struc)){
            modeindx <- modejj
            theta <- z[modeindx]
          } else {
            theta <- z[modejj - min(wstruc) + 1L]
          }
          if(theta == 0){
            alpha <- 1
            beta <- runif(1, min = 1, max = 20)
          } else if (theta == 1){
            alpha <- runif(1, min = 1, max = 20)
            beta <- 1
          } else {
            alpha <- runif(1, min = 1, max = 20)
            beta <- (alpha * (1 - theta) + 2 * theta - 1) / theta
          }
          if(is.null(struc)){
            B[,jj] <- dbeta(z, alpha, beta)
          } else {
            B[struc[,jj],jj] <- dbeta(z, alpha, beta)
          }
        } # end for(jj in 1:nfac)
      } # end while(any(colSums(abs(B)) <= 0))
      
      if(standardize) B <- scale(B, center = FALSE)
      return(matrix(B, nrow = nobs, ncol = nfac))
    } # end if(any(const == const.unimod))
    
  }