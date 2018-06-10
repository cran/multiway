print.mcr <-
  function(x, ...){
    if(x$model == "tucker"){
      nfac <- c(ncol(x$A), ncol(x$B), ncol(x$C))
      nfac <- paste(nfac, collapse = " x ")
    } else {
      nfac <- ncol(x$B)
    }
    cap <- function(xx){
      ss <- strsplit(xx, "", fixed = TRUE)[[1]]
      paste0(c(toupper(ss[1]), ss[-1]), collapse = "")
    }
    fixedID <- which(x$fixed)
    if(length(fixedID) > 0) x$const[fixedID] <- "fixed"
    strucID <- which(x$struc)
    if(length(strucID) > 0) x$const[strucID] <- "struct"
    cvec <- data.frame(t(x$const), row.names="")
    names(cvec) <- LETTERS[1:4]
    cat(paste0("\nMCR-",cap(x$model)," with ", nfac, " factors"), "\n")
    cat("\nConstraints:\n")
    print(cvec)
    cat("\nFit Information:\n")
    cat("  R^2.X =", x$Rsq[1], "\n")
    cat("  R^2.Y =", x$Rsq[2], "\n")
    cat("   LOSS =", x$LOSS, "\n")
    cat("  alpha =", x$alpha, "\n")
    cat("\nConverged: ", (x$cflag == 0), " (", x$iter, " iterations)", "\n ", sep = "")
  }