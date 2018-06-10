print.parafac <- 
  function(x, ...){
    
    nway <- ifelse(is.null(x$D), 3, 4)
    nfac <- ncol(x$B)
    fixedID <- which(x$fixed)
    if(length(fixedID) > 0) x$const[fixedID] <- "fixed"
    strucID <- which(x$struc)
    if(length(strucID) > 0) x$const[strucID] <- paste(x$const[strucID], "struct", sep = "+")
    cvec <- data.frame(t(x$const), row.names="")
    names(cvec) <- LETTERS[1:nway]
    cat(paste0("\n",nway,"-way Parafac with ",nfac," factors"),"\n")
    cat("\nConstraints:\n")
    print(cvec)
    cat("\nFit Information:\n")
    cat("  SSE =", x$SSE, "\n")
    cat("  R^2 =", x$Rsq, "\n")
    cat("  GCV =", x$GCV, "\n")
    cat("  EDF =", sum(x$edf), "\n ")
    cat("\nConverged: ",(x$cflag == 0)," (",x$iter," iterations)","\n ",sep="")
    
  }