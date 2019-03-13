print.cpd <- 
  function(x, ...){
    
    nway <- length(x$A)
    nfac <- ncol(x$A[[1]])
    cat(paste0("\n",nway,"-way Canonical Polyadic Decomposition with ",nfac," factors"),"\n")
    cat("\nFit Information:\n")
    cat("  SSE =", x$SSE, "\n")
    cat("  R^2 =", x$Rsq, "\n")
    cat("  GCV =", x$GCV, "\n")
    cat("  EDF =", sum(x$edf), "\n ")
    cat("\nConverged: ",(x$cflag == 0)," (",x$iter," iterations)","\n ",sep="")
    
  }