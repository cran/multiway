print.tucker <- 
  function(x, ...){
    
    nway <- ifelse(is.null(x$D),3,4)
    nfac <- dim(x$G)
    cat(paste0("\n",nway,"-way Tucker with ",paste(nfac,collapse=" x ")," factors"),"\n")
    cat("\nFit Information:\n")
    cat("  SSE =", x$SSE, "\n")
    cat("  R^2 =", x$Rsq, "\n")
    cat("  GCV =", x$GCV, "\n")
    cat("  EDF =", sum(x$edf), "\n ")
    cat("\nConverged: ",(x$cflag == 0)," (",x$iter," iterations)","\n ",sep="")
    
  }