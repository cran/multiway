print.parafac <- 
  function(x, ...){
    
    conlab <- c("none", "orthogonal", "nonnegative", "unimodal", "monotonic", "periodic", "smooth", "fixed", "structure")
    nway <- ifelse(is.null(x$D),3,4)
    nfac <- ncol(x$B)
    fixedID <- which(x$fixed)
    if(length(fixedID) > 0) x$const[fixedID] <- 7
    strucID <- which(x$struc)
    if(length(strucID) > 0) x$const[strucID] <- 8
    cvec <- data.frame(t(conlab[x$const + 1]), row.names="")
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