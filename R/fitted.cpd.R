fitted.cpd <- function(object, ...){
  nmode <- length(object$A)
  xdims <- sapply(object$A, nrow)
  Z <- krprod(object$A[[3]], object$A[[2]])
  if(nmode > 3L) for(n in 4:nmode) Z <- krprod(object$A[[n]], Z)
  array(tcrossprod(object$A[[1]], Z), dim = xdims)
}