imputeAFDM <- function (X, ncp = 2, method="Regularized",row.w=NULL,coeff.ridge=1,threshold = 1e-6,seed = NULL,nb.init=1,maxiter=1000,...){

  type = rep("s",ncol(X))
  type[!unlist(lapply(X,is.numeric))]="n"
  res <- imputeMFA(X=X,group=rep(1,ncol(X)),type=type, ncp = ncp, method=method,row.w=row.w,coeff.ridge=coeff.ridge,threshold = threshold,seed = seed,nb.init=nb.init,maxiter=maxiter)
  res$call <- NULL
  return(res)
}

