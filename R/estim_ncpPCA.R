estim_ncpPCA <- function(X,ncp.min=0,ncp.max=5,method="Regularized",scale=TRUE,method.cv="loo",nbsim=100,pNA=0.05,threshold=1e-4){

## method = "em" or "Regularized"
## method.cv = "loo" (for leave-one-out) or "Kfold" (a percentage of pNA missing values is added and nbsim are done)

method <- tolower(method)
method.cv <- tolower(method.cv)
auxi = NULL
for (j in 1:ncol(X)) if (!is.numeric(X[,j])) auxi = c(auxi,colnames(X)[j])
if (!is.null(auxi)) stop(paste("\nThe following variables are not quantitative: ", auxi))
ncp.max <- min(ncp.max,ncol(X)-1,nrow(X)-2)
res <- NULL

if (method.cv=="loo"){
 for (nbaxes in ncp.min:ncp.max){
   Xhat <- X
   for (i in 1:nrow(X)){
    for (j in 1:ncol(X)){
     if (!is.na(X[i,j])){
      XNA <- as.matrix(X)
      XNA[i,j] <- NA
      if (nbaxes==0) Xhat[i,j] <- mean(XNA[,j],na.rm=TRUE)
      else Xhat[i,j] <- imputePCA(XNA,ncp=nbaxes,threshold=threshold,method=method,scale=scale)$completeObs[i,j]
    }
   }
  }
  res <- c(res,mean((Xhat-X)^2,na.rm=TRUE))
 }
 names(res) <- c(ncp.min:ncp.max)
 result = list(ncp = which.min(res)+ncp.min-1,criterion=res)
}


if (method.cv=="kfold"){
  res <- matrix(NA,ncp.max-ncp.min+1,nbsim)
  for (sim in 1:nbsim){
   XNA <- as.matrix(X)
   XNA[sample(1:(nrow(XNA)*ncol(XNA)),round(pNA*nrow(XNA)*ncol(XNA),0))] <- NA
   for (nbaxes in ncp.min:ncp.max){
    if (nbaxes==0) {
       Xhat <- XNA
       for (j in 1:ncol(X)) Xhat[,j] <- replace(XNA[,j],is.na(XNA[, j]),mean(XNA[,j],na.rm=TRUE))
    } else Xhat <- imputePCA(XNA,ncp=nbaxes,threshold=threshold,method=method,scale=scale)$completeObs  
   res[nbaxes-ncp.min+1,sim] <- sum((Xhat-X)^2,na.rm=TRUE)
  }
 }
 resu <- apply(res,1,mean)
 result <- list(ncp = which.min(resu)+ncp.min-1,criterion=resu)
}
return(result)
}
