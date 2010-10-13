imputeMFA <- function (X, group, ncp = 2, scale=TRUE, method="Regularized",threshold = 1e-6,seed = NULL,nb.init=1,maxiter=1000,...){

impute <- function (X, group, ncp = 4, scale=TRUE, method=NULL,threshold = 1e-6,seed = NULL,init=1,maxiter=1000,...){

   if (!is.null(seed)) set.seed(seed)
   X <- as.matrix(X)
   nbr <- nrow(X)
   nbc <- ncol(X)
   ncp <- min(ncp,nbc,nbr-1)
   Obs <- !is.na(X)
   missing <- which(is.na(X))
   moy.p <- apply(X, 2, mean,na.rm=TRUE)
   et <- apply(X, 2, sd,na.rm=TRUE)
   Xhat <- sweep(X, 2,moy.p,FUN="-")
   if (scale) Xhat <- sweep(Xhat, 2,et,FUN="/")
   if (any(is.na(X))) Xhat[missing] <- 0
   if (init>1) Xhat[missing] <- rnorm(length(missing)) ## random initialization
   recon <- Xhat
 ponderation = rep(nbr/svd(scale(Xhat[,1:group[1]],scale=FALSE),nu=1,nv=1)$d[1]^2,group[1])
       for (i in 2:length(group)) ponderation = c(ponderation,rep(nbr/svd(scale(Xhat[,(sum(group[1:(i-1)])+1):sum(group[1:i])],scale=FALSE),nu=1,nv=1)$d[1]^2,group[i]))
 Xhat <- sweep(Xhat, 2, sqrt(ponderation), FUN = "*")
 
   nb.iter <- 1
   old <- Inf
   while (nb.iter > 0) {
       Xhat[missing] <- recon[missing]
       Xhat <- sweep(Xhat, 2, sqrt(ponderation), FUN = "/")
       if (scale) Xhat <- sweep(Xhat,2,et, FUN="*")
       Xhat <- sweep(Xhat,2,moy.p, FUN="+")
       moy.p <- apply(Xhat, 2, mean,na.rm=TRUE)
       et <- apply(Xhat, 2, sd,na.rm=TRUE)
       Xhat=sweep(Xhat,2,moy.p, FUN="-")
       if (scale) Xhat=sweep(Xhat,2,et, FUN="/")
       ponderation = rep(nbr/svd(scale(Xhat[,1:group[1]],scale=FALSE),nu=1,nv=1)$d[1]^2,group[1])
       for (i in 2:length(group)) ponderation = c(ponderation,rep(nbr/svd(scale(Xhat[,(sum(group[1:(i-1)])+1):sum(group[1:i])],scale=FALSE),nu=1,nv=1)$d[1]^2,group[i]))
       Xhat = sweep(Xhat,2,sqrt(ponderation),FUN="*")

       svd.res <- svd.triplet(scale(Xhat,scale=FALSE))
       sigma2 <- mean(svd.res$vs[-(1:ncp)]^2)
       if (method=="em") sigma2 <-0

       if (ncp==1) {
         lambda.shrinked=(svd.res$vs[1]^2-sigma2)/sqrt(svd.res$vs[1]^2)
         recon=(svd.res$U[,1]%*%diag(lambda.shrinked[1],1)%*%(t(svd.res$V[,1])))
       } else {
         lambda.shrinked=diag((svd.res$vs[1:ncp]^2-sigma2)/sqrt(svd.res$vs[1:ncp]^2))
         recon=(svd.res$U[,1:ncp]%*%lambda.shrinked%*% (t(svd.res$V[,1:ncp])))
       }

       objective <- mean((Xhat[-missing]-recon[-missing])^2)
       criterion <- abs(1 - objective/old)
       old <- objective
       nb.iter <- nb.iter + 1
       if (!is.nan(criterion)) {
         if ((criterion < threshold) && (nb.iter > 5))  nb.iter <- 0
       }
       if (nb.iter>maxiter) {
         nb.iter <- 0
         warning(paste("Stopped after ",maxiter," iterations"))
       }
   }
   Xhat = sweep(Xhat,2,sqrt(ponderation),"/")
   if (scale) Xhat <- sweep(Xhat,2,et, FUN="*")
   Xhat <- sweep(Xhat,2,moy.p, FUN="+")
   completeObs <- X
   completeObs[missing] <- Xhat[missing]
   if (scale) recon <- sweep(recon,2,et, FUN="*")
   recon <- sweep(recon,2,moy.p, FUN="+")

   result <- list()
   result$completeObs <- completeObs
   result$objective <- objective
   result$recon <- recon
   return(result) 
}

#### Main program
 obj=Inf
 method <- tolower(method)
 if (ncp>=min(nrow(X)-2,ncol(X)-1)) stop("ncp is too large")
 for (i in 1:nb.init){
  if (!any(is.na(X))) return(X)
  res.impute=impute(X, group=group, ncp=ncp, scale=scale, method=method, threshold = threshold,seed=seed,init=i,maxiter=maxiter)
  if (mean((res.impute$recon[!is.na(X)]-X[!is.na(X)])^2) < obj){
    res <- res.impute
    obj <- mean((res.impute$recon[!is.na(X)]-X[!is.na(X)])^2)
  }
 }
return(res)
}

