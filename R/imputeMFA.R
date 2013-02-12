imputeMFA <- function (X, group, ncp = 2, type=rep("s",length(group)), method="Regularized",row.w=NULL,coeff.ridge=1,threshold = 1e-6,seed = NULL,maxiter=1000,...){

impute <- function (X, group, ncp = 2, type=rep("s",length(group)), method=NULL,threshold = 1e-6,seed = NULL,maxiter=1000,row.w=NULL,coeff.ridge=1,...){
    moy.p <- function(V, poids) {
        res <- sum(V * poids,na.rm=TRUE)/sum(poids[!is.na(V)])
    }
    ec <- function(V, poids) {
        res <- sqrt(sum(V^2 * poids,na.rm=TRUE)/sum(poids[!is.na(V)]))
    }
    tab.disjonctif.NA <- function(tab) {
        tab <- as.data.frame(tab)
        modalite.disjonctif <- function(i) {
            moda <- tab[, i]
            nom <- names(tab)[i]
            n <- length(moda)
            moda <- as.factor(moda)
            x <- matrix(0, n, length(levels(moda)))
            ind <- (1:n) + n * (unclass(moda) - 1)
            indNA <- which(is.na(ind))
            x[(1:n) + n * (unclass(moda) - 1)] <- 1
            x[indNA, ] <- NA
            if ((ncol(tab) != 1) & (levels(moda)[1] %in% c(1:nlevels(moda), 
                "n", "N", "y", "Y"))) 
                dimnames(x) <- list(row.names(tab), paste(nom, 
                  levels(moda), sep = "."))
            else dimnames(x) <- list(row.names(tab), levels(moda))
            return(x)
        }
        if (ncol(tab) == 1) 
            res <- modalite.disjonctif(1)
        else {
            res <- lapply(1:ncol(tab), modalite.disjonctif)
            res <- as.matrix(data.frame(res, check.names = FALSE))
        }
        return(res)
    }

find.category <- function (X,tabdisj){
  nbdummy <- rep(1,ncol(X))
  is.quali <- which(!unlist(lapply(X,is.numeric)))
  nbdummy[is.quali] <- unlist(lapply(X[,is.quali,drop=FALSE],nlevels))
  vec = c(0,cumsum(nbdummy))
  Xres <- X
  for (i in 1:ncol(X)){
    if (i %in% is.quali) Xres[,i] <- as.factor(levels(X[,i])[apply(tabdisj[,(vec[i]+1):vec[i+1]],1,which.max)])
	else Xres[,i] <- tabdisj[,vec[i]+1]
  }
  return(Xres)
}

############# main program
   if (!is.null(seed)) set.seed(seed)
   if ("n"%in%type){
      niveau = NULL
      for (j in 1:ncol(X)){
        if (!is.numeric(X[,j])) niveau = c(niveau,levels(X[,j]))
      }
      for (j in 1:ncol(X)) {
        if (!is.numeric(X[,j])){
          if (sum(niveau%in%levels(X[,j]))!=nlevels(X[,j])) levels(X[,j]) = paste(colnames(X)[j],levels(X[,j]),sep="_")
        }
      }
    }
  if (is.null(row.w)) row.w=rep(1/nrow(X),nrow(X))
	group.mod=group
    Xhat <- matrix(0,nrow(X),0)
    Xhat2 <- matrix(0,nrow(X),0)
	MM <- vector(mode = "list", length = length(group)) ## vecteur de moyenne ou de marge selon quanti ou quali
	ET <- vector(mode = "list", length = length(group))
	tab.disj.comp <- vector(mode = "list", length = length(group))
	ponderation <- rep(1,length(group))
	for (g in 1:length(group)){
	    if (g==1) aux.base <- X[,1:group[1],drop=FALSE]
		else aux.base <- X[,(cumsum(group)[g-1]+1):cumsum(group)[g],drop=FALSE]
        if (type[g] == "s") {
          Xhat2 <- cbind.data.frame(Xhat2, aux.base)
		  MM[[g]] <- apply(as.data.frame(aux.base), 2, moy.p, row.w)
          aux.base <- as.matrix(sweep(as.data.frame(aux.base), 2, MM[[g]], FUN = "-"))
          ET[[g]] <- apply(as.data.frame(aux.base), 2, ec, row.w)
#          ET[[g]][ET[[g]] <= 1e-08] <- 1
          aux.base <- as.matrix(sweep(as.data.frame(aux.base), 2, ET[[g]], FUN = "/"))
          missing <- which(is.na(as.matrix(aux.base)))
          if (any(is.na(aux.base))) aux.base[missing] <- 0
          ponderation[g] <- svd.triplet(aux.base,ncp=1,row.w=row.w)$vs[1]
		  Xhat <- cbind.data.frame(Xhat, aux.base/ponderation[g])
          if (!is.null(seed)&(length(missing)!=0)) Xhat[missing,] <- rnorm(length(missing)) ## random initialization
        }
        if (type[g] == "c") {
          Xhat2 <- cbind.data.frame(Xhat2, aux.base)
		  MM[[g]] <- apply(as.data.frame(aux.base), 2, moy.p,row.w)
          aux.base <- as.matrix(sweep(as.data.frame(aux.base), 2, MM[[g]], FUN = "-"))
          missing <- which(is.na(as.matrix(aux.base)))
          if (any(is.na(aux.base))) aux.base[missing] <- 0       
          ponderation[g] = svd.triplet(aux.base,ncp=1,row.w=row.w)$vs[1]		 
		  Xhat <- cbind.data.frame(Xhat, aux.base/ponderation[g])
          if (!is.null(seed)&(length(missing)!=0)) Xhat[missing,] <- rnorm(length(missing)) ## random initialization   
        }
        if (type[g] == "n") {
		   tab.disj = tab.disjonctif.prop(aux.base,seed,row.w=row.w)
		   tab.disj.comp[[g]] = tab.disj
		   group.mod[g] <- ncol(tab.disj)
           MM[[g]] = apply(tab.disj, 2, moy.p,row.w)/ncol(aux.base)
           Z = sweep(tab.disj, 2, apply(tab.disj, 2, moy.p,row.w), FUN = "/")
           Z = sweep(Z, 2,apply(Z,2,moy.p,row.w),FUN="-")
	       Zscale = sweep(Z, 2, sqrt(MM[[g]]), FUN = "*")
           ponderation[g] <- svd.triplet(Zscale,row.w=row.w)$vs[1]
		   Xhat <- cbind.data.frame(Xhat, Zscale/ponderation[g])
           Xhat2 <- cbind.data.frame(Xhat2, as.data.frame(tab.disjonctif.NA(aux.base)))
        }
	}
   ind.var <- vector(mode = "list", length = length(group))
   ind.var[[1]] <- 1:group.mod[1]
   for (g in 2:length(group)) ind.var[[g]] <- (cumsum(group.mod)[g-1]+1):cumsum(group.mod)[g]
   recon <- Xhat <- as.matrix(Xhat)
   if (ncp>=min(nrow(Xhat)-2,ncol(Xhat)-1)) stop("ncp is too large")
   ncp <- min(ncp,ncol(X),nrow(X)-1)
   missing <- which(is.na(as.matrix(Xhat2)))
   nb.iter <- 1
   old <- Inf
   while (nb.iter > 0) {
     Xhat[missing] <- recon[missing]
 	 for (g in 1:length(group)){
       if (g==1) aux.base <- Xhat[,1:group.mod[1],drop=FALSE]
	   else aux.base <- Xhat[,(cumsum(group.mod)[g-1]+1):cumsum(group.mod)[g],drop=FALSE]
 	   aux.base <- aux.base*ponderation[g]
       if (type[g] == "s") {
		  aux.base <- sweep(aux.base,2,ET[[g]],FUN="*")
		  aux.base <- sweep(aux.base,2,MM[[g]],FUN="+")
		  MM[[g]] <- apply(aux.base, 2, moy.p,row.w)
		  aux.base <- sweep(aux.base,2,MM[[g]], FUN="-")
		  ET[[g]] <- apply(aux.base, 2, ec,row.w)
		  aux.base <- sweep(aux.base,2,ET[[g]], FUN="/")
		  ponderation[g] = svd.triplet(aux.base,ncp=1,row.w=row.w)$vs[1]
       }
       if (type[g] == "c") {
		  aux.base <- sweep(aux.base,2,MM[[g]],FUN="+")
		  MM[[g]] <- apply(aux.base, 2, moy.p,row.w)
		  aux.base <- sweep(aux.base,2,MM[[g]], FUN="-")
		  ponderation[g] = svd.triplet(aux.base,ncp=1,row.w=row.w)$vs[1]
       }
       if (type[g] == "n") {
          tab.disj = sweep(aux.base,2,sqrt(MM[[g]]),FUN="/") + matrix(1,nrow(aux.base),ncol(aux.base)) 
          tab.disj = sweep(tab.disj,2,apply(tab.disj.comp[[g]],2,moy.p,row.w),FUN="*")
          tab.disj.comp[[g]] = tab.disj
          MM[[g]] = apply(tab.disj, 2, moy.p,row.w)/ncol(aux.base)
		  if (any(MM[[g]]<0)) stop(paste("The algorithm fails to converge. Choose a number of components (ncp) less or equal than ",ncp-1," or a number of iteratsions (maxiter) less or equal than ",maxiter-1,sep=""))
          Z = sweep(tab.disj, 2, apply(tab.disj, 2, moy.p,row.w), FUN = "/")
          Z = sweep(Z, 2,apply(Z,2,moy.p,row.w),FUN="-")
	      aux.base = sweep(Z, 2, sqrt(MM[[g]]), FUN = "*")
          ponderation[g] <- svd.triplet(aux.base,row.w=row.w)$vs[1]
       }
	   if (g==1) Xhat[,1:group.mod[1]] <- aux.base / ponderation[g]
	   else Xhat[,(cumsum(group.mod)[g-1]+1):cumsum(group.mod)[g]] <- aux.base / ponderation[g]
	}
	 svd.res <- svd.triplet(Xhat,row.w=row.w)
     sigma2 <- mean(svd.res$vs[-(1:ncp)]^2)
     sigma2 <- min(sigma2*coeff.ridge,svd.res$vs[ncp+1]^2)
	 
     if (method=="em") sigma2 <-0
     lambda.shrinked=(svd.res$vs[1:ncp]^2-sigma2)/svd.res$vs[1:ncp]
     if (ncp==1) recon=tcrossprod(sweep(svd.res$U[,1,drop=FALSE],1,row.w,FUN="*")*lambda.shrinked,svd.res$V[,1,drop=FALSE])
     else recon=tcrossprod(sweep(sweep(svd.res$U[,1:ncp],1,row.w,FUN="*"),2,lambda.shrinked,FUN="*"),svd.res$V[,1:ncp])
     recon <- sweep(recon,1,row.w,FUN="/")

	 diff <- Xhat-recon
	 diff[missing] <- 0
     objective <- sum(sweep(diff^2,1,row.w,FUN="*"))
     criterion <- abs(1 - objective/old)
     old <- objective
     nb.iter <- nb.iter + 1
     if (!is.nan(criterion)) {
       if ((criterion < threshold) && (nb.iter > 5))  nb.iter <- 0
       if ((objective < threshold) && (nb.iter > 5))  nb.iter <- 0
     }
     if (nb.iter>maxiter) {
       nb.iter <- 0
       warning(paste("Stopped after ",maxiter," iterations"))
     }
   }
### end of the loop
	for (g in 1:length(group)){
	    if (g==1) aux.base <- Xhat[,1:group.mod[1],drop=FALSE]
		else aux.base <- Xhat[,(cumsum(group.mod)[g-1]+1):cumsum(group.mod)[g],drop=FALSE]
		aux.base <- aux.base*ponderation[g]
        if (type[g] == "s") {
		  aux.base <- sweep(aux.base,2,ET[[g]],FUN="*")
		  aux.base <- sweep(aux.base,2,MM[[g]],FUN="+")
        }
        if (type[g] == "c")  aux.base <- sweep(aux.base,2,MM[[g]],FUN="+")
        if (type[g] == "n") {
          tab.disj = sweep(aux.base,2,sqrt(MM[[g]]),FUN="/") + matrix(1,nrow(aux.base),ncol(aux.base)) 
##          aux.base = sweep(sweep(tab.disj,1,row.w,FUN="*"),2,apply(tab.disj.comp[[g]],2,sum),FUN="*")
          aux.base = sweep(tab.disj,2,apply(tab.disj.comp[[g]],2,moy.p,row.w),FUN="*")
		  }
		if (g==1) Xhat[,1:group.mod[1]] <- aux.base
		else Xhat[,(cumsum(group.mod)[g-1]+1):cumsum(group.mod)[g]] <- aux.base
	}
   completeObs <- as.matrix(Xhat2)
   completeObs[missing] <- Xhat[missing]

   result <- list()
   result$tab.disj <- completeObs
   result$completeObs <- find.category(X,completeObs)
   result$call$group.mod = group.mod
   result$call$ind.var = ind.var
   return(result) 
}

#### Main program
 obj=Inf
 method <- tolower(method)
 if (is.null(row.w)) row.w = rep(1,nrow(X))/nrow(X)
# for (i in 1:nb.init){
# i=1
  if (!any(is.na(X))) stop("no missing values in X, this function is not useful. Perform MFA on X.")
  res.impute <- impute(X, group=group,ncp=ncp, type=type, method=method, threshold = threshold,seed=seed,maxiter=maxiter,row.w=row.w,coeff.ridge=coeff.ridge)
return(res.impute)
#  if (mean((res.impute$recon[!is.na(X)]-X[!is.na(X)])^2) < obj){
#    res <- res.impute
#    obj <- mean((res.impute$recon[!is.na(X)]-X[!is.na(X)])^2)
#  }
# }
#return(res)
}

