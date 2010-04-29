imputeMCA <- function(don,ncp=2,threshold=1e-6,seed=NULL,maxiter=1000){   
 
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

tab.disjonctif.prop<-function (tab,seed=NULL) 
{
    tab <- as.data.frame(tab)
    modalite.disjonctif <- function(i) {
        moda <- tab[, i]
        nom <- names(tab)[i]
        n <- length(moda)
        moda <- as.factor(moda)
        x <- matrix(0, n, length(levels(moda)))
          ind<-(1:n) + n * (unclass(moda) - 1)
          indNA<-which(is.na(ind))
                
        x[(1:n) + n * (unclass(moda) - 1)] <- 1
        if (length(indNA)!=0){
          if (is.null(seed)) {
           x[indNA,]<- matrix(rep(apply(x,2,sum)/sum(x),each=length(indNA)),nrow=length(indNA))
          } else {
           for (k in 1:length(indNA)) {
            aux <- runif(length(levels(moda)))
            x[indNA[k],]=aux/sum(aux)
           }
          }
         }
          if ((ncol(tab) != 1) & (levels(moda)[1] %in% c(1:nlevels(moda),"n", "N", "y", "Y"))) 
            dimnames(x) <- list(row.names(tab), paste(nom, levels(moda),sep = "."))
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

########## Debut programme principal
if (ncp==0) return(tab.disjonctif.prop(don,NULL))

tab.disj.NA = tab.disjonctif.NA(don)
hidden = which(is.na(tab.disj.NA))
tab.disj.comp=tab.disjonctif.prop(don,seed)
tab.disj.rec.old=tab.disj.comp


D_I = diag(1/nrow(don),nrow(don))
inv_D_I = diag(nrow(don),nrow(don))
continue=TRUE
nbiter=0

while (continue){
  nbiter=nbiter+1
  M = apply(tab.disj.comp,2,sum)/(nrow(don)*ncol(don))
  Z=nrow(don)*sweep(tab.disj.comp,2,apply(tab.disj.comp,2,sum),FUN="/")
  Z=scale(Z,scale=FALSE)
  inv_racine_M = sqrt((nrow(don)*ncol(don))/apply(tab.disj.comp,2,sum))
  Zscale=sweep(Z,2,sqrt(M),FUN="*")

  svd.Zscale=svd.triplet(Zscale)
  moyeig=0
  if (ncp>0){
    if (nrow(don)>ncol(Zscale)) moyeig=mean(svd.Zscale$vs[-c(1:ncp,(length(svd.Zscale$vs)-ncol(don)+1):length(svd.Zscale$vs))]^2)
    else moyeig=mean(svd.Zscale$vs[-c(1:ncp)]^2)
  }
  if (ncp==1) eig.shrunk=diag(((svd.Zscale$vs[1]^2-moyeig)/svd.Zscale$vs[1]),1)
  else eig.shrunk=diag(((svd.Zscale$vs[1:ncp]^2-moyeig)/svd.Zscale$vs[1:ncp])) 
        
  rec=svd.Zscale$U[,1:ncp]%*%eig.shrunk%*% (t(svd.Zscale$V[,1:ncp]))
        
  tab.disj.rec = sweep(rec,2,inv_racine_M,FUN="*") + matrix(1,nrow(rec),ncol(rec)) 
  tab.disj.rec=D_I%*%sweep(tab.disj.rec,2,apply(tab.disj.comp,2,sum),FUN="*")
  relch=sum((tab.disj.rec[hidden] - tab.disj.rec.old[hidden])^2)
  tab.disj.rec.old=tab.disj.rec
  tab.disj.comp[hidden] = tab.disj.rec[hidden]
  continue=(relch > threshold)&(nbiter<maxiter)
  if (ncp==0) continue=FALSE
}

return(tab.disj.comp)
}

