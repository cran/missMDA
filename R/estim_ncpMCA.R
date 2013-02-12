estim_ncpMCA <- function(don,ncp.min=0,ncp.max=5,method.cv="gcv",nbsim=100,pNA=0.05,threshold=1e-4){

#### Debut tab.disjonctif.NA
tab.disjonctif.NA<-function (tab) {
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
        x[indNA,]<-NA 
	  if ((ncol(tab) != 1) & (levels(moda)[1] %in% c(1:nlevels(moda), "n", "N", "y", "Y"))) dimnames(x) <- list(row.names(tab), paste(nom, levels(moda), sep = "."))
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
#### Fin tab.disjonctif.NA

########## Debut programme principal

method.cv <- tolower(method.cv)
auxi = NULL
for (j in 1:ncol(don)) if (is.numeric(don[,j])) auxi = c(auxi,colnames(don)[j])
if (!is.null(auxi)) stop(paste("\nAll variables are not categorical, the following ones are numeric: ", auxi))
vrai.tab=tab.disjonctif.NA(don)

if (method.cv=="gcv"){
p=ncol(don)
n=nrow(don)
crit <- NULL
tabX <- tab.disjonctif.NA(don)
if (is.null(ncp.max)) ncp.max <- ncol(tabX)-2
ncp.max <- min(nrow(tabX)-2,ncol(tabX)-1,ncp.max)
pquali <- ncol(tabX)
if (ncp.min == 0) crit = mean((tabX - rep(colMeans(tabX, na.rm = TRUE), each = nrow(tabX)))^2, na.rm = TRUE)
    for (q in max(ncp.min,1):ncp.max) {
        Z <- imputeMCA(don,ncp=q)$tab.disj
        Z <- scale(Z)*sqrt(nrow(Z)/(nrow(Z)-1))/sqrt(ncol(Z))
        ponder <- 1-apply(Z/nrow(Z), 2, sum)
        Z <- sweep(Z,2,sqrt(ponder),FUN="*")
        Z <- scale(Z,scale=FALSE)
        res.pca = PCA(Z, scale.unit = FALSE, graph = FALSE, ncp = max(q, 2))
        rec = reconst(res.pca, ncp = q)
        crit=c(crit,mean(( (n*(p-pquali))*(tabX-rec)/ (n*(p-pquali)- q*(n+p-pquali-q)))^2,na.rm=T))
    }
  if (any(diff(crit)>0)) { ncp = which(diff(crit)>0)[1]
  } else ncp <- which.min(crit)
 names(crit) <- c(ncp.min:ncp.max)
  return(list(ncp = as.integer(ncp+ncp.min-1),criterion=crit))
}

if (method.cv=="cv"){
res = matrix(NA,ncp.max-ncp.min+1,nbsim)

for (sim in 1:nbsim){
 donNA <- as.matrix(don)
 donNA[sample(1:(nrow(donNA)*ncol(donNA)),round(pNA*nrow(donNA)*ncol(donNA),0))] <- NA
 for (i in 1:ncol(don)) donNA[,i]=as.factor(as.character(donNA[,i]))

 for (nbaxes in ncp.min:ncp.max){
  tab.disj.comp <- imputeMCA(as.data.frame(donNA),ncp=nbaxes,threshold=threshold)$tab.disj
  res[nbaxes-ncp.min+1,sim] <- sum((tab.disj.comp-vrai.tab)^2,na.rm=TRUE)
 }
}
crit=apply(res,1,mean)
 names(crit) <- c(ncp.min:ncp.max)
result = list(ncp = as.integer(which.min(crit)+ncp.min-1),criterion=crit)
return(result)
}
}
