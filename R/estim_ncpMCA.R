estim_ncpMCA <- function(don,ncp.min=0,ncp.max=5,nbsim=100,pNA=0.05,threshold=1e-4){

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
threshold=1e-4
auxi = NULL
for (j in 1:ncol(don)) if (is.numeric(don[,j])) auxi = c(auxi,colnames(don)[j])
if (!is.null(auxi)) stop(paste("\nAll variables are not categorical, the following ones are numeric: ", auxi))
vrai.tab=tab.disjonctif.NA(don)
res = matrix(NA,ncp.max-ncp.min+1,nbsim)

for (sim in 1:nbsim){
 donNA <- as.matrix(don)
 donNA[sample(1:(nrow(donNA)*ncol(donNA)),round(pNA*nrow(donNA)*ncol(donNA),0))] <- NA
 for (i in 1:ncol(don)) donNA[,i]=as.factor(as.character(donNA[,i]))

 for (nbaxes in ncp.min:ncp.max){
  tab.disj.comp <- imputeMCA(donNA,ncp=nbaxes,threshold=threshold)
  res[nbaxes-ncp.min+1,sim] <- sum((tab.disj.comp-vrai.tab)^2,na.rm=TRUE)
 }
}
tab.disj=apply(res,1,mean)

result = list(ncp = which.min(tab.disj)+ncp.min-1,tab.disj=tab.disj)
return(result)

}

