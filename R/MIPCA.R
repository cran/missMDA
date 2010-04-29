MIPCA <- function(X,ncp=2,scale=TRUE,method="Regularized",threshold=1e-4,nboot=100){

## Initialization
  method <- tolower(method)
  missing <- which(is.na(X))
  impute.data <- imputePCA(X,scale=scale,ncp=ncp,method=method,threshold=threshold)$completeObs
  reference <- PCA(impute.data,scale=scale,graph=FALSE,ncp=ncp)
  rec <- reconst(reference,ncp)
  rec.pca <- as.matrix(X)
  rec.pca[missing] <- rec[missing]
  resid <- rec.pca-rec #residuals are centred
  if (scale) resid <- sweep(resid,2,apply(rec,2,sd),FUN="/")  ###
  sigma <- sqrt((sum((resid[-missing])^2))/ (nrow(X)*ncol(X)-(length(missing)+ncol(X)+ncp*(nrow(X)-1+ncol(X)-ncp))))

  rownames(rec.pca) <- rownames(X)
  res.MI <- array(NA,dim=c(nrow(X),ncol(X),nboot))

for(i in 1:nboot){
### Sampling variability
 resid.star <- matrix(rnorm(nrow(X)*ncol(X),0,sigma),ncol=ncol(X))
 if (scale) resid.star <- sweep(resid.star,2,apply(rec,2,sd),FUN="*") ###
 Xstar <- rec+resid.star
## 2 rows to add some NA values
 missing2 <- sample(1:(nrow(X)*ncol(X)),length(missing))
 Xstar[missing2] <- NA
 acpboot <- PCA(imputePCA(Xstar,scale=scale,ncp=ncp,method=method,threshold=threshold)$completeObs,scale=scale,ncp=ncp,graph=FALSE)

###Drawing from the predictive distribution
 residstar2 <- matrix(rnorm(nrow(X)*ncol(X),0,sigma),ncol=ncol(X))
 if (scale) residstar2 <- sweep(residstar2,2,apply(rec,2,sd),FUN="*")
 rec.pca[missing] <- (reconst(acpboot,ncp)+residstar2)[missing]
 res.MI[,,i] <- rec.pca
}
 result=list(res.imputePCA=impute.data,res.MI=res.MI,call=list(X=X,ncp=ncp,missing=missing,nboot=nboot,scale=scale))
 class(result) <- c("MIPCA", "list")
 return(result)
}
