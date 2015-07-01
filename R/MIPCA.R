MIPCA <- function(X,ncp=2,scale=TRUE,method=c("Regularized","EM"),threshold=1e-4,nboot=100){

## Initialization
  method <- match.arg(method,c("Regularized","regularized","EM","em"),several.ok=T)[1]
  method <- tolower(method)
  missing <- which(is.na(X))
  inputeData <- imputePCA(X,scale=scale,ncp=ncp,method=method,threshold=threshold)$completeObs
  reference <- PCA(inputeData,scale.unit=scale,graph=FALSE,ncp=ncp)
  rec <- reconst(reference,ncp)
  rec.pca <- as.matrix(X)
  rec.pca[missing] <- rec[missing]
  resid <- rec.pca-rec #residuals are centred
  sdResid <- apply(inputeData,2,sd)
##  sdResid <- apply(rec,2,sd)
  if (scale) resid <- t(t(resid)/sdResid)  ###
  sigma <- sqrt((sum((resid[-missing])^2))/ (nrow(X)*ncol(X)-(length(missing)+ncol(X)+ncp*(nrow(X)-1+ncol(X)-ncp))))

  rownames(rec.pca) <- rownames(X)
  res.MI <- array(NA,dim=c(nrow(X),ncol(X),nboot))

for(i in 1:nboot){
### Sampling variability
 resid.star <- matrix(rnorm(nrow(X)*ncol(X),0,sigma),ncol=ncol(X))
 if (scale) resid.star <- t(t(resid.star)*sdResid) ###
 resid.star[missing]<- NA
 Xstar <- rec+resid.star-matrix(mean(resid.star,na.rm=TRUE),ncol=ncol(resid.star),nrow=nrow(resid.star))

# Xstar <- rec+resid.star
## 2 rows to add some NA values
# missing2 <- sample(1:(nrow(X)*ncol(X)),length(missing))
##missing2=missing
# Xstar[missing2] <- NA

 acpboot <- PCA(imputePCA(Xstar,scale=scale,ncp=ncp,method=method,threshold=threshold)$completeObs,scale.unit=scale,ncp=ncp,graph=FALSE)

###Drawing from the predictive distribution
 residstar2 <- matrix(rnorm(nrow(X)*ncol(X),0,sigma),ncol=ncol(X))
 if (scale) residstar2 <- t(t(residstar2)*sdResid)
 rec.pca[missing] <- (reconst(acpboot,ncp)+residstar2)[missing]
 res.MI[,,i] <- rec.pca
}
 dimnames(res.MI)=list(rownames(X),colnames(X),NULL)
 result=list(res.imputePCA=inputeData,res.MI=res.MI,call=list(X=X,ncp=ncp,missing=missing,nboot=nboot,scale=scale))
 class(result) <- c("MIPCA", "list")
 return(result)
}
