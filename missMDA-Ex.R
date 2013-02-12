pkgname <- "missMDA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('missMDA')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("MIPCA")
### * MIPCA

flush(stderr()); flush(stdout())

### Name: MIPCA
### Title: Multiple Imputation with PCA
### Aliases: MIPCA
### Keywords: multivariate

### ** Examples

data(orange)
## First the number of components has to be chosen 
##   (for the reconstruction step)
## nb <- estim_ncpPCA(orange,ncp.max=5) ## Time consuming, nb = 2

## Multiple Imputation
resMI <- MIPCA(orange,ncp=2)

## Visualization on the PCA map
plot(resMI)



cleanEx()
nameEx("estim_ncpMCA")
### * estim_ncpMCA

flush(stderr()); flush(stdout())

### Name: estim_ncpMCA
### Title: Estimate the number of dimensions for the Multiple
###   Correspondence Analysis by cross-validation
### Aliases: estim_ncpMCA
### Keywords: multivariate

### ** Examples

## Not run: 
##D data(vnf)
##D result <- estim_ncpMCA(vnf,ncp.min=0, ncp.max=3)
## End(Not run)



cleanEx()
nameEx("estim_ncpPCA")
### * estim_ncpPCA

flush(stderr()); flush(stdout())

### Name: estim_ncpPCA
### Title: Estimate the number of dimensions for the Principal Component
###   Analysis by cross-validation
### Aliases: estim_ncpPCA
### Keywords: multivariate

### ** Examples

## Not run: 
##D data(orange)
##D nb <- estim_ncpPCA(orange,ncp.min=0,ncp.max=4) 
## End(Not run)



cleanEx()
nameEx("imputeFAMD")
### * imputeFAMD

flush(stderr()); flush(stdout())

### Name: imputeFAMD
### Title: Impute dataset with mixed data
### Aliases: imputeFAMD
### Keywords: models multivariate

### ** Examples

data(ozone)
res.comp <- imputeFAMD(ozone, ncp=3)
res.afdm <- AFDM(ozone,tab.comp=res.comp)



cleanEx()
nameEx("imputeMCA")
### * imputeMCA

flush(stderr()); flush(stdout())

### Name: imputeMCA
### Title: Impute missing values in categorical variables with Multiple
###   Correspondence Analysis
### Aliases: imputeMCA
### Keywords: models multivariate

### ** Examples

data(vnf)
## First the number of components has to be chosen 
##   (for the reconstruction step)
## nb <- estim_ncpMCA(vnf,ncp.max=5) ## Time-consuming, nb = 4

## Impute indicator matrix and perform a MCA
tab.disj.impute <- imputeMCA(vnf, ncp=4)
res.mca <- MCA(vnf,tab.disj=tab.disj.impute$tab.disj)



cleanEx()
nameEx("imputeMFA")
### * imputeMFA

flush(stderr()); flush(stdout())

### Name: imputeMFA
### Title: Impute dataset with MFA
### Aliases: imputeMFA
### Keywords: models multivariate

### ** Examples

data(orange)
res.comp <- imputeMFA(orange,group=c(5,3),type=rep("s",2),ncp=2)
## Note that MFA is performed on the completed matrix
res.mfa <- MFA(res.comp$completeObs,group=c(5,3),type=rep("s",2))

## Not run: 
##D data(vnf)
##D res.comp <- imputeMFA(vnf,group=c(6,5,3),type=c("n","n","n"),ncp=2)
##D res.mfa <- MFA(vnf,group=c(6,5,3),type=c("n","n","n"),tab.comp=res.comp)
## End(Not run)



cleanEx()
nameEx("imputePCA")
### * imputePCA

flush(stderr()); flush(stdout())

### Name: imputePCA
### Title: Impute dataset with PCA
### Aliases: imputePCA
### Keywords: models multivariate

### ** Examples

data(orange)
## First the number of components has to be chosen 
##   (for the reconstruction step)
## nb <- estim_ncpPCA(orange,ncp.max=5) ## Time consuming, nb = 2

## Imputation
res.comp <- imputePCA(orange,ncp=2)

## A PCA can be performed
res.pca <- PCA(res.comp$completeObs)



cleanEx()
nameEx("orange")
### * orange

flush(stderr()); flush(stdout())

### Name: orange
### Title: Sensory description of 12 orange juices by 8 attributes.
### Aliases: orange
### Keywords: datasets

### ** Examples

data(orange)
## Not run: 
##D nb <- estim_ncpPCA(orange,ncp.min=0,ncp.max=5,method.cv="Kfold",nbsim=20,pNA=0.05)
##D res.comp <- imputePCA(orange,ncp=nb$ncp)
##D res.pca <- PCA(res.comp$completeObs)
##D resMI <- MIPCA(orange,ncp=nb$ncp)
##D plot(resMI)
## End(Not run)



cleanEx()
nameEx("ozone")
### * ozone

flush(stderr()); flush(stdout())

### Name: ozone
### Title: Daily measurements of meteorological variables and ozone
###   concentration
### Aliases: ozone
### Keywords: datasets

### ** Examples

data(ozone)
res.comp <- imputeFAMD(ozone, ncp=3)
res.afdm <- AFDM(ozone,tab.comp=res.comp)



cleanEx()
nameEx("plot.MIPCA")
### * plot.MIPCA

flush(stderr()); flush(stdout())

### Name: plot.MIPCA
### Title: Plot the graphs for the Multiple Imputation in PCA
### Aliases: plot.MIPCA
### Keywords: dplot

### ** Examples

data(orange)
## nb <- estim_ncpPCA(orange,ncp.max=5) ## Time consuming, nb = 2
resMI <- MIPCA(orange,ncp=2)
plot(resMI)



cleanEx()
nameEx("vnf")
### * vnf

flush(stderr()); flush(stdout())

### Name: vnf
### Title: Questionnaire done by 1232 individuals who answered 14 questions
### Aliases: vnf
### Keywords: datasets

### ** Examples

data(vnf)
tab.disj.impute <- imputeMCA(vnf, ncp=4)$tab.disj
res.mca <- MCA(vnf,tab.disj=tab.disj.impute)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
