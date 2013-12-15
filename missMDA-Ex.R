pkgname <- "missMDA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('missMDA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("MIPCA")
### * MIPCA

flush(stderr()); flush(stdout())

### Name: MIPCA
### Title: Multiple Imputation with PCA
### Aliases: MIPCA
### Keywords: multivariate

### ** Examples

## Not run: 
##D data(orange)
##D ## First the number of components has to be chosen 
##D ##   (for the reconstruction step)
##D ## nb <- estim_ncpPCA(orange,ncp.max=5) ## Time consuming, nb = 2
##D 
##D ## Multiple Imputation
##D resMI <- MIPCA(orange,ncp=2)
##D 
##D ## Visualization on the PCA map
##D plot(resMI)
## End(Not run)



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
##D result <- estim_ncpMCA(vnf,ncp.min=0, ncp.max=5)
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
nameEx("gene")
### * gene

flush(stderr()); flush(stdout())

### Name: gene
### Title: Gene expression
### Aliases: gene
### Keywords: datasets

### ** Examples

data(gene)
res.impute <- imputeMFA(gene[,-1], group = c(76,356), 
    type = rep("s",2), ncp = 2) 
res.mfa <- MFA(cbind.data.frame(gene[,1], res.impute$completeObs), 
      group = c(1,76,356), type=c("n",rep("s",2)), 
	  name.group = c("WHO","CGH","expr"), num.group.sup = 1)
plot.MFA(res.mfa, habillage = 1, lab.ind = FALSE)
plot.MFA(res.mfa, habillage = "group", invisible = "ind", partial = "all")
plot.MFA(res.mfa, habillage = "group", lab.ind = FALSE, partial = "all")
plot.MFA(res.mfa, choix = "var", habillage = "group", lab.var = FALSE)
plot.MFA(res.mfa, choix = "group", habillage = "group")



cleanEx()
nameEx("geno")
### * geno

flush(stderr()); flush(stdout())

### Name: geno
### Title: Genotype-environment data set with missing values
### Aliases: geno
### Keywords: datasets

### ** Examples

data(geno)

res.ncp.gcv <- estim_ncpPCA(geno)
res.imp <- imputePCA(geno, ncp= res.ncp.gcv$ncp)
res.pca <- PCA(res.imp$completeObs)



cleanEx()
nameEx("imputeFAMD")
### * imputeFAMD

flush(stderr()); flush(stdout())

### Name: imputeFAMD
### Title: Impute mixed dataset
### Aliases: imputeFAMD
### Keywords: models multivariate

### ** Examples

data(ozone)
res.impute <- imputeFAMD(ozone, ncp=3) 
## The output can be used as an input of the FAMD function of the FactoMineR package 
##to perform the FAMD on the incomplete data ozone 
res.afdm <- FAMD(ozone,tab.comp=res.impute) 



cleanEx()
nameEx("imputeMCA")
### * imputeMCA

flush(stderr()); flush(stdout())

### Name: imputeMCA
### Title: Impute categorical dataset
### Aliases: imputeMCA
### Keywords: models multivariate

### ** Examples

data(vnf)
## First the number of components has to be chosen 
##   (for the reconstruction step)
## nb <- estim_ncpMCA(vnf,ncp.max=5) ## Time-consuming, nb = 4

## Impute the indicator matrix and perform a MCA
res.impute <- imputeMCA(vnf, ncp=4)

## The imputed indicator matrix can be used as an input of the MCA function of the
## FactoMineR package to perform the MCA on the incomplete data ozone 
res.mca <- MCA(vnf,tab.disj=res.impute$tab.disj) 



cleanEx()
nameEx("imputeMFA")
### * imputeMFA

flush(stderr()); flush(stdout())

### Name: imputeMFA
### Title: Impute dataset with variables structured into groups of
###   variables (groups of continuous or categorical variables)
### Aliases: imputeMFA
### Keywords: models multivariate

### ** Examples

data(orange)
## Impute the data and perform a MFA
## with groups of continuous variables only
res.impute <- imputeMFA(orange, group=c(5,3), type=rep("s",2),ncp=2) 
res.mfa <- MFA(res.impute$completeObs,group=c(5,3),type=rep("s",2)) 

## Not run: 
##D data(vnf)
##D ## Impute the indicator matrix and perform a MFA 
##D ## with groups of categorical variables only
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
##   (for the imputation step)
## nb <- estim_ncpPCA(orange,ncp.max=5) ## Time consuming, nb = 2

## Imputation
res.comp <- imputePCA(orange,ncp=2)

## A PCA can be performed on the imputed data 
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
res.afdm <- FAMD(ozone,tab.comp=res.comp)



cleanEx()
nameEx("plot.MIPCA")
### * plot.MIPCA

flush(stderr()); flush(stdout())

### Name: plot.MIPCA
### Title: Plot the graphs for the Multiple Imputation in PCA
### Aliases: plot.MIPCA
### Keywords: dplot

### ** Examples

## Not run: 
##D data(orange)
##D ## nb <- estim_ncpPCA(orange,ncp.max=5) ## Time consuming, nb = 2
##D resMI <- MIPCA(orange,ncp=2)
##D plot(resMI)
## End(Not run)



cleanEx()
nameEx("snorena")
### * snorena

flush(stderr()); flush(stdout())

### Name: snorena
### Title: Characterization of people who snore
### Aliases: snorena
### Keywords: datasets

### ** Examples

data(snorena)
res.comp <- imputeFAMD(snorena, ncp=3)
res.afdm <- FAMD(snorena, tab.comp = res.comp)



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
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
