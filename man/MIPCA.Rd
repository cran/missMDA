\name{MIPCA}
\alias{MIPCA}
\title{Multiple Imputation with PCA}
\description{
MIPCA performs Multiple Imputation with a PCA model. Can be used as a preliminary step to perform Multiple Imputation in PCA}
\usage{
MIPCA(X, ncp = 2, scale = TRUE, method = "Regularized", 
   threshold = 1e-04, nboot = 100)
}
\arguments{
  \item{X}{a data.frame with continuous variables containing missing values}
  \item{ncp}{integer corresponding to the number of components used to reconstruct data with the PCA reconstruction formulae}
  \item{scale}{boolean. By default TRUE leading to a same weight for each variable}
  \item{method}{"Regularized" by default or "EM"}
  \item{threshold}{the threshold for the criterion convergence}
  \item{nboot}{the number of imputed datasets}
}
\details{
MIPCA generates nboot imputed datasets from a PCA model. The observed values are the same from one dataset to the others whereas the imputed values change. The variation among the imputed values reflects the variability with which missing values can be predicted.
The multiple imputation is proper in the sense of Little and Rubin (2002) since it takes into account the variability of the parameters.
}
\value{
\item{res.imputePCA}{A matrix corresponding to the imputed dataset obtained with the function imputePCA (the completed dataset)}
\item{res.MI }{An array corresponding to nboot imputed dataset. The dimensions of the array are: the number of row of X, the number of column of X and nboot}
}

\references{Josse, J., Husson, F. (2010). Multiple Imputation in PCA.}

\author{Francois Husson  \email{husson@agrocampus-ouest.fr} and Julie Josse  \email{Julie.Josse@agrocampus-ouest.fr}
}

\seealso{\code{\link{imputePCA}},\code{\link{plot.MIPCA}}}
\examples{
data(orange)
## First the number of components has to be chosen 
##   (for the reconstruction step)
## nb <- estim_ncpPCA(orange,ncp.max=5) ## Time consuming, nb = 2

## Multiple Imputation
resMI <- MIPCA(orange,ncp=2)

## Visualization on the PCA map
plot(resMI)
}

\keyword{multivariate}
