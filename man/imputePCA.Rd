\name{imputePCA}
\alias{imputePCA}
\title{Impute dataset with PCA}
\description{Impute the missing values of a dataset with the Principal Components Analysis model. Can be used as a preliminary step before performing a PCA on an incomplete dataset.}
\usage{imputePCA(X, ncp = 2, scale = TRUE, method = "Regularized", 
       row.w=NULL, coeff.ridge=1, threshold = 1e-06, seed = NULL, nb.init = 1,  
	   maxiter = 1000, ...)
}

\arguments{
  \item{X}{a data.frame with continuous variables containing missing values}
  \item{ncp}{integer corresponding to the number of components used to reconstruct data with the PCA reconstruction formulae}
  \item{scale}{boolean. By default TRUE leading to a same weight for each variable}
  \item{method}{"Regularized" by default or "EM"}
  \item{row.w}{an optional row weights (by default, a vector of 1 over the number of rows for uniform row weights)}
  \item{coeff.ridge}{a positive coefficient that permits to shrink the eigenvalues more than by the mean of the last eigenvalues 
  (by default, 1 the eigenvalues are shrunk by the mean of the last eigenvalues; a coefficient between 1 and 2 is required)}
  \item{threshold}{the threshold for assessing convergence}
  \item{seed}{a single value, interpreted as an integer for the set.seed function (if seed = NULL, missing values are initially imputed by the mean of each variable)}
  \item{nb.init}{integer corresponding to the number of random initializations; the first initialization is the mean of each variable}
  \item{maxiter}{integer, maximum number of iteration for the algorithm}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
Impute the missing entries of a data frame using the iterative PCA algorithm (EM) or a regularized iterative PCA algorithm. The iterative PCA algorithm first imputes the missing values with initial values (the means of each variable), then performs PCA on the completed dataset, imputes the missing values with the reconstruction formulae of order ncp and iterates until convergence. The regularized version allows to avoid overfitting problems, especially important when there are many missing values.
}
\value{
  \item{completeObs}{the imputed dataset; the observed values for non-missing entries and the imputed values for missing values}
  \item{objective}{the criterion, the mean square error between the reconstructed data and the observed data}
  \item{recon}{the reconstructed data}
}

\references{J. Josse, F. Husson et J. Pagès (2009) Gestion des données manquantes en Analyse en Composantes Principales. Journal de la SFdS. 150 (2), pp. 28-51.\cr
Josse, J., Husson, F. (2010). Multiple Imputation in PCA.}

\author{Francois Husson  \email{husson@agrocampus-ouest.fr} and Julie Josse  \email{Julie.Josse@agrocampus-ouest.fr}
}

\seealso{\code{\link{estim_ncpPCA}},\code{\link{MIPCA}}}
\examples{
data(orange)
## First the number of components has to be chosen 
##   (for the reconstruction step)
## nb <- estim_ncpPCA(orange,ncp.max=5) ## Time consuming, nb = 2

## Imputation
res.comp <- imputePCA(orange,ncp=2)

## A PCA can be performed
res.pca <- PCA(res.comp$completeObs)
}
\keyword{models}
\keyword{multivariate}
