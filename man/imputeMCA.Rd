\name{imputeMCA}
\alias{imputeMCA}

\title{Impute missing values in categorical variables with Multiple Correspondence Analysis}
\description{Impute the missing values of a categorical dataset (in the indicator matrix) with Multiple Correspondence Analysis}
\usage{
imputeMCA(don, ncp=2, threshold=1e-06, seed=NULL, maxiter=1000)
}

\arguments{
  \item{don}{a data.frame with categorical variables containing missing values}
  \item{ncp}{integer corresponding to the number of dimensions used to reconstruct data with the reconstruction formulae}
  \item{threshold}{the threshold for assessing convergence}
  \item{seed}{an integer to specify the seed for the initialization for the regularized iterative MCA algorithm (if seed = NULL the initialization step corresponds to the imputation of the proportion of each category)}
  \item{maxiter}{integer, maximum number of iterations for the regularized iterative MCA algorithm}
}
\details{
Use a Regularized Iterative Multiple Correspondence Analysis to impute missing values. 
The regularized iterative MCA algorithm first imputes the missing values in the indicator matrix with initial values (the proportion of each category), then performs MCA on the completed dataset, imputes the missing values with the reconstruction formulae of order ncp and iterates until convergence. 

If ncp=0, the Average method (imputation with the proportion) is performed.
}

\value{
Return the imputed indicator matrix. The imputed valued are real numbers and may be seen as degree of membership to the corresponding category.
}

\author{
Francois Husson  \email{husson@agrocampus-ouest.fr} and Julie Josse  \email{Julie.Josse@agrocampus-ouest.fr}
}

\references{Josse, J., Chavent, M., Liquet, B. and Husson, F. (2010). Handling missing values with Regularized Iterative Multiple Correspondence Analysis.}

\seealso{\code{\link{estim_ncpMCA}}}

\examples{
\dontrun{
data(vnf)
## First the number of components has to be chosen 
##   (for the reconstruction step)
## nb <- estim_ncpMCA(vnf,ncp.max=5) ## Time-consuming, nb = 4

## Impute indicator matrix
tab.disj <- imputeMCA(vnf, ncp=4)

## A MCA can be performed
res.mca <- MCA(vnf,tab.disj=tab.disj)
}
}

\keyword{models}
\keyword{multivariate}
