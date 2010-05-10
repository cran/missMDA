\name{estim_ncpPCA}
\alias{estim_ncpPCA}
\title{Estimate the number of dimensions for the Principal Component Analysis by cross-validation}
\description{Estimate the number of dimensions for the Principal Component Analysis by cross-validation}
\usage{
estim_ncpPCA(X, ncp.min = 0, ncp.max = 5, method = "Regularized", scale = TRUE, method.cv = "loo", nbsim = 100, pNA = 0.05, threshold=1e-4)
}
\arguments{
  \item{X}{a data.frame with continuous variables; with missing entries or not}
  \item{ncp.min}{integer corresponding to the minimum number of components to test}
  \item{ncp.max}{integer corresponding to the maximum number of components to test}
  \item{method}{"Regularized" by default or "EM"}
  \item{scale}{boolean. By default TRUE leading to a same weight for each variable}
  \item{method.cv}{string with the values "loo" for leave-one-out or "Kfold" cross-validation}
  \item{nbsim}{number of simulations, useful only if method.cv="Kfold"}
  \item{pNA}{percentage of missing values added in the data set, useful only if method.cv="Kfold"}
  \item{threshold}{the threshold for assessing convergence}
}

\details{
For leave-one-out (loo) cross-validation, each value is alternatively removed and predicted with a PCA model using ncp.min to ncp.max dimensions. The number of components which leads to the smallest MSEP is retained. Each cell is predicted using the imputePCA function, it means using the regularized iterative PCA algorithm or the iterative PCA (EM cross-validation).\cr
For the Kfold cross-validation, pNA percentage of missing values is removed and predicted with a PCA model using ncp.min to ncp.max dimensions. This process is repeated nbsim times.
The leave-one-out method is time-consuming (method.cv="loo") when the number of cells is important in the data.frame.\cr
The regularized version is more appropriate when there are many missing values in the dataset (to avoid overfitting).
}
\value{
  \item{ncp}{the number of components retained for the PCA}
  \item{criterion}{the criterion (the MSEP) calculated for each number of components}
}
\references{Bro, R., Kjeldahl, K. Smilde, A. K. and Kiers, H. A. L. (2008) Cross-validation of component models: A critical look at current methods. 	Analytical and Bioanalytical Chemistry, 5, 1241-1251.\cr
J. Josse, F. Husson et J. Pagès (2009) Gestion des données manquantes en Analyse en Composantes Principales. Journal de la SFdS. 150 (2), pp. 28-51.}

\author{Francois Husson  \email{husson@agrocampus-ouest.fr} and Julie Josse  \email{Julie.Josse@agrocampus-ouest.fr}
}

\examples{
\dontrun{
data(orange)
nb <- estim_ncpPCA(orange,ncp.min=0,ncp.max=4) ## Time consuming, nb = 2
}
}
\keyword{multivariate}
