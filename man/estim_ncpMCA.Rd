\name{estim_ncpMCA}
\alias{estim_ncpMCA}
\title{Estimate the number of dimensions for the Multiple Correspondence Analysis by cross-validation}
\description{Estimate the number of dimensions for the Multiple Correspondence Analysis by cross-validation}
\usage{
estim_ncpMCA(don, ncp.min=0, ncp.max=5, nbsim=100, pNA=0.05, threshold=1e-4)
}
\arguments{
  \item{don}{a data.frame with categorical variables; with missing entries or not}
  \item{ncp.min}{integer corresponding to the minimum number of components to test}
  \item{ncp.max}{integer corresponding to the maximum number of components to test}
  \item{nbsim}{number of simulations}
  \item{pNA}{percentage of missing values added in the data set}
  \item{threshold}{the threshold for assessing convergence}
}

\details{
For the cross-validation, pNA percentage of missing values are removed at random and predicted with a MCA model using ncp.min to ncp.max dimensions. This process is repeated nbsim times.  The number of components which leads to the smallest MSEP is retained. Each cell is predicted using the imputeMCA function, it means using the regularized iterative MCA algorithm.
}
\value{
  \item{ncp}{the number of components retained for the MCA}
  \item{criterion}{the criterion (the MSEP) calculated for each number of components}
}
\references{Josse, J., Chavent, M., Liquet, B. and Husson, F. (2010). Handling missing values with Regularized Iterative Multiple Correspondence Analysis.}

\author{Francois Husson  \email{husson@agrocampus-ouest.fr} and Julie Josse  \email{Julie.Josse@agrocampus-ouest.fr}
}

\examples{
\dontrun{
data(vnf)
result <- estim_ncpMCA(vnf,ncp.min=0, ncp.max=3, nbsim=100)
}
}
\keyword{multivariate}
