\name{imputeAFDM}
\alias{imputeAFDM}
\title{Impute dataset with mixed data}
\description{Impute the missing values of a dataset with the quantitative and categorical variables using the principal component method AFDM. 
Can be used as a preliminary step before performing a AFDM on an incomplete dataset.}

\usage{imputeAFDM(X, ncp = 2, method = "Regularized", row.w = NULL,
      coeff.ridge=1,threshold = 1e-06, seed = NULL, nb.init = 1, maxiter = 1000,...)
	   }

\arguments{
  \item{X}{a data.frame with continuous and categorical variables containing missing values}
  \item{ncp}{integer corresponding to the number of components used to reconstruct data with the PCA reconstruction formulae}
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
The output can be used as an input in the AFDM function.
}
\value{
  \item{completeObs}{the imputed dataset; the observed values for non-missing entries and the imputed values for missing values}
  \item{objective}{the criterion, the mean square error between the reconstructed data and the observed data}
  \item{recon}{the reconstructed data}
}

\author{Francois Husson  \email{husson@agrocampus-ouest.fr} and Julie Josse  \email{Julie.Josse@agrocampus-ouest.fr}
}

\seealso{\code{\link{imputePCA}}}

\keyword{models}
\keyword{multivariate}
