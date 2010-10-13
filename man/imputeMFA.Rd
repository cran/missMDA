\name{imputeMFA}
\alias{imputeMFA}
\title{Impute dataset with MFA}
\description{Impute the missing values of a dataset with the Multiple Factor Analysis model. Can be used as a preliminary step before performing a MFA on an incomplete dataset with coninuous variables.}

\usage{imputeMFA(X, group, ncp = 2, scale = TRUE, method = "Regularized", 
       threshold = 1e-06, seed = NULL, nb.init = 1, maxiter = 1000, ...)
}

\arguments{
  \item{X}{a data.frame with continuous variables containing missing values}
  \item{group}{a list indicating the number of variables in each group}
  \item{ncp}{integer corresponding to the number of components used to reconstruct data with the PCA reconstruction formulae}
  \item{scale}{boolean. By default TRUE leading to a same weight for each variable}
  \item{method}{"Regularized" by default or "EM"}
  \item{threshold}{the threshold for assessing convergence}
  \item{seed}{a single value, interpreted as an integer for the set.seed function (if seed = NULL, missing values are initially imputed by the mean of each variable)}
  \item{nb.init}{integer corresponding to the number of random initializations; the first initialization is the mean of each variable}
  \item{maxiter}{integer, maximum number of iteration for the algorithm}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
Impute the missing entries of a data frame using the iterative MFA algorithm (EM) or a regularized iterative MFA algorithm. The iterative MFA algorithm first imputes the missing values with initial values (the means of each variable), then performs MFA on the completed dataset, imputes the missing values with the reconstruction formulae of order ncp and iterates until convergence. The regularized version allows to avoid overfitting problems, especially important when there are many missing values.
}
\value{
  \item{completeObs}{the imputed dataset; the observed values for non-missing entries and the imputed values for missing values}
  \item{objective}{the criterion, the mean square error between the reconstructed data and the observed data}
  \item{recon}{the reconstructed data}
}

\references{PhD thesis of J. Josse or HDR of F. Husson}

\author{Francois Husson  \email{husson@agrocampus-ouest.fr} and Julie Josse  \email{Julie.Josse@agrocampus-ouest.fr}
}

\seealso{\code{\link{imputePCA}}}
\examples{
data(orange)
res.comp <- imputePCA(orange,group=c(5,3),scale=TRUE,ncp=2)
res.pca <- MFA(res.comp$completeObs,group=c(5,3),type=rep("s",2))
}
\keyword{models}
\keyword{multivariate}
