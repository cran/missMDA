\name{imputeFAMD}
\alias{imputeFAMD}
\title{Impute dataset with mixed data}
\description{Impute the missing values of a dataset with the quantitative and categorical variables using the principal component method AFDM. 
Can be used as a preliminary step before performing a AFDM on an incomplete dataset.}

\usage{imputeFAMD(X, ncp = 2, method = "Regularized", row.w = NULL,
      coeff.ridge=1,threshold = 1e-06, seed = NULL, maxiter = 1000,...)
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
  \item{maxiter}{integer, maximum number of iteration for the algorithm}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
Impute the missing entries of a data frame using the iterative FAMD algorithm (EM) or a regularised iterative FAMD algorithm. The iterative FAMD algorithm first imputes the missing values with initial values: the mean of the variable for the continuous variables and the proportion of the category for each category using the non-missing entries. Then performs FAMD on the completed dataset, imputes the missing values with the reconstruction formulae of order ncp and iterates until convergence. The regularized version allows to avoid overfitting problems, especially important when there are many missing values. The output can be used as an input in the AFDM function.
}
\value{
  \item{completeObs}{the imputed dataset; the observed values for non-missing entries and the imputed values for missing values}
  \item{tab.disj}{the imputed matrix in which categorical variables are coding by dummy variables. In the dummy variables, the imputed values are real numbers and may be seen as degree of membership to the corresponding category.}
  \item{call}{the matched call}
}

\author{Francois Husson  \email{husson@agrocampus-ouest.fr} and Julie Josse  \email{Julie.Josse@agrocampus-ouest.fr}
}

\references{Audigier, V., Husson, F. & Josse, J. (2013). A principal components method to impute mixed data.}

\seealso{\code{\link{imputePCA}}}
\examples{
data(ozone)
res.comp <- imputeFAMD(ozone, ncp=3)
res.afdm <- AFDM(ozone,tab.comp=res.comp)
}
\keyword{models}
\keyword{multivariate}
