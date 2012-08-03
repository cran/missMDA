\name{imputeMFA}
\alias{imputeMFA}
\title{Impute dataset with MFA}
\description{Impute the missing values of a dataset with the Multiple Factor Analysis model. Can be used as a preliminary step before performing a MFA on an incomplete dataset with coninuous variables.}

\usage{imputeMFA(X, group, ncp = 2, type=rep("s",length(group)), method = "Regularized", 
       row.w = NULL,coeff.ridge=1,threshold = 1e-06, seed = NULL, nb.init = 1,  
	   maxiter = 1000,...)
	   }

\arguments{
  \item{X}{a data.frame with continuous variables containing missing values}
  \item{group}{a list indicating the number of variables in each group}
  \item{ncp}{integer corresponding to the number of components used to reconstruct data with the PCA reconstruction formulae}
  \item{type}{the type of variables in each group; three possibilities: "c" or "s" for quantitative variables
  (the difference is that for "s" variables are scaled to unit variance), "n" for categorical variables; by default,
  all variables are quantitative and scaled to unit variance}
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
If all the groups are quantitative, impute the missing entries of a data frame using the iterative MFA algorithm (EM) or a regularized iterative MFA algorithm. 
The iterative MFA algorithm first imputes the missing values with initial values (the means of each variable), then performs MFA on the completed dataset, imputes the missing values with the reconstruction formulae of order ncp and iterates until convergence. The regularized version allows to avoid overfitting problems, especially important when there are many missing values.\cr
If some groups are qualitative, impute the missing entries of the disjunctive tables for the qualitative groups and the missing entries of the quantitative variables. 
The output can be used as an input in the MFA function using the argument tab.comp.
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
res.comp <- imputeMFA(orange,group=c(5,3),type=rep("s",2),ncp=2)
## Note that MFA is performed on the completed matrix
res.mfa <- MFA(res.comp$completeObs,group=c(5,3),type=rep("s",2))

\dontrun{
data(vnf)
res.comp <- imputeMFA(vnf,group=c(6,5,3),type=c("n","n","n"),ncp=2)
res.mfa <- MFA(vnf,group=c(6,5,3),type=c("n","n","n"),tab.comp=res.comp)
}
}
\keyword{models}
\keyword{multivariate}
