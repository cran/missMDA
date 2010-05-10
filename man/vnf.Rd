\name{vnf}
\alias{vnf}
\docType{data}
\title{
Questionnaire done by 1232 individuals who answered 14 questions
}
\description{
A user satisfaction survey of pleasure craft operators on the ``Canal des Deux Mers'', located in South of France, was carried out by the public corporation ``Voies Navigables de France'' responsible for managing and developing the largest network of navigable waterways in Europe
}
\usage{data(vnf)}
\format{
  A data frame with 1232 observations on the following 14 categorical variables.
}

\source{Josse, J., Chavent, M., Liquet, B. and Husson, F. (2010). Handling missing values with Regularized Iterative Multiple Correspondence Analysis.}

\examples{
data(vnf)
tab.disj <- imputeMCA(vnf,ncp=2)

res.mca <- MCA(vnf,tab.disj=tab.disj)
}
\keyword{datasets}
