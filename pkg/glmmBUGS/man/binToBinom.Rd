\name{binToBinom}
\alias{binToBinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert Bernoulli observations to Binomial}
\description{
  Combines multiple Bernoulli observations with the same covariates into one Binomial response
}
\usage{
binToBinom(obs, covariates)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{logical vector of observations}
  \item{covariates}{Data frame or matrix of covariates}
}

\value{
  A data frame with one row for each unique value for the covariates, including the covariates and the following additional columns:

  \item{y}{Number of positive observations for the corresponding covariate values}
  \item{N}{Total number of observations for these covariates}
}
 
\author{Patrick Brown}
 
 
\examples{
  thedata = data.frame(sex = rep(c("m", "f"), 10), age=rep(c(20,30), c(10, 10)))
  y = rbinom(dim(thedata)[1], 1, 0.5)
  bindata = binToBinom(y, thedata)
  bindata$zeros = bindata$N - bindata$y
  glm(as.matrix(bindata[,c("y", "zeros")]) ~ sex, data=bindata, family=binomial)
}

