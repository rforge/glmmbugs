\name{getDesignMatrix}
\alias{getDesignMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes a design matrix from factors and interactions}
\description{
  Converts all factors and interactions to indicator variables, suitable for passing to WinBUGS.  
}
\usage{
getDesignMatrix(formula, data, effects = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A \code{formula} object specifying the fixed effects for the model}
  \item{data}{A data frame containing the covariates and factors for random effects}
\item{effects}{A vector of character strings containing the grouping levels, from most general to most specific}
}
\details{
    The most populous level of a factor is made the baseline.
}
\value{
A matrix containing the covariates, the response(s), and the random effect factors.  Also attributes
  \item{covariates}{A list giving the covariates which apply at each level, suitable for passing to \code{\link{winBugsRaggedArray}}}
  \item{response}{A vector of character strings giving the responses}
}
\author{Patrick Brown}
\seealso{\code{\link{winBugsRaggedArray}}, \code{\link{glmmBUGS}}}
\examples{
library(nlme)
data(Muscle)
muscleDesign = getDesignMatrix(conc ~ length, data=Muscle, effects="Strip" )
attributes(muscleDesign)$covariates
attributes(muscleDesign)$response
}

