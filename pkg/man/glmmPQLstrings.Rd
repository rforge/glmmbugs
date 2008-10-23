\name{glmmPQLstrings}
\alias{glmmPQLstrings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{An alternat interface to glmmPQL}
\description{
  Calls glmmPQL in the MASS library, with the model being specified in the same manner as \code{\link{writeBugsModel}}
}
\usage{
glmmPQLstrings(effects, covariates, observations, data = NULL, family=c("bernoulli", "binomial", "poisson", "gaussian"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   
  \item{effects}{A vector of character strings containing the grouping levels, from most general to most specific}
  \item{covariates}{A list with names corresponding to effects and each element being a vector of covariates applicable at that level}
  \item{observations}{A character string giving the column of observations, or a vector where the first element is the observations and the remaning are offsets.  For binomial responses, the first element is the counts (of successes), and the second element is the total number of trials.  Note this differs from glmmPQL and glm's notation, but is consistent with WinBUGS.}
   \item{data}{A data frame containing the response, covariates, and group membership.}
  \item{family}{The distribution to use.  Either using \code{\link[MASS:glmmPQL]{glmmPQL}}'s specifications or \code{\link{writeBugsModel}}}
  \item{\dots}{further arguments to \code{\link[MASS:glmmPQL]{glmmPQL}}}
}
\details{
  This function is useful for generating starting values for an MCMC chain.
}
\value{
In addition to the output from glmmPQL, the following are returned
  \item{effects,covariates,observations}{As input}
}
\references{}
\author{Patrick Brown, \code{patrick.brown@utoronto.ca}}
\note{}
\seealso{\code{\link{getStartingValues}},\code{\link[MASS:glmmPQL]{glmmPQL}}}
\examples{
library(nlme)
data(Muscle)
glmmPQLstrings(effects="Strip",  observations="conc",
  covariates=list(observations="length") ,
  data=Muscle, family="gaussian")
}
