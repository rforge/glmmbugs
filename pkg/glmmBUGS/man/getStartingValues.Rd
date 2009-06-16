\name{getStartingValues}
\alias{getStartingValues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract starting values for an MCMC chain from glmmPQL results}
\description{
  Parameter estimates and random effect predictions are extracted from a glmmPQL model fit, and formatted to correspond to the levels in the supplied ragged array.
}
\usage{
getStartingValues(pql, ragged,prefix=NULL, reparam=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pql}{output from the \code{\link{glmmPQLstrings}} function }
  \item{ragged}{a ragged array, from \code{\link{winBugsRaggedArray}}}
  \item{prefix}{string to append to object names}
  \item{reparam}{vector of random effect names, subtract covariates at this level from the intercept. }
  }
\details{
  This function produces a list suitable for passing to \code{\link{startingFunction}} to generate random starting values for use with \code{\link[R2WinBUGS:bugs]{bugs}}.
  If \code{ragged} has a spatial component, starting values for a spatial random effect will also be computed.
}
\value{
  A list of vectors, one for each set of parameters or random effects, and a list of estimated standard deviations.  
}
\author{Patrick Brown \email{patrick.brown@utoronto.ca}}
\seealso{\code{\link{glmmPQLstrings}}, \code{\link{startingFunction}}, \code{\link[R2WinBUGS:bugs]{bugs}}, \code{\link{glmmBUGS}}
}


