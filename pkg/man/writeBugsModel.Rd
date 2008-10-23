\name{writeBugsModel}
\alias{writeBugsModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Write a bugs model file for a Generalised Linear Mixed Model}
\description{
  Given a list of effect groups, and the covariates associated with each level, a bugs model file is written using ragged arrays corresponding to output from \code{\link{winBugsRaggedArray}}
}
\usage{
writeBugsModel(file, effects, covariates, observations, family = c("bernoulli", "binomial", "poisson", "normal", "other"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{a character string denoting the name of the bugs model file writen.}
  \item{effects}{vector of effect groups}
  \item{covariates}{A list with names corresponding to effects and each element being a vector of covariates applicable at that level}
  \item{observations}{A character string giving the column of observations, or a vector where the first element is the observations and the remaning are offsets.}
  \item{family}{Response distribution}
}
\details{
  The arguments to the function specify a generalised linear mixed model.  A file containing code for a corresponding bugs model is written.  The model uses ragged arrays to specify grouping factors, and includes covariates at the appropriate levels to aid in chain convergence.  It is intended that the user will edit this file before it's use.  The prior distributions in particular may not be appropriate.
}
\value{
A file, suitable for passing to the \code{\link[R2WinBUGS:bugs]{bugs}} function in \code{R2WinBUGS}.
}
\references{"Handling unbalanced datasets" in the "Tricks: Advanced Use of the BUGS Language" section of the bugs manual, at \code{\url{http://mathstat.helsinki.fi/openbugs/data/Docu/Tricks.html#HandlingUnbalancedDatasets}}}
\author{Patrick Brown, \email{patrick.brown@utoronto.ca}}
\note{}
\section{Warning }{You are strongly encouraged to modify the model file prior to using it. 
}
\seealso{}
\examples{
writeBugsModel("model.bug", effects="Strip",  observations="conc",
  covariates=list(observations="length"),
  family="normal") 
}

